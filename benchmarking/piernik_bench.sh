#!/bin/bash

export LC_ALL=C

# The defaults
DO_MAKE=1  # measure the compilation time scaling
COMPILER_CONFIG="benchmarking"  # for the setup script it means "./compilers/../benchmarking/${COMPILER_CONFIG}.in"
N_PROC_LIST=""

usage() {
    echo "Usage:   $0 [options] [n_threads_1 [n_threads_2 ...]]"
    echo "         -h | --help : this message"
    echo "         -f | --fast : skip compilation test if possible"
    echo "         -c | --config file : use ./benchmarking/\${file}.in as the compiler config (default: ./benchmarking/benchmarking.in)"
    echo "default: $0 \$( seq number_of_logical_CPUs )"
    exit 1
}

parse_arguments() {
    while (( "$#" )); do
        case "$1" in
            -f|--fast)
                DO_MAKE=0
                shift
                ;;
            -h|--help)
                usage
                ;;
            -c|--config)
                if [ -n "$2" ] && [ "${2:0:1}" != "-" ]; then
                    COMPILER_CONFIG="$2"
                    shift 2
                else
                    echo "Error: Argument for $1 is missing" >&2
                    exit 1
                fi
                ;;
            -*|--*=)
                echo "Error: Unsupported flag $1" >&2
                exit 1
                ;;
            *)
                N_PROC_LIST="${N_PROC_LIST} $1"
                shift
                ;;
        esac
    done
}

check_dependencies() {
    command -v make > /dev/null 2>&1 || { echo >&2 "make is required but it's not installed. Aborting."; exit 1; }
    command -v mpirun > /dev/null 2>&1 || { echo >&2 "mpirun is required but it's not installed. Aborting."; exit 1; }
}

print_system_info() {
    MEMG=$( LC_ALL=C free -m | awk '/Mem/ {print $2/1024.}' )
    echo "## "$( grep 'model name' /proc/cpuinfo | uniq )
    echo "## Memory : $MEMG GB"
}

prepare_directories() {
    touch Makefile  # run correctly also when the clock is messed up
    for problem in "${PROBLEM_LIST[@]}"; do
        rm -rf runs/${problem}_B_${problem}
    done
    rm -rf obj_B_*
}

awkfor() {
    awk '{printf("%19s %8s %8s %8s %8s\n", $1, $2, $3, $4, $5)}'
}

make_objects() {
    local obj_list=$1
    make $obj_list CL=1 > /dev/null
    ( time make -j $obj_list > /dev/null ) 2>&1 | awkfor
}

setup_flood_scaling() {
    local threads=$1
    for j in $( seq $threads ); do
        [ -e $j ] && rm -rf $j
        mkdir $j
        cp piernik problem.par $j
    done
}

run_strong_weak_scaling() {
    local scaling=$1
    local threads=$2
    local nx=$3
    local mpirun_cmd=$4
    local max_mem=$5
    local xmul=$6

    case $scaling in
        weak)
            $mpirun_cmd -np $threads ./piernik -n '&BASE_DOMAIN n_d = '$(( $threads * $nx ))', 2*'$nx' xmin = -'$(( $threads * $xmul ))' xmax = '$(( $threads * $xmul ))' / &MEMORY max_mem = '$max_mem'/' 2> /dev/null
            ;;
        strong)
            $mpirun_cmd -np $threads ./piernik -n '&BASE_DOMAIN n_d = 3*'$nx' / &MEMORY max_mem = '$max_mem'/' 2> /dev/null
            ;;
    esac
}

run_piernik() {
    local problem=$1
    local scaling=$2
    local threads=$3
    local mpirun_cmd=$4
    local max_mem=$5

    case $problem in
        sedov)
            local basesize=64 ;;
        crtest)
            local basesize=32 ;;
        maclaurin)
            if [ $scaling == strong ] ; then
                local basesize=128
            else
                local basesize=64
            fi ;;
    esac
    local nx=$( echo $basesize $SCALE | awk '{print int($1*$2)}' )

    if [ $scaling == flood ]; then
        for j in $( seq $threads ); do
            ( cd $j && ./piernik -n '&BASE_DOMAIN n_d = 3*'$nx' / &MEMORY max_mem = '$max_mem'/' > _stdout_ 2> /dev/null ) &
        done
        wait
        sleep 1
        for j in $( seq $threads ); do
            process_output $problem $threads $j
        done
    else
        case $problem in
            sedov)
                local xmul=1
                run_strong_weak_scaling $scaling $threads $nx "$mpirun_cmd" $max_mem $xmul | grep "dWallClock" | awk 'BEGIN {t=0; n=0;} {if ($12 != 0.) {printf(" %8.4f ", $12); t+=$12; n++;}} END {printf(" %9.5f ", t/n)}'
                ;;
            crtest)
                local xmul=512
                run_strong_weak_scaling $scaling $threads $nx "$mpirun_cmd" $max_mem $xmul | grep "p+-cycles" | awk '{if (NR==1) printf("%8.4f %8.4f ", $5, $8)}'
                awk '/Spent/ { printf("%s ", $5) }' *log
                ;;
            maclaurin)
                local xmul=2
                run_strong_weak_scaling $scaling $threads $nx "$mpirun_cmd" $max_mem $xmul | grep "cycles" | awk '{printf(" %8.4f %8.4f ", $5, $8)}'
                awk '/Spent/ { printf("%s ", $5) }' *log
                grep -q Spent *log || return 1  # exception: some piernik threads have returned prematurely
                ;;
        esac
    fi
    return 0
}

format_output() {
    awk '{
        if ($0 ~ "##") {
            print "## ", $0
        } else {
            if (substr($1, 0, 1) == "#") split($0, form)
            if (NF > 0) {
                for (i = 1; i <= NF; i++) printf("%-*s ", length(form[i]), $i)
                print ""
                fflush()
            }
        }
    }'
}

run_benchmark() {
    local problem=$1
    local scaling=$2
    local rundir="runs/${problem}_B_${problem}"
    local config_file="$(dirname $0)/problem.par.${problem}"

    echo "## "$(date)
    echo "Benchmarking ${problem}, ${scaling} scaling"

    (
        cp $config_file $rundir/problem.par
        cd $rundir || exit

        print_header $problem

        for threads in $N_PROC_LIST; do
            local mpirun_cmd="mpirun"
            mpirun -np $threads echo > /dev/null 2>&1 || mpirun_cmd="mpirun --use-hwthread-cpus"
	    # OpenMPI refuses to run more jobs than CPU cores but with --use-hwthread-cpus its performance is poor on 2 and more threads
            local max_mem=$( echo $MEMM $threads | awk '{print int(0.90*$1*1024/$2)}' )

            rm *log 2> /dev/null

            if [ $scaling == flood ]; then
                setup_flood_scaling $threads  # refresh the subdirectories for flood scaling runs
            else
                echo -n $threads
            fi

            run_piernik $problem $scaling $threads "$mpirun_cmd" $max_mem
            OOM=$?
            [ $scaling != flood ] && echo
            echo

            # Check if we got any exception
            if [ $OOM != 0 ]; then
                echo "## Warning: Invalid output detected (possible OOM). Skipping higher thread counts."
                break
            fi
        done
    ) | format_output
    echo
}

print_header() {
    local problem=$1
    case $problem in
        sedov)     echo "#Threads dWallClock1 dWallClock2 dWallClock3 dWallClock4 dWallClock5 dWallClock_Average";;
        crtest)    echo "#Threads MG_prepare MG_cycle Total_MG";;
        maclaurin) echo "#Threads MG_prepare MG_i-cycle MG_multipole MG_o-cycle Total_MG";;
    esac
}

process_output() {
    local problem=$1
    local core=$2
    local dir=${3:-.}
    case $problem in
        sedov)
            ( grep "dWallClock" $dir/_stdout_ || echo "" ) | awk 'BEGIN {t=0; n=0; printf("%3d",'$core');} {if ($3 != 0) {printf(" %8.4f ", $12); t+=$12; n++;}} END {printf(" %9.5f\n", t/n)}'
            ;;
        crtest)
            grep "p+-cycles" $dir/_stdout_ | awk '{if (NR==1) printf("%d %8.4f %8.4f ", '$core', $5, $8)}'
            awk '/Spent/ { printf("%s\n", $5) }' $dir/*log
            ;;
        maclaurin)
            grep cycles $dir/_stdout_ | awk 'BEGIN {printf("%d", '$core');} {printf(" %8.4f %8.4f ", $5, $8)}'
            grep -q cycles $dir/_stdout_ || echo ""
            awk '/Spent/ { printf("%s\n", $5) }' $dir/*log
            ;;
    esac
}

run_on_cores() {
    local problem=$1
    local nx=$2
    local cores=$3
    local sequential=${4:-0}
    local header=$5

    setup_flood_scaling "$N_CPU"

    echo "## "$(date)
    echo $header
    ( print_header $problem
    for j in $cores; do
        if [ $sequential -eq 1 ]; then
            ( cd $j && taskset -c $(( $j - 1 )) ./piernik -n '&BASE_DOMAIN n_d = 3*'$nx' /' > _stdout_ 2> /dev/null
              process_output $problem $j
              sleep 1 )
        else
            ( cd $j && taskset -c $(( $j - 1 )) ./piernik -n '&BASE_DOMAIN n_d = 3*'$nx' /' > _stdout_ 2> /dev/null ) &
        fi
    done
    [ $sequential -eq 0 ] && wait && sleep 1
    [ $sequential -eq 0 ] && for j in $cores; do
        process_output $problem $j $j
    done ) | format_output
    echo
}

profile_cores() {
    local problem=$1
    local rundir="runs/${problem}_B_${problem}"
    local config_file="$(dirname $0)/problem.par.${problem}"

    case $problem in
        sedov)     local basesize=64 ;;
        crtest)    local basesize=32 ;;
        maclaurin) local basesize=64 ;;
    esac

    local nx=$( echo $basesize $SCALE | awk '{print int($1*$2)}' )

    (
        cp $config_file $rundir/problem.par
        cd $rundir || exit

        # Run on each core sequentially
        run_on_cores $problem $nx "$( seq $N_CPU )" 1 "## Core profiling on ${problem} problem, single core on cores 1 to $N_CPU"

        # Run on physical cores in parallel
        # Trick: +1 to "processor" number to be consistent with the 1-based numbering elsewhere
        PHYS_CORES=$(awk '
            function min(a,b) { return (a<b)?a:b; }
            # Process lines containing "processor", "core id", or "physical id"
            /^processor|^core id|^physical id/ {
                if ($1 ~ /processor/) p=$NF;
                if ($1 ~ /physical/) phys=$NF;
                if ($1 ~ /core/) {
                    core_id = phys "-" $NF;  # Combine physical id and core id
                    if (core_id in cores) cores[core_id] = min(p, cores[core_id]);
                    else cores[core_id] = p;
                }
            }
            # Print the physical cores
            END {
                for (cpu in cores) printf "%d ", cores[cpu] + 1;
            }
        ' /proc/cpuinfo | sed 's/ $/\n/')
        run_on_cores $problem $nx "$PHYS_CORES" 0 "## Core profiling on ${problem} problem, all physical cores ($PHYS_CORES)"

        # Run on all threads in parallel
        run_on_cores $problem $nx "$( seq $N_CPU )" 0 "## Core profiling on ${problem} problem, all threads (1 to $N_CPU)"
    )
}

# Main script execution
parse_arguments "$@"
check_dependencies
print_system_info

SCALE=${BIG:=1}

MEMM=$( LC_ALL=C free -m | awk '/Mem/ {print $2}' )
# alternatively (should give value closer to the amount of physical RAM installed):
#MEMG=0
#for mem in /sys/devices/system/memory/memory*; do
#    [[ "$(cat ${mem}/online)" == "1" ]] && MEMG=$(( MEMG+$((0x$(cat /sys/devices/system/memory/block_size_bytes)))))
#done
#MEMG=$( echo $MEMG | awk '{print $1 / 1024.^3 }' )

[ "$SCALE" != "1" ] && echo "# test domains are scaled by factor of $SCALE"

PROBLEM_LIST=("maclaurin" "advection_test" "crtest" "jeans" "tearing" "sedov" "otvortex" "2body")
B_PROBLEM_LIST=("sedov" "crtest" "maclaurin")

# create list of thread count to be tested
N_CPU=$( awk 'BEGIN {c=0} /processor/ {if ($NF > c) c=$NF} END {print c+1}' /proc/cpuinfo )
if [ -z "$N_PROC_LIST" ]; then
    N_PROC_LIST=$( seq $N_CPU )
else
    N_PROC_LIST=$( echo $N_PROC_LIST | awk '{for (i=1; i<=NF; i++) printf("%d\n",1*$i); print ""}' | sort -n | uniq )
fi

for i in $N_PROC_LIST; do
    if [ $i -le 0 ]; then
        echo "N_PROC_LIST: ${N_PROC_LIST}"
        usage
    fi
done

TIMEFORMAT='real/user/sys/CPU%%: %1R %1U %1S %P%%'

if [ $DO_MAKE == 0 ] ; then
    for p in "${B_PROBLEM_LIST[@]}"; do
        [ -x runs/${p}_B_${p}/piernik ] || DO_MAKE=1
    done
fi

if [ $DO_MAKE == 1 ]; then
    prepare_directories
    {
        echo -n "Preparing objects                "
        SETUP_PARAMS="-c ../benchmarking/${COMPILER_CONFIG} -n --linkexe -d BENCHMARKING_HACK"
        ( time for problem in "${PROBLEM_LIST[@]}"; do
            ./setup $problem $SETUP_PARAMS -o "B_$problem" > /dev/null
        done ) 2>&1 | awkfor

        get_n_problems() {
            echo $1 "${PROBLEM_LIST[@]}" | awk '{for (i=0; i<$1; i++) printf("obj_B_%s ",$(2+i)); print ""}'
        }

        echo -n "Single-thread make object        "
        ( time make $( get_n_problems 1 ) > /dev/null ) 2>&1 | awkfor

        echo -n "Multi-thread make object         "
        make_objects "$(get_n_problems 1)"

        echo -n "Multi-thread make two objects    "
        make_objects "$(get_n_problems 2)"

        echo -n "Multi-thread make four objects   "
        make_objects "$(get_n_problems 4)"

        echo -n "Multi-thread make eight objects  "
        make_objects "$(get_n_problems 8)"
    }
else
    echo "## compilation skipped"
fi
echo

# for problem in "${B_PROBLEM_LIST[@]}"; do
for problem in crtest ; do
    profile_cores $problem
done

for problem in "${B_PROBLEM_LIST[@]}"; do
    for scaling in flood strong weak; do
        run_benchmark $problem $scaling
    done
done

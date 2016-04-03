#!/bin/bash

# create list of thread count to be tested
if [ $# -lt 1 ] ; then
    N=$( awk 'BEGIN {c=0} /processor/ {if ($NF > c) c=$NF} END {print c+1}' /proc/cpuinfo )
    N_PROC_LIST=$( seq $N ) 
else
    N_PROC_LIST=$( echo $* | awk '{for (i=1; i<=NF; i++) printf("%d ",1*$i); print ""}' )
fi

for i in $N_PROC_LIST ; do
    if [ $i -le 0 ] ; then
	echo "Usage:   $0 [n_threads_1 [n_threads_2 ...]]"
	echo "default: $0 \$( seq number_of_logical_CPUs )"
	exit 1
    fi
done

TIMEFORMAT='real/user/sys/CPU%%: %1R %1U %1S %P%%'
MAKE=make

# some cleanup
PROBLEM_LIST="maclaurin advection_test crtest jeans tearing sedov otvortex 3body"
rm -rf obj_B_*
for i in $PROBLEM_LIST; do
    rm -rf runs/${i}_B_${i}
done

{
# create object and run directories
echo -n "Preparing objects "
SETUP_PARAMS="-c ../benchmarking/benchmarking -n"
time for i in $PROBLEM_LIST; do
    ./setup $i $SETUP_PARAMS -o "B_"$i > /dev/null
done

get_n_problems() {
    echo $1 $PROBLEM_LIST | awk '{for (i=0; i<$1; i++) printf("obj_B_%s ",$(2+i)); print ""}'
}

#
# Benchmarking: make objects
#

OBJ_LIST=$( get_n_problems 1 )

echo -n "Single-thread make object "
time $MAKE $OBJ_LIST > /dev/null
$MAKE $OBJ_LIST CL=1 > /dev/null

echo -n "Multi-thread make object "
time $MAKE -j $OBJ_LIST > /dev/null
$MAKE $OBJ_LIST CL=1 > /dev/null

OBJ_LIST=$( get_n_problems 2 )

echo -n "Multi-thread make two objects "
time $MAKE -j $OBJ_LIST > /dev/null
$MAKE $OBJ_LIST CL=1 > /dev/null

OBJ_LIST=$( get_n_problems 4 )

echo -n "Multi-thread make four objects "
time $MAKE -j $OBJ_LIST > /dev/null
$MAKE $OBJ_LIST CL=1 > /dev/null

OBJ_LIST=$( get_n_problems 8 )

echo -n "Multi-thread make eight objects "
time $MAKE -j $OBJ_LIST > /dev/null
} 2>&1 | $( dirname $0 )"/pretty_time_form.awk"
echo

#
# Benchmarking: Piernik
#

for p in sedov crtest maclaurin ; do
    for t in weak strong ; do
	echo "Benchmarking $p, $t scaling"
	(
	    RUNDIR=runs/${p}_B_${p}
	    cp  $( dirname $0 )"/problem.par."${p} ${RUNDIR}/problem.par
	    cd $RUNDIR
	    case $p in
		sedov)     echo "#Threads dWallClock1 dWallClock2 dWallClock3 dWallClock4 dWallClock5 dWallClock_Average";;
		crtest)    echo "#Threads MG_prepare MG_cycle Total_MG";;
		maclaurin) echo "#Threads MG_prepare MG_i-cycle MG_multipole MG_o-cycle Total_MG";;
	    esac
	    for i in $N_PROC_LIST ; do
		rm *log 2> /dev/null
		echo -n $i
		case $p in
		    sedov)
			case $t in
			    weak)
				mpirun -np $i ./piernik -n '&BASE_DOMAIN n_d = '$(( $i * 64 ))', 2*64 xmin = -'$(( $i * 1 ))' xmax = '$(( $i * 1 ))'/' 2> /dev/null ;;
			    strong)
				mpirun -np $i ./piernik -n '&BASE_DOMAIN n_d = 3*64 /' 2> /dev/null ;;
			esac | grep "dWallClock" | awk 'BEGIN {t=0; n=0;} {if ($12 != 0.) {printf("%7.2f ", $12); t+=$12; n++;} } END {printf("%7.3f\n", t/n)}' ;;
		    crtest)
			case $t in
			    weak)
				mpirun -np $i ./piernik -n '&BASE_DOMAIN n_d = '$(( $i * 32 ))', 2*32 xmin = -'$(( $i * 512 ))' xmax = '$(( $i * 512 ))'/' 2> /dev/null ;;
			    strong)
				mpirun -np $i ./piernik -n '&BASE_DOMAIN n_d = 3*32 /' 2> /dev/null ;;
			esac | grep "C1-cycles" | awk '{if (NR==1) printf("%7.3f %7.3f ", $5, $8)}'
			awk '/Spent/ { print $5 }' *log ;;
		    maclaurin)
			case $t in
			    weak)
				mpirun -np $i ./piernik -n '&BASE_DOMAIN n_d = '$(( $i * 64 ))', 2*64 xmin = -'$(( $i * 2 ))' xmax = '$(( $i * 2 ))' / &MPI_BLOCKS AMR_bsize = 3*32 /' 2> /dev/null ;;
			    strong)
				mpirun -np $i ./piernik -n '&BASE_DOMAIN n_d = 3*128 / &MPI_BLOCKS AMR_bsize = 3*32 /' 2> /dev/null ;;
			esac | grep cycles | awk '{printf("%7.3f %7.3f ", $5, $8)}'
			awk '/Spent/ { print $5 }' *log ;;
	        esac
	    done
	) | column -t
	echo
    done
done

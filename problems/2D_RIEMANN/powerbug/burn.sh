#!/bin/bash

# Sample script that may help to stress-test a computer with Piernik.  This
# variant is intended for stress-testing RAM problems.  For stressing
# PSU/VRM/cooling, see the other variant in problems/jeans/powerbug .

# By using the right size of the 2D Riemann problem and smart decomposition
# into blocks, this setup may test nearly whole RAM while reaching quite high
# usage of CPU power.

# Tested only on Linux. Will require changes to run on MacOS.

MODE="ramspam"
if [ $# -ge 1 ] ; then
    MODE=$1
    if [ "$MODE" == "help" ] ; then
	echo "Usage: $0 [mode]"
	echo "  where the mode is:"
	echo "  burn    – for huge CPU heating (may need adjustments for particular CPU)"
	echo "  ramspam – for high CPU heating and a lot of RAM usage (default)"
	echo "  scan    – tries few combinations of adjustable parameters to help to determine optimum block size for the burn mode"
	echo "  user    – user-defined: provide bs and k as arguments"
	echo "  help    – this help"
	exit 0
    fi
fi

# Count available CPU threads.
NTHR=$( lscpu -p | grep -cv '^#' )
VENDOR=$( grep vendor_id /proc/cpuinfo | uniq | awk '{print $3}' )
k=1
#SL="Riemann"
SL="RTVD"
case $MODE in
    "burn")
	# Maximum power draw depends on CPU model and problem configuration:
	# Ryzen 7 2700: bs=300 solver=RTVD
	# i5-11600: bs=116 solver=RTVD
	case $VENDOR in
	    "GenuineIntel")
		BSL=116
		;;
	    "AuthenticAMD")
		BSL=300
		;;
	    *)
		echo "Unknown vendor '$VENDOR'. Do a scan first, then make an update."
		exit 1
		;;
	esac
	;;
    "user")
	# Get $bs and $k from arguments. Bogus values will crash Piernik thus no error checking	here.
	if [ $# -ge 2 ] ; then
	    BSL=$2
	else
	    echo "You need to provide at least bs, e.g.:"
	    echo "$0 user 300 1  # equivalent to '$0 burn' on AMD"
	    exit 5
	fi
	[ $# -ge 3 ] && k=$3
	;;
    "ramspam")
	BSL=256
	# Use all the available memory.
	# Assume an estimate for RES = 19450 + 11960 * k**2 kB
	k=$( LC_ALL=C free | awk '/Mem/ {print int(sqrt(($NF/'"$NTHR"'-19450.)/11960.))}' )
	if [ "$k" -lt 1 ] ; then
	    echo "k=$k means that there is not enough free RAM"
	    exit 6
	fi
	;;
    "scan")
	# Scan for most energy-consuming setup parameters.
	# It seems that n_d has the biggest impact. Its optimal size may change with size and organization of the CPU cache.

	BSL="92 116 148 188 236 300 380 480"
	SL="Riemann RTVD"
	;;
    *)
	echo "Unknown mode '$MODE'."
	exit 2
	;;
esac

DURATION=604800
[ "$MODE" == "scan" ] && DURATION=30

echo "Mode=$MODE Vendor=$VENDOR Nthr=$NTHR bs=$BSL solver=$SL k=$k duration=${DURATION}s"

# If there is mprime  running for warm-up, stop it.
MPRIME=1
killall -STOP mprime 2> /dev/null || MPRIME=0

SHOW_TEMP=.__st__
# A "daemon" for reporting temperatures. Adjust the grep filter to particular machine.
case $VENDOR in
    "GenuineIntel")
	SENS='(Package)'
	;;
    "AuthenticAMD")
	SENS='(Tctl|TSI0_TEMP)'
	;;
    *)
	echo "Unknown vendor '$VENDOR'. Figure out sensors first, then make an update."
	exit 3
	;;
esac
[ $( sensors | grep -E $SENS | wc -l ) == 0 ] && SENS='(:.*°C)'  # Fallback
touch $SHOW_TEMP
sleep 5 && while [ -e $SHOW_TEMP ] ; do
    ( cat $SHOW_TEMP; sensors | grep -E $SENS | sed 's/\([^°]*\)°C.*/\1°C/;s/  */ /' ) | tr '\n' ' ' ; echo
    sleep 5
    killall -CONT piernik 2> /dev/null || rm $SHOW_TEMP
done &

SHOW_STEP=.__step__
touch $SHOW_STEP
# Check if all instances are reaching the same values.
OSTEP=0
sleep 5 && while [ -e $SHOW_TEMP ] ; do
    [ -e 1/out ] && \
	STEP=$( tail -n 1 [1-9]*/out | \
		    grep nstep | \
		    sort -n -k 3 | \
		    head -n 1 | \
		    awk '{print $3}' )
    [ "${STEP:-0}" -gt 0 ] && \
        if [ "$STEP" != "$OSTEP" ] ; then
            OSTEP=$STEP
            grep "nstep = *$STEP" [1-9]*/out | \
		sort -n | \
		column -t | \
		awk '{if (NR==1) {dt=$9; t=$12; print} else if (dt != $9 || t != $12) print "Differs: ", $0}'
        fi
    sleep 10
done &

# Create separate directories for each run.
for i in $( seq "$NTHR" ) ; do
    [ -d "$i" ] && rm -r "$i"
    mkdir "$i"
    cp problem.par "$i"
    ln -s ../piernik "${i}/piernik"
done

for s in $SL ; do
    for bs in $BSL ; do
	nd=$(( "$bs * $k" ))
	echo "${MODE}: n_d= $nd  solver= $s" > $SHOW_TEMP
	# Run $NTHR independent single-threaded Piernik copies as this draws most electrical power.
	# Parallel Piernik runs do spend some time on communication, so we choose flooding until we implement shared-memory :)
	for d in $( seq "$NTHR" ) ; do
	    (
		cd "$d" || exit 4
		rm ./*log 2> /dev/null
		./piernik -n '&BASE_DOMAIN n_d = 2*'"$nd"', 1/ &AMR bsize = 2*'"$bs"', 1/ &NUMERICAL_SETUP solver_str="'"$s"'"' | grep --line-buffered . > out &
	    )
	done

        # ToDo also check if all $NTHR are alive and not stuck.

	# Run each Piernik configuration for limited amount of time.
	sleep $DURATION && killall piernik
	echo "kill $nd $k"
	sleep 1 && killall -9 piernik
	echo
    done
done

rm $SHOW_TEMP $SHOW_STEP

# Continue mprime, if it was there.
[ $MPRIME == 1 ]  && killall -CONT mprime

# This may help to visualize CPU temperature changes.
# gnuplot
# se sty dat lp
# p "< grep Tctl outfile" u 7
# while (1) {rep; pause 3}

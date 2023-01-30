#!/bin/bash

# Sample script that may help to stress-test a computer with Piernik.  This
# variant is intended for stress-testing power supply (PSU, VRM, connectors)
# an thermal (cooling) problems.  For stressing RAM, see the other variant
# in problems/2D_RIEMANN/powerbug .

# By using the right size of the 2D Jeans problem, in a configuration
# limited to red-black Gauss-Seidel relaxation on a single level, this setup
# may reach extremely high usage of CPU resources.  On a Ryzen 7 an average
# of 1.95 IPC was reached, whereas mprime 28.5 reached only 1.87 IPC.  This
# resulted in about 10% higher CPU power demand.

# Tested only on Linux. Will require changes to run on MacOS.

MODE="burn"
if [ $# -ge 1 ] ; then
    MODE=$1
    if [ "$MODE" == "help" ] ; then
	echo "Usage: $0 [mode]"
	echo "  where the mode is:"
	echo "  burn    – for maximum CPU heating (default, may need adjustments for particular CPU)"
	echo "  scan    – tries few combinations of adjustable parameters to help to determine optimum block size for the burn mode"
	echo "  user    – user-defined: provide n_d and nb as arguments"
	echo "  help    – this help"
	exit 0
    fi
fi

# Count available CPU threads.
NTHR=$( lscpu -p | grep -cv '^#' )
VENDOR=$( grep vendor_id /proc/cpuinfo | uniq | awk '{print $3}' )
NB=32
ND=128
case $MODE in
    "burn")
	# Maximum power draw depends on CPU model and problem configuration:
	# Ryzen 7 2700: ND=144 NB=32
	case $VENDOR in
	    "GenuineIntel")
		# ND=128  # ToDo: run a test on some Intel CPU
		;;
	    "AuthenticAMD")
		ND=144
		;;
	    *)
		echo "Unknown vendor '$VENDOR'. Do a scan first, then make an update."
		exit 1
		;;
	esac
	;;
    "user")
	# Get $nd and $nb from arguments. Bogus values will crash Piernik thus no error checking here.
	if [ $# -ge 2 ] ; then
	    ND=$2
	else
	    echo "You need to provide at least n_d, e.g.:"
	    echo "$0 user 144 32  # equivalent to '$0 burn' on AMD"
	    exit 5
	fi
	[ $# -ge 3 ] && NB=$3
	;;
    "scan")
	# Scan for most energy-consuming setup parameters.
	# It seems that n_d has the biggest impact. Its optimal size may change with size and organization of the CPU cache.

	ND="64 96 128 160 192 224 256"
	NB="16 32 48"
	;;
    *)
	echo "Unknown mode '$MODE'."
	exit 2
	;;
esac

DURATION=604800
[ "$MODE" == "scan" ] && DURATION=30

echo "Mode=$MODE Vendor=$VENDOR Nthr=$NTHR n_d=$ND nb=$NB duration=${DURATION}s"

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

for nb in $NB ; do
    for nd in $ND ; do
	echo "${MODE}: n_d= $nd  nb= $nb" > $SHOW_TEMP
	# Run $NTHR independent single-threaded Piernik copies as this draws most electrical power.
	# Parallel Piernik runs do spend some time on communication, so we choose flooding until we implement shared-memory :)
	for d in $( seq "$NTHR" ) ; do
	    (
		cd "$d" || exit 4
		rm ./*log 2> /dev/null
		./piernik -n '&BASE_DOMAIN nb = '"$nb"' n_d = 2*'"$nd"', 1/' | grep --line-buffered . > out &
	    )
	done

        # ToDo also check if all $NTHR are alive and not stuck.

	# Run each Piernik configuration for limited amount of time.
	sleep $DURATION && killall piernik
	echo "kill $nd $nb"
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

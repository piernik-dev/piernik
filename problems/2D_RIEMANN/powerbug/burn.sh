#!/bin/bash

# Sample script that may help to determine parameters for stress-testing a computer with Piernik.

# Count available CPU threads.
NTHR=$( lscpu -p | grep -v '^#' | wc -l )

# If there is mprime  running for warm-up, stop it.
MPRIME=1
killall -STOP mprime 2> /dev/null || MPRIME=0

SHOW_TEMP=.__st__
# A "daemon" for reporting temperatures. Adjust the grep filter to particular machine.
touch $SHOW_TEMP
sleep 5 && while [ -e $SHOW_TEMP ] ; do
    ( cat $SHOW_TEMP; sensors | grep -E '(Tctl|TSI0_TEMP)' ) | tr '\n' ' ' ; echo
    #( cat $SHOW_TEMP; sensors | grep -E '(Package)' ) | tr '\n' ' ' ; echo
    sleep 5
done &

# Create separate directories for each run.
for i in `seq $NTHR` ; do
    [ -d $i ] && rm -r $i
    mkdir $i
    cp problem.par $i
    ln -s ../piernik ${i}/piernik
done

# Scan for most energy-consuming setup parameters.
# It seems that n_d has the biggest impact. Its optimal size may change with size and organization of the CPU cache.

# Maximum power draw depends on CPU model and problem configuration:
# Ryzen 7 2700: bs=300 solver=RTVD
# i5-11600: bs=116 solver=RTVD

for bs in 116 300 ; do
    k=1
    for s in Riemann RTVD ; do
	nd=$(( $bs * $k ))
	echo "flood: n_d= $nd  solver= $s" > $SHOW_TEMP
	# Run $NTHR independent single-threaded Piernik copies as this draws most electrical power.
	# Parallel Piernik runs do spend some time on communication, so we choose flooding until we implement shared-memory :)
	for d in `seq $NTHR` ; do
	    cd $d
	    rm *log 2> /dev/null
	    ./piernik -n '&BASE_DOMAIN n_d = 2*'$nd', 1/ &AMR bsize = 2*'$bs', 1/ &NUMERICAL_SETUP solver_str="'$s'"' > /dev/null &
	    cd - > /dev/null
	done

	# Check if all instances are reaching the same values.
        OSTEP=0
	while [ -e $SHOW_TEMP ] ; do
	    [ -e 1/*log ] && \
		STEP=$( tail -n 1 [1-9]*/*log | \
			    grep nstep | \
			    sort -n -k 5 | \
			    head -n 1 | \
			    awk '{print $5}' )
            [ ${STEP:-0} -gt 0 ] && \
                if [ $STEP != $OSTEP ] ; then
                    OSTEP=$STEP
                    grep "nstep = *$STEP" [1-9]*/*log | \
			sort -n | \
			column -t | \
			awk '{if (NR==1) {dt=$9; t=$12; print} else if (dt != $9 || t != $12) print "Differs: ", $0}'
                fi
            sleep 3
        done &

        # ToDo also check if all $NTHR are alive and not stuck.

	# Run each Piernik configuration for 60 seconds.
	sleep 60 && killall piernik
	echo "kill $nd $k"
	sleep 1 && killall -9 piernik
	echo
    done
done

rm $SHOW_TEMP

# Continue mprime, if it was there.
[ $MPRIME == 1 ]  && killall -CONT mprime

# This may help to visualize CPU temperature changes.
# gnuplot
# se sty dat lp
# p "< grep Tctl outfile" u 7
# while (1) {rep; pause 3}

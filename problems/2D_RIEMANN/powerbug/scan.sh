#!/bin/bash

# Sample script that may help to determine parameters for stress-testing a computer with Piernik.

# Count available CPU threads.
NTHR=$( lscpu -p | grep -v '^#' | wc -l )

# If there is mprime  running for warm-up, stop it.
killall -STOP mprime

SHOW_TEMP=.__st__
# A "daemon" for reporting temperatures. Adjust the grep filter to particular machine.
touch $SHOW_TEMP
sleep 5 && while [ -e $SHOW_TEMP ] ; do
    ( cat $SHOW_TEMP; sensors | grep -E '(Tctl|TSI0_TEMP)' ) | tr '\n' ' ' ; echo
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

# For a Ryzen 7 2700 maximum power draw occurs at bs=300 and RTVD solver

for bs in 240 300 360 ; do
    k=1
    for s in Riemann RTVD ; do
	nd=$(( $bs * $k ))
	echo "flood: n_d= $nd  solver= $s" > $SHOW_TEMP
	# Run $NTHR independent single-threaded Piernik copies as this draws most electrical power.
	# Parallel Piernik runs do spend some time on communication, so we choose flooding until we implement shared-memory :)
	for d in `seq $NTHR` ; do
	    cd $d
	    ./piernik -n '&BASE_DOMAIN n_d = 2*'$nd', 1/ &AMR bsize = 2*'$bs', 1/ &NUMERICAL_SETUP solver_str="'$s'"' > /dev/null &
	    cd - > /dev/null
	done
	# Run each Piernik configuration for 120 seconds.
	sleep 120 && killall piernik
	echo "kill $nd $k"
	killall -9 piernik
	echo
    done
done

rm $SHOW_TEMP

# Continue mprime, if it was there.
killall -CONT mprime

# This may help to visualize CPU temperature changes.
# gnuplot
# se sty dat lp
# p "< grep Tctl outfile" u 7
# while (1) {rep; pause 3}

#!/bin/bash

if [ $# -lt 1 ] ; then
	echo "Usage: $0 results_directory"
	exit 1
fi

if [ ! -d $1 ] ; then
	echo "trying to do mkdir -p $1"
	mkdir -p $1
fi

if [ ! -d $1 ] ; then
	echo "there is no directory $1"
	exit 2
fi

N=$( ( echo 0 ; \ls -1 $1 ) | grep ^[0-9]*$ | sort -n | tail -n 1 | awk '{print 1*$1}' )

NN=0
while [ $NN == 0 ] ; do
	N=$(( $N + 1 ))
	if [ $N -ge 1000 ] ; then
		echo "Too busy directory"
		exit 3
	fi
	N3=$( printf "%03d" $N )
	[ ! -e $1/$N3 -a ! -e $1/${N3}_big2 ] && NN=$N3
done

echo "starting from $NN"

OPT=""
which mpif90 | grep -q intel && OPT="-c benchmarking_ifx"


for i in $( seq $N $(( $N + 2 )) ) ; do
  ./benchmarking/piernik_bench.sh ${OPT} | tee $1/$( printf "%03d" $i )
done
for b in 1.5 2 ; do
  BIG=$b ./benchmarking/piernik_bench.sh ${OPT} | tee $1/$( printf "%03d" $N )_big$b
done

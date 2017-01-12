#!/bin/bash

if [ $# -lt 1 ] ; then
	echo "Usage: $0 results_directory"
	exit 1
fi

if [ ! -d $1 ] ; then
	echo "trying to do mkdir $1"
	mkdir $1
fi

if [ ! -d $1 ] ; then
	echo "there is no directory $1"
	exit 2
fi

N=$( ( echo 0 ; \ls -1 $1 ) | grep ^[0-9]*$ | sort -n | tail -n 1 )

NN=0
while [ $NN == 0 ] ; do
	N=$(( $N + 1 ))
	if [ $N -ge 1000 ] ; then
		echo "Too busy directory"
		exit 3
	fi
	[ ! -e $1/$N -a ! -e $1/${N}_big2 ] && NN=$N
done

echo "starting from $N"

for i in $( seq $N $(( $N + 2 )) ) ; do
  ./benchmarking/piernik_bench.sh | tee $1/${i}
done
BIG=2 ./benchmarking/piernik_bench.sh  | tee $1/${N}_big2 

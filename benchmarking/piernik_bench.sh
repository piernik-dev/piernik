#!/bin/bash

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


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

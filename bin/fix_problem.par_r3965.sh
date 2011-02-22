#!/bin/bash

# remove this script as soon as everuone adopts r3965

if [ $# -lt 1 ] ; then
	echo "Usage: $0 problem.par [ problem.par [...] ]"
	exit 1
fi

for i in $* ; do
        cp "$i" "$i".backup
	sed -i 's/\$CONSTANTS/$UNITS/;s/constants_set/units_set/' "$i" 
done

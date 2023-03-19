#!/bin/bash

# The .par file extracted from .h5 file matches the problem.par file used for the run
# but it ignores parameters that were overidden by commandline (provided via '-n' option).

# The .par file extracted from .log file should list the parameters that were altered both
# by original problem.par and commandline add-ons.
# It usually differs from the original problem.par due to different sorting of the namelists.

if [ $# -lt 1 ] ; then
	echo "Usage: $0 H5_FILE|LOGFILE [...]"
	exit 1
fi

H5DUMP=h5dump
which $H5DUMP > /dev/null 2>&1 || { echo "Cannot find h5dump executable"; exit 2; }
EPPLOG="${0/.sh/.py}"
echo "$EPPLOG"

for i in $* ; do
	if [ -f "$i" ] ; then
		o=${i}".par"
		if [ -e "$o" ] ; then
			echo "$o already exists. Aborting."
			exit 3
		fi
		if [[ $i == *.log ]]; then
			$EPPLOG "$i" > "$o"
		elif [[ $i == *.h5 ]]; then
			$H5DUMP -d "problem.par" "$i" |  grep ":" | sed 's/ *([0-9]*): "\(.*\)",*/\1/;s/ *$//' > "$o"
		else
			echo "Cannot extract from $i (only .log or .h5 files)"
		fi
		if [ -s "$o" ]; then
			echo "Successfully extracted parameters to '$o' file."
		else
			rm "$o"
		fi
	else
		echo "$i is not a regular file"
	fi
done

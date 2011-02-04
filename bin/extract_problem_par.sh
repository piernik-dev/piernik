#!/bin/bash

if [ $# -lt 1 ] ; then
	echo "Usage: $0 Piernik_HDF5_file.h5 [...]"
	exit 1
fi

H5DUMP=h5dump
which $H5DUMP > /dev/null 2>&1 || { echo "Cannot find h5dump executable"; exit 2; }

for i in $* ; do
	if [ -f "$i" ] ; then
		o=${i}".par"
		if [ -e "$o" ] ; then
			echo "$o already exists. Aborting."
			exit 3
		fi
		$H5DUMP -d "problem.par" "$i" |  grep ":" | sed 's/ *([0-9]*): "\(.*\)",*/\1/;s/ *$//' > "$o" 
		if [ -s "$o" ]; then
			echo "Succesfully extracted parameters to '$o' file."
		else
			rm "$o"
		fi
	else
		echo "$i is not a regular file"
	fi
done

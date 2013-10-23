#!/bin/bash

SEP="===================================================="

for i in `find src -name "*.F90"` ; do
	echo $SEP $i
	grep -rn '\[[a-zA-Z_0-9]*:[a-zA-Z_0-9]*]' $i | grep -v $( basename ${i/.F90/} )":"
done | awk '{ if ($1 !~"'$SEP'") {if (NR == lastsep+1) print "\n"fname; print;} else {lastsep=NR; fname=$2;}}'

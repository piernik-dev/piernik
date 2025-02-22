#!/bin/bash

SEP="===================================================="

NONMATCH=$( for i in `git ls-files | grep -vE "^(compilers/tests|doc/general)" | grep "\.F90$"` ; do
	echo $SEP $i
	grep -rn '\[[a-zA-Z_][a-zA-Z_0-9]*:[a-zA-Z_][a-zA-Z_0-9]*]' $i | grep -v $( basename ${i/.F90/} )":"
done | awk '{ if ($1 !~"'$SEP'") {if (NR == lastsep+1) print "\n"fname; print;} else {lastsep=NR; fname=$2;}}' )

if [ ! -z "$NONMATCH" ]; then
	echo "$NONMATCH"
	exit 1
fi

#!/bin/bash
for i in `find src -name "*.F90"` ; do
	echo "====================================================" $i
	grep -rn '\[[a-zA-Z_0-9]*:[a-zA-Z_0-9]*]' $i | grep -v $( basename ${i/.F90/} )":"
done

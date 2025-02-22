#!/bin/bash

NL=26

md5_head() {
    for i in $( git ls-files | grep -vE "^(compilers/tests|doc/general)" | grep "\.F90$" ) doc/general/new_source_file.F90 ; do
	echo -n $i " "
	head -n $NL $i | md5sum
    done | awk '{print $2,$1}'
}

COUNT=$( md5_head | awk '{print $1}' | sort | uniq -c | wc -l )

if [ $COUNT != 1 ]; then
    DOMINATING=$( md5_head | awk '{print $1}' | sort | uniq -c | sort -n | tail -n 1 | awk '{print $2}' )
    md5_head | grep -v $DOMINATING
    exit 1
fi

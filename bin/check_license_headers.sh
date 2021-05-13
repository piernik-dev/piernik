#!/bin/bash

NL=26

md5_head() {
    for i in $( find src problems  -name *.F90 ) doc/general/new_source_file.F90 ; do
	echo -n $i " "
	head -n $NL $i | md5sum
    done | awk '{print $2,$1}'
}

COUNT=$( md5_head | awk '{print $1}' | sort | uniq -c | wc -l )

if [ $COUNT == 1 ]; then
    echo "All $NL-line license headers are the same"
else
    DOMINATING=$( md5_head | awk '{print $1}' | sort | uniq -c | head -n 1 | awk '{print $2}' )
    echo "Exceptional headers found in:"
    md5_head | grep -v $DOMINATING
    exit 1
fi

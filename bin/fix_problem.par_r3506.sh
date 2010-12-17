#!/bin/bash

# remove this script as soon as everuone adopts r3506

if [ $# -lt 1 ] ; then
	echo "Usage: $0 problem.par [ problem.par [...] ]"
	exit 1
fi

for i in $* ; do
	PP=`mktemp ppXXXXXX`
	awk '{\
  l[NR]=$0;\
  if ($0 ~ "problem_name") pn=NR;\
  if ($0 ~ "run_id") ri=NR;\
  if ($0 ~ "OUTPUT_CONTROL") oc=NR;\
} END {\
  for (i=1;i<=NR;i++) {\
    if (i != pn && i != ri) print l[i];\
    if (i == oc) {\
      print l[pn];\
      print l[ri];\
    }\
  }\
  if (ri==0 || pn ==0 || oc ==0) print "Error: ","'$i'",pn,ri,oc > "/dev/stderr";\
}' "$i" > $PP
	mv "$i" "$i".backup
	mv $PP "$i"
done

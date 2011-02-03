#!/bin/bash

# remove this script as soon as everuone adopts r3871

if [ $# -lt 1 ] ; then
	echo "Usage: $0 problem.par [ problem.par [...] ]"
	exit 1
fi

for i in $* ; do
	PP=`mktemp ppXXXXXX`
	awk '\
BEGIN {\
  px=1; py=1; pz=1; nset=0;\
} {\
  if ($0 ~ "pxsize") { px=$3; nset++; }\
  else if ($0 ~ "pysize") { py=$3; nset++; }\
  else if ($0 ~ "pzsize") { pz=$3; nset++; }\
  else {\
    if (nset == 0 || $0 ~ "reorder") print;\
    else if ($0 ~ "/") {\
      printf "    psize = %d, %d, %d\n", px, py, pz;\
      print;\
      nset = 0;\
    }\
  }\
}' "$i" > $PP
	mv "$i" "$i".backup
	mv $PP "$i"
done

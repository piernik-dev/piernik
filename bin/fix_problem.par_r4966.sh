#!/bin/bash

# remove this script as soon as everuone adopts r4966

if [ $# -lt 1 ] ; then
	echo "Usage: $0 problem.par [ problem.par [...] ]"
	exit 1
fi

for i in $* ; do
	PP=`mktemp ppXXXXXX`
	awk '\
BEGIN {\
  nx=1; ny=1; nz=1; nset=0;\
} {\
  if      ($0 ~ "nxd" && NF==3) { nx=$3; nset++; }\
  else if ($0 ~ "nyd" && NF==3) { ny=$3; nset++; }\
  else if ($0 ~ "nzd" && NF==3) { nz=$3; nset++; }\
  else {\
    if (nset == 3) {\
      printf "    n_d = %d, %d, %d\n", nx, ny, nz;\
      nset = 0;\
    }\
    print;\
  }\
}' "$i" > $PP
	mv "$i" "$i".backup
	mv $PP "$i"
done

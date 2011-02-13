#!/bin/bash

# remove this script as soon as everuone adopts r3912

if [ $# -lt 1 ] ; then
	echo "Usage: $0 problem.par [ problem.par [...] ]"
	exit 1
fi

for i in $* ; do
	PP=`mktemp ppXXXXXX`
	awk '\
BEGIN {\
  dom=0;\
  ld=0;\
  ln=0;\
} {\
  if (dom >0) {\
    if ($1 == "/") dom=0;\
    else domnml[++ld] = $0;\
  }\
  else {\
    if ($0 ~ "DOMAIN_SIZES" || $0 ~ "BOUNDARIES" || $0 ~ "DOMAIN_LIMITS") dom=1;\
    else nml[++ln] = $0;\
  }\
} END {\
  print " $DOMAIN";\
  for (i=1;i<=ld;i++) print domnml[i];\
  print " /";\
  for (i=1;i<=ln;i++)
    if (length(nml[i])>0) {\
      if (nml[i] ~ "\\$") print "";\
      print nml[i];\
    }\
}' "$i" > $PP
	mv "$i" "$i".backup
	mv $PP "$i"
done

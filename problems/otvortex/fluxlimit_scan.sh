#!/bin/bash

for f in vanleer minmod moncen superbee ; do
	for b in vanleer minmod moncen superbee ; do
		d=${f}_${b}
		mkdir $d
		cd $d
		cp ../problem.par .
		ln -s ../piernik
		./piernik -n '&NUMERICAL_SETUP limiter="'$f'" limiter_b="'$b'"/' &
		cd -
	done
done
wait
for i in *_*/*log ; do
	echo -n $i
	grep nstep $i | tail -n 1
done | sed 's-/- -'  | sort | column -t

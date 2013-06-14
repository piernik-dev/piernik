#!/bin/bash

LOG=_log_
JEANS=jeans_C_jmg

rm $JEANS/{verify.gpl,jeans-*.png} 2> /dev/null
for i in *_C_* ; do
	cd $i
	rm *.h5 *.res *.tsl *.log *.plt 2> /dev/null
	mpirun -np 13 ./piernik 2>&1 | grep -v "Fortran runtime warning: IEEE 'denormal number' exception not supported."
	cd - > /dev/null
done | tee $LOG

cd $JEANS
gnuplot verify.gpl
display jeans-mg.png &
cd - > /dev/null

grep -E "(Starting problem|L2 error norm)" $LOG
rm $LOG

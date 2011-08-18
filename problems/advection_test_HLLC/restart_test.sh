#!/bin/bash

PP="problem.par"
PPRT=${PP}".restart_test"

OUTB="moving_pulse_ts1_"
OUT2=${OUTB}"0002.h5"
OUT1=${OUTB}"0001.h5"
OUT2_1=${OUTB}"1002.h5"
OUT1_1=${OUTB}"1001.h5"

if [ ! -e $PPRT ] ; then
    echo "Cannot find $PPRT. Trying to find one up to ../../."
    find ../../ -name $PPRT
    echo "You must manually copy an appropriate $PPRT to ."
    exit 1
fi

PRG="./piernik"

if [ ! -x $PRG ] ; then
    echo "Cannot find Piernik executable"
    exit 2
fi

rm *.res 2> /dev/null

cp $PPRT $PP
$PRG
cp $OUT1 $OUT1_1
sed '/tend/s/.*/ tend = 2.0/' $PPRT > $PP
$PRG
cp $OUT2 $OUT2_1

rm *.res
cp $PPRT $PP
mpirun -np 5 $PRG
echo "Comparing $OUT1 and $OUT1_1"
h5diff $OUT1 $OUT1_1

sed '/tend/s/.*/ tend = 2.0/' $PPRT > $PP
mpirun -np 9 $PRG

echo "Comparing $OUT2 and $OUT2_1"
h5diff $OUT2 $OUT2_1

#rm *.res moving_pulse_ts1_????.h5

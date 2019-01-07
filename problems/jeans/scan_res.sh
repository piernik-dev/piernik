#!/bin/bash

# This script performs a number of 1D Jeans tests with varying resolution and CFL.
# For each run we determine the amplitude, frequency and damping factor for the oscillations.
# Then we can determine how these parameters depend on resolution and timestep by fitting simple polynomial function across the runs

SCAN=$( mktemp scan_XXXXXX )

echo "#res cfl aa bb a b c maxval maxres maxdiff" >> $SCAN

for cfl in 0.05 0.07 0.1 0.15 0.2 0.3 0.4 0.6 ; do
    for i in `seq 14 23` ; do
        r=$(( 2 ** ($i / 2) * ( 2 + ($i % 2)) / 2 ))
        mpirun -np $(cat /proc/cpuinfo | grep "cpu cores" | head -n 1 | sed 's/.*: //') ./piernik \
            -n '&BASE_DOMAIN n_d = '$r', 1, 1 / &OUTPUT_CONTROL run_id = "'$( printf "t%02d" $i )'"/ &NUMERICAL_SETUP cfl = '$cfl'/'
        gnuplot jeans.gnuplot 2>&1 | awk 'BEGIN {printf("%5s %-5s ", '$r','$cfl')} /#/ {for (i=3;i<=NF;i++) printf("%s ",$i)} END {print ""}' | tee -a $SCAN
        mv jeans.png jeans_${r}_${cfl}.png
    done
    echo "" >> $SCAN
done
echo "Logged to $SCAN"

fit_poly() {
rm fit.log 2> /dev/null
sleep .1
gnuplot << EOF
f0(x,y) = a
f1(x,y) = f0(x,y) + ax*x       + ay*y
f2(x,y) = f1(x,y) + axx*x**2   + axy*x*y      + ayy*y**2
f3(x,y) = f2(x,y) + axxx*x**3  + axxy*x**2*y  + axyy*x*y**2     + ayyy*y**3
f4(x,y) = f3(x,y) + axxxx*x**4 + axxxy*x**3*y + axxyy*x**2*y**2 + axyyy*x*y**3 + ayyyy*y**4
fit f4(x,y) "$SCAN" u (1/\$1):2:(column($c)) via a, ax,ay, axx,axy,ayy, axxx,axxy,axyy,ayyy, axxxx,axxxy,axxyy,axyyy,ayyyy
se te png size 1024, 768
se sty dat lp
se vi 60, 330
se ou "${SCAN}.${c}.png"
sp "$SCAN" u (1/\$1):2:(column($c)), f4(x,y)
se ou
EOF

mv fit.log ${SCAN}"."${c}".fit"
}

FIT=$( mktemp fit_XXXXXX )
(
    echo "# f(dx,cfl) = a           +"
    echo "#             ax*dx       + ay*cfl"
    echo "#             axx*dx**2   + axy*dx*cfl      + ayy*cfl**2         +"
    echo "#             axxx*dx**3  + axxy*dx**2*cfl  + axyy*dx*cfl**2     + ayyy*cfl**3     +"
    echo "#             axxxx*dx**4 + axxxy*dx**3*cfl + axxyy*dx**2*cfl**2 + axyyy*dx*cfl**3 + ayyyy*cfl**4"
    echo "   "
    echo "#column a ax ay axx axy ayy axxx axxy axyy axyy axxxx axxxy axxyy axyyy ayyyy"
) >> $FIT
for c in 5 6 7 8 9 10 ; do
    fit_poly $SCAN $c
    awk 'BEGIN {printf("%s ",'$c')} /a[xy]* *=/ {printf("%s ",$3)} END {print ""}' ${SCAN}"."${c}".fit" >> $FIT
done

column -t $FIT
echo "Fit logged to $FIT"

# Try to find which are the leading terms
# awk 'BEGIN {x=1/64.; y=.8} {if (NF==16) if ($0 ~ "column") print; else print $1, $2, $3*x, $4*y, $5*x**2, $6*x*y, $7*y**2, $8*x**3, $9*x**2*y, $10*x*y**2, $11*y**3, $12*x**4, $13*x**3*y, $14*x**2*y**2, $15*x*y**3, $16*y**4}' $FIT | column -t

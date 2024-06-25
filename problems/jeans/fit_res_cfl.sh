#!/bin/bash

# This script performs a number of 1D Jeans tests with varying resolution and CFL.
# For each run we determine the amplitude, frequency and damping factor for the oscillations.
# Then we can determine how these parameters depend on resolution and timestep by fitting simple polynomial function across the runs

SCAN=$( mktemp scan_XXXXXX )

echo "#res cfl aa bb a b c maxval maxres maxdiff" >> $SCAN

for cfl in 0.1 0.2 0.3 0.5 0.7; do
    for i in `seq 12 18` ; do
        r=$(( 2 ** ($i / 2) * ( 2 + ($i % 2)) / 2 ))  # resolution: 64 96 128 192 256 384 512
        ryz=8
        xmax=$( grep xmax problem.par | awk '{print $3}'  )
        yzmax=$( awk 'BEGIN {print '$xmax' * '$ryz'/'$r'}' )
        mpirun -np $( cat /proc/cpuinfo | grep "cpu cores" | head -n 1 | sed 's/.*: //') ./piernik \
            -n '&BASE_DOMAIN n_d = '$r', '$ryz', '$ryz' ymax = '$yzmax' zmax = '$yzmax' / &OUTPUT_CONTROL run_id = "'$( printf "t%02d" $i )'"/ &NUMERICAL_SETUP cfl = '$cfl'/ &AMR bsize = 32, '$ryz', '$ryz' /'
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
se xla "dx"
se yla "CFL"
se vi 60, 330
se ou "${SCAN}.${c}.png"
sp "$SCAN" u (1/\$1):2:(column($c)) t "column($c): $vname", f4(x,y) t "4th order fit"
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
for c in 5 6 7 ; do
    case $c in
	(5) vname="amplitude";;
	(6) vname="period";;
	(7) vname="damping coefficient";;
	(*) vname="???";;
    esac
    fit_poly $SCAN $c $vname
    awk 'BEGIN {printf("%s ",'$c')} /a[xy]* *=/ {printf("%s ",$3)} END {print ""}' ${SCAN}"."${c}".fit" >> $FIT
done

column -t $FIT
echo "Fit logged to $FIT"

# Try to find which are the leading terms
awk 'BEGIN {\
    x = 1/64.;\
    y = .8;\
    white = "\033[97m";\
    gray = "\033[90m";\
    normal = "\033[0m";\
    high = 0.9;\
    low = 0.01;\
    printf("# Relative importance of coefficients (y = CFL = %s, x = Δx = %s = 1/%d)\n", y, x, int(1/x+.1));\
    val[5] = "5:amplitude";\
    val[6] = "6:period";\
    val[7] = "7:damping";\
} function abs(x) {\
  return (x > 0) ? x : -x;\
} {\
  if (NF == 16)\
     if ($0 ~ "column") {\
        printf("%-11s ", $1);
	for (i=2; i<=NF; i++)\
	    printf("%13s ", $i);\
        printf("\n");\
     }\
     else {\
          a[0] = $1;\
	  a[1] = $2;\
	  a[2] = $3*x;\
	  a[3] = $4*y;\
	  a[4] = $5*x**2;\
	  a[5] = $6*x*y;\
	  a[6] = $7*y**2;\
	  a[7] = $8*x**3;\
	  a[8] = $9*x**2*y;\
	  a[9] = $10*x*y**2;\
	  a[10] = $11*y**3;\
	  a[11] = $12*x**4;\
	  a[12] = $13*x**3*y;\
	  a[13] = $14*x**2*y**2;\
	  a[14] = $15*x*y**3;\
	  a[15] = $16*y**4;\
	  max = 0.;
	  second = 0.;\
	  for (i=1 ;i<=15; i++)\
	      if (abs(a[i]) > max)\
	          max = a[i];\
	      else if (abs(a[i]) > second)\
	          second = a[i];\
          printf("%-11s ", val[a[0]]);\
	  for (i=1 ;i<=15; i++) {\
	      color = normal;\
	      if (abs(a[i]) >= high * max)\
	          color = white;\
	      if (abs(a[i]) < low * max && abs(a[i]) < high * second)\
	          color = gray;\
	      printf("%s%13s%s ", color, a[i], normal);\
          }\
	  printf("\n");\
      }\
 }' $FIT

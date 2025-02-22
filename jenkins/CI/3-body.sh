echo "" > 3body.csv

[ ! -e problems/2body/3body ] && exit 0
CONF=lothlorien
[ -e compilers/jenkins.in ] && CONF=jenkins

python setup 2body/3body -c $CONF
cd runs/3body

#restart test
RES=restart
touch $RES
rm -r $RES
mkdir -p $RES
(
	cd $RES
    cp ../problem.par .
    sed -i '/bgfile/s/..\/../..\/..\/../' problem.par
    ln -s ../piernik
	touch ___dummy___.res
	rm *.res
	./piernik -n '$END_CONTROL nend = 10/ $OUTPUT_CONTROL run_id = "rs1"/'
	./piernik -n '$END_CONTROL nend = 20/ $OUTPUT_CONTROL run_id = "rs1"/'
	./piernik -n '$END_CONTROL nend = 20/ $OUTPUT_CONTROL run_id = "rs2"/'

	../../../bin/gdf_distance leapfrog_rs{1,2}_0001.h5 | tee compare.log
) &

touch ___dummy___.res
rm *.res
./piernik
../../problems/2body/3body/particle_error.py leapfrog_tst_0001.h5 > 3body.csv

wait # restart
[ $( grep "^Total difference between" ${RES}/compare.log | awk '{print $NF}' ) == 0 ] || exit 1

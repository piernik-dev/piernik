if [ -e runs/jeans ] ; then find runs/jeans/ -name "*.png" -exec \rm {} \; ; fi
CONF=lothlorien
[ -e compilers/jenkins.in ] && CONF=jenkins

# This should go away once we upgrade Jenkins
[ -e compilers/tests/mpi_f08.F90 ] && sed -i 's/use mpi_f08/use mpi_f08, only: MPI_DOUBLE_INT/' compilers/tests/mpi_f08.F90

python setup jeans -c $CONF
python make_test.py --test=jeans

#restart test
cd runs/jeans
touch ___dummy___.res
rm *.res
./piernik -n '$END_CONTROL nend = 10/ $OUTPUT_CONTROL dt_hdf = 1.0 dt_res = 1.0 run_id = "rs1"/'
./piernik -n '$END_CONTROL nend = 20/ $OUTPUT_CONTROL dt_hdf = 1.0 dt_res = 1.0 run_id = "rs1"/'
./piernik -n '$END_CONTROL nend = 20/ $OUTPUT_CONTROL dt_hdf = 1.0 dt_res = 1.0 run_id = "rs2"/'

../../bin/gdf_distance jeans_rs{1,2}_0001.h5 | tee compare.log
[ $( grep "^Total difference between" ./compare.log | awk '{print $NF}' ) == 0 ] || exit 1

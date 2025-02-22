git log -3
. /etc/profile.d/modules.sh
module load mpi
CONF=lothlorien
[ -e compilers/jenkins.in ] && CONF=jenkins

# This should go away once we upgrade Jenkins
[ -e compilers/tests/mpi_f08.F90 ] && sed -i 's/use mpi_f08/use mpi_f08, only: MPI_DOUBLE_INT/' compilers/tests/mpi_f08.F90

python setup maclaurin -c $CONF
python make_test.py --test=maclaurin
echo "L2 error norm,min. error,max. error" > norm.csv


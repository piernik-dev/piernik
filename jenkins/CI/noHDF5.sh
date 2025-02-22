git log -3
. /etc/profile.d/modules.sh
module load mpi
python setup crtest -c jenkins --param problem.par.build -d I_KNOW_WHAT_I_AM_DOING


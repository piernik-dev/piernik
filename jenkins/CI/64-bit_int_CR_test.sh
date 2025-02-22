git log -3
. /etc/profile.d/modules.sh
module load mpi
CONF=lothlorien
[ -e compilers/jenkins.in ] && CONF=jenkins

if [ -e problems/testing_and_debuging/chimaera ] ; then
    I64=1 python setup testing_and_debuging/chimaera -c $CONF
else
    I64=1 python setup crtest -c $CONF --param problem.par.build -d VERBOSE
fi


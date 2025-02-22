git log -3
. /etc/profile.d/modules.sh
module load mpi
CONF=lothlorien
[ -e compilers/jenkins.in ] && CONF=jenkins

echo "-c $CONF" > .setuprc

bash problems/advection_test/restart_test_v2_jenkins.sh
if [ -e compare.log ] ; then
	[ $( grep "^Total difference between" compare.log | awk '{print $NF}' ) == 0 ] || exit 1
fi


git log -3
. /etc/profile.d/modules.sh
module load mpi

./jenkins/gold_test.sh ./jenkins/gold_configs/mcrtest_new.config


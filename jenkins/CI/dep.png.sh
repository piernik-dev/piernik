git log -3
ls -l

PROB=crtest
[ -e problems/testing_and_debuging/chimaera ] && PROB=testing_and_debuging/chimaera

./setup $PROB -n
make -C obj dep.png


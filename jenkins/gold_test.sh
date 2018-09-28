#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 config_file"
    echo "or somethink like: RUN_COMMAND=\"mpirun -np 2\" $0 config_file"
    exit 1
fi

for i in GOLD_COMMIT PROBLEM_NAME SETUP_PARAMS GOLD_PARAMS OUTPUT; do
    eval $( sed -n '/^'"$i"'=/p' $1 )
done

if [ -z ${GOLD_COMMIT+x} ] ; then
    echo "GOLD_COMMIT not set"
    exit 2
fi

if [ -z ${PROBLEM_NAME+x} ] ; then
    echo "PROBLEM_NAME not set"
    exit 3
fi

if [ -z ${SETUP_PARAMS+x} ] ; then
    echo "SETUP_PARAMS not set"
    exit 4
fi
SETUP_PARAMS=$SETUP_PARAMS" -n --copy"

if [ -z ${GOLD_PARAMS+x} ] ; then
    echo "GOLD_PARAMS not set"
    exit 5
fi

if [ -z ${OUTPUT+x} ] ; then
    echo "OUTPUT not set"
    exit 5
fi

PIERNIK_REPO="http://github.com/piernik-dev/piernik"
PIERNIK=piernik
GOLD_DIR=gold_dir
OBJ_PREFIX=obj_
GOLD_OBJ=${PROBLEM_NAME}_gold
TEST_OBJ=${PROBLEM_NAME}_test
RIEM_OBJ=${PROBLEM_NAME}_Rtest
GOLD_LOG=gold_log
RIEM_LOG=riem_log
RIEM_CSV=riem.csv
TMP_DIR=/tmp/jenkins_gold/
RUNS_DIR=$TMP_DIR
GOLD_SHA_FILE=__sha__

[ ! -d $TMP_DIR ] && mkdir -p $TMP_DIR
cp Makefile $TMP_DIR

rm -rf $GOLD_DIR ${OBJ_PREFIX}$GOLD_OBJ ${OBJ_PREFIX}$TEST_OBJ ${RUNS_DIR}/${PROBLEM_NAME}_$TEST_OBJ $GOLD_LOG

for i in $GOLD_OBJ $TEST_OBJ $RIEM_OBJ ; do
    mkdir ${RUNS_DIR}/${PROBLEM_NAME}_$i
done

python setup $PROBLEM_NAME $SETUP_PARAMS -o $TEST_OBJ
python setup $PROBLEM_NAME $SETUP_PARAMS -d RIEMANN -o $RIEM_OBJ
for i in $TEST_OBJ $RIEM_OBJ ; do
    rsync -Icvxa --no-t ${OBJ_PREFIX}$i $TMP_DIR
done

git clone $PIERNIK_REPO $GOLD_DIR
[ -e .setuprc ] && cp .setuprc $GOLD_DIR
cp python/piernik_setup.py ${GOLD_DIR}/python
(
    cd $GOLD_DIR
    git fetch $PIERNIK_REPO +refs/pull/*:refs/remotes/origin/pr/*
    git checkout $GOLD_COMMIT
    rsync -avx --delete ../compilers/ ./compilers
    python setup $PROBLEM_NAME $SETUP_PARAMS -o $GOLD_OBJ
)
rsync -Icvxa --no-t ${GOLD_DIR}/${OBJ_PREFIX}$GOLD_OBJ $TMP_DIR

cp ${GOLD_DIR}/runs/${PROBLEM_NAME}_${GOLD_OBJ}/problem.par ${RUNS_DIR}/${PROBLEM_NAME}_${GOLD_OBJ}/
cp runs/${PROBLEM_NAME}_${TEST_OBJ}/problem.par ${RUNS_DIR}/${PROBLEM_NAME}_${TEST_OBJ}/
cp runs/${PROBLEM_NAME}_${RIEM_OBJ}/problem.par ${RUNS_DIR}/${PROBLEM_NAME}_${RIEM_OBJ}/

(
    cd $TMP_DIR
    rm ${OBJ_PREFIX}{$GOLD_OBJ,$TEST_OBJ,$RIEM_OBJ}/piernik
    make -j ${OBJ_PREFIX}{$GOLD_OBJ,$TEST_OBJ,$RIEM_OBJ}
)

for i in $GOLD_OBJ $TEST_OBJ $RIEM_OBJ ; do
    rm ${RUNS_DIR}/${PROBLEM_NAME}_${i}/${PIERNIK} 2> /dev/null || echo "rundir clean"
    cp ${TMP_DIR}/${OBJ_PREFIX}${i}/${PIERNIK} ${RUNS_DIR}/${PROBLEM_NAME}_$i
done

for i in $TEST_OBJ $RIEM_OBJ ; do
    (
	cd ${RUNS_DIR}/${PROBLEM_NAME}_$i
	rm *.h5 2> /dev/null
	eval $RUN_COMMAND ./${PIERNIK}
    )
done

(
    cd ${RUNS_DIR}/${PROBLEM_NAME}_$GOLD_OBJ
    # skip recalculating gold outpuf if and only if we know for sure that $GOLD_COMMIT did not change
    if [ ! -e $GOLD_SHA_FILE ] ; then
	[ -e ${OUTPUT} ] && rm ${OUTPUT}
    else
	if [ $( cat $GOLD_SHA_FILE  ) != $GOLD_COMMIT ] ; then
	   [ -e ${OUTPUT} ] && rm ${OUTPUT}
	fi
    fi
    if [ ! -e ${OUTPUT} ] ; then
       eval $RUN_COMMAND ./${PIERNIK} $GOLD_PARAMS
       echo $GOLD_COMMIT > $GOLD_SHA_FILE
    fi
)

[ ! -z $YT ] && source $YT
./bin/gdf_distance ${RUNS_DIR}/${PROBLEM_NAME}_{${TEST_OBJ},${GOLD_OBJ}}/${OUTPUT} 2>&1 | tee $GOLD_LOG

if [ -e ${RUNS_DIR}/${PROBLEM_NAME}_${RIEM_OBJ}/${OUTPUT} ] ; then
   ./bin/gdf_distance ${RUNS_DIR}/${PROBLEM_NAME}_{${TEST_OBJ},${RIEM_OBJ}}/${OUTPUT} 2>&1 | tee $RIEM_LOG

   grep 'Difference of datafield `' $RIEM_LOG |\
       sed 's/.*`\([^ ]*\)[^ ] *: \(.*\)/\1 \2/' |\
       awk '{a[$1]=$2} END {for (i in a) printf("%s,", i); print ""; for (i in a) printf("%s,", a[i]) ; print ""}' |\
       sed 's/,$//' |\
       tee $RIEM_CSV
else
    (
	echo "Riemann failed"
	echo 0
    ) | tee $RIEM_CSV
fi

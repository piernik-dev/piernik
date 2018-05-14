#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 config_file"
    echo "or somethink like: RUN_COMMAND=\"mpirun -np 2\" $0 config_file"
    exit 1
fi

for i in GOLD_COMMIT PROBLEM_NAME SETUP_PARAMS GOLD_PARAMS OUTPUT; do
    eval $( sed -n '/^'"$i"'=/p' $1 )
    if [ -z ${!i+x} ] ; then
	echo "$i not set"
	exit 2
    else
	echo $i"  =  "${!i}
    fi
done

PIERNIK_REPO="http://github.com/piernik-dev/piernik"
PIERNIK=piernik
GOLD_DIR=${PROBLEM_NAME}_gold_dir
OBJ_PREFIX=obj_
GOLD_OBJ=${PROBLEM_NAME}_gold
TEST_OBJ=${PROBLEM_NAME}_test
GOLD_LOG=${PROBLEM_NAME}.gold_log
TMP_DIR=/tmp/jenkins_gold/
RUNS_DIR=$TMP_DIR
GOLD_SHA_FILE=__sha__

[ ! -d $TMP_DIR ] && mkdir -p $TMP_DIR
cp Makefile $TMP_DIR

rm -rf $GOLD_DIR ${OBJ_PREFIX}$GOLD_OBJ ${OBJ_PREFIX}$TEST_OBJ ${RUNS_DIR}/${PROBLEM_NAME}_$TEST_OBJ $GOLD_LOG

for i in ${RUNS_DIR}/${PROBLEM_NAME}_${TEST_OBJ} ${RUNS_DIR}/${PROBLEM_NAME}_${GOLD_OBJ} ; do
    [ ! -d "$i" ] && mkdir "$i"
done

python setup $PROBLEM_NAME $SETUP_PARAMS -n --copy -o $TEST_OBJ
rsync -Icvxaq --no-t ${OBJ_PREFIX}$TEST_OBJ $TMP_DIR

git clone -q $PIERNIK_REPO $GOLD_DIR
[ -e .setuprc ] && cp .setuprc $GOLD_DIR
cp python/piernik_setup.py ${GOLD_DIR}/python
(
    cd $GOLD_DIR
    git fetch -q $PIERNIK_REPO +refs/pull/*:refs/remotes/origin/pr/*
    git checkout -q $GOLD_COMMIT
    rsync -avxq --delete ../compilers/ ./compilers
    python setup $PROBLEM_NAME $SETUP_PARAMS -n --copy -o $GOLD_OBJ
)
rsync -Icvxaq --no-t ${GOLD_DIR}/${OBJ_PREFIX}$GOLD_OBJ $TMP_DIR

cp ${GOLD_DIR}/runs/${PROBLEM_NAME}_${GOLD_OBJ}/problem.par ${RUNS_DIR}/${PROBLEM_NAME}_${GOLD_OBJ}/
cp runs/${PROBLEM_NAME}_${TEST_OBJ}/problem.par ${RUNS_DIR}/${PROBLEM_NAME}_${TEST_OBJ}/

(
    cd $TMP_DIR
    sed -i 's/-fcheck=all/& -fcheck=no-array-temps/' ${OBJ_PREFIX}$GOLD_OBJ"/Makefile" ${OBJ_PREFIX}$TEST_OBJ"/Makefile"
    make -j ${OBJ_PREFIX}$GOLD_OBJ ${OBJ_PREFIX}$TEST_OBJ > ${PROBLEM_NAME}.make_stdout 2>&1
)

for i in $GOLD_OBJ $TEST_OBJ ; do
    rm ${RUNS_DIR}/${PROBLEM_NAME}_${i}/${PIERNIK} 2> /dev/null || echo "${PROBLEM_NAME}: rundir clean"
    cp ${TMP_DIR}/${OBJ_PREFIX}${i}/${PIERNIK} ${RUNS_DIR}/${PROBLEM_NAME}_$i
done

(
    cd ${RUNS_DIR}/${PROBLEM_NAME}_$TEST_OBJ
    eval $RUN_COMMAND ./${PIERNIK} > ${PROBLEM_NAME}.test_stdout
) &

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
       eval $RUN_COMMAND ./${PIERNIK} $GOLD_PARAMS > ${PROBLEM_NAME}.gold_stdout
       echo $GOLD_COMMIT > $GOLD_SHA_FILE
    fi
) &

wait

[ ! -z $YT ] && source $YT
./bin/gdf_distance ${RUNS_DIR}/${PROBLEM_NAME}_{${TEST_OBJ},${GOLD_OBJ}}/${OUTPUT} 2> /dev/null | tee $GOLD_LOG

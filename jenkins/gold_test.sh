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

if [ -z ${GOLD_PARAMS+x} ] ; then
    echo "GOLD_PARAMS not set"
    exit 5
fi

if [ -z ${OUTPUT+x} ] ; then
    echo "OUTPUT not set"
    exit 5
fi

PIERNIK=piernik
GOLD_DIR=gold_dir
RUNS_DIR=runs
OBJ_PREFIX=obj_
GOLD_OBJ=gold
TEST_OBJ=test
GOLD_LOG=gold_log

rm -rf $GOLD_DIR ${OBJ_PREFIX}* ${RUNS_DIR}/* $GOLD_LOG

python setup $PROBLEM_NAME $SETUP_PARAMS -n --copy -o $TEST_OBJ

git clone http://github.com/piernik-dev/piernik $GOLD_DIR
[ -e .setuprc ] && cp .setuprc $GOLD_DIR
(
    cd $GOLD_DIR
    git checkout $GOLD_COMMIT
    rsync -avx --delete ../compilers/ ./compilers
    python setup $PROBLEM_NAME $SETUP_PARAMS -n --copy -o $GOLD_OBJ
)
mv ${GOLD_DIR}/${OBJ_PREFIX}$GOLD_OBJ .
mv ${GOLD_DIR}/${RUNS_DIR}/* ${RUNS_DIR}

make -j

for i in $GOLD_OBJ $TEST_OBJ ; do
    cp ${OBJ_PREFIX}${i}/${PIERNIK} ${RUNS_DIR}/${PROBLEM_NAME}_$i
done

(
    cd ${RUNS_DIR}/${PROBLEM_NAME}_$TEST_OBJ
    eval $RUN_COMMAND ./${PIERNIK}
)

(
    cd ${RUNS_DIR}/${PROBLEM_NAME}_$GOLD_OBJ
    eval $RUN_COMMAND ./${PIERNIK} $GOLD_PARAMS
)

[ ! -z $YT ] && source $YT
./bin/gdf_distance ${RUNS_DIR}/${PROBLEM_NAME}_{${TEST_OBJ},${GOLD_OBJ}}/${OUTPUT} 2>&1 | tee $GOLD_LOG

#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage: $0 config_file"
    exit 1
fi

NTHR=1  # serial run by default

# Read parameters from given config file
for i in GOLD_COMMIT PROBLEM_NAME SETUP_PARAMS GOLD_PARAMS OUTPUT NTHR; do
    eval $( sed -n '/^'"$i"'=/p' $1 )
    if [ -z ${!i+x} ] ; then
	echo "$i not set"
	exit 2
    else
	echo $i"  =  "${!i}
    fi
done
FLAT_PROBLEM_NAME=${PROBLEM_NAME//\//___}  # handle subproblems safely

OVERSUB=$( mpirun --oversubscribe test true 2> /dev/null && echo "--oversubscribe" || echo "" )
RUN_COMMAND="mpirun -np $NTHR $OVERSUB"
SETUP_PARAMS=$SETUP_PARAMS" -n --copy --linkexe"
PIERNIK_REPO="http://github.com/piernik-dev/piernik"
PIERNIK=piernik
BASE_DIR=$( pwd )
OUT_DIR=jenkins/workspace/${FLAT_PROBLEM_NAME}/
GOLD_DIR=${OUT_DIR}gold/
GOLD_CLONE=${GOLD_DIR}piernik_gold
TEST_DIR=${OUT_DIR}test/
OBJ=obj_${FLAT_PROBLEM_NAME}
RUN_TEST_DIR=${TEST_DIR}runs/$(basename ${PROBLEM_NAME})/
RUN_GOLD_DIR=${GOLD_DIR}runs/$(basename ${PROBLEM_NAME})/

GOLD_LOG=${OUT_DIR}gold_log
GOLD_CSV=${OUT_DIR}gold.csv
RIEM_LOG=${OUT_DIR}riem_log
RIEM_CSV=${OUT_DIR}riem.csv
GOLD_SHA_FILE=${OUT_DIR}__sha__

# Determine whether gold version has to be recreated
# Check the content of the $GOLD_SHA_FILE and existence of the output file to be compared to
if [ ! -e $GOLD_SHA_FILE ] ; then
    [ -e ${RUN_GOLD_DIR}${OUTPUT} ] && rm ${RUN_GOLD_DIR}${OUTPUT}
else
    if [ $( cat $GOLD_SHA_FILE  ) != $GOLD_COMMIT ] ; then
	[ -e ${RUN_GOLD_DIR}${OUTPUT} ] && rm ${RUN_GOLD_DIR}${OUTPUT}
    fi
fi

#Prevent old problem.par from messing up the tests
HAS_KEEPPAR=0
if [ -e .setuprc ] ; then
    grep -q keeppar .setuprc && HAS_KEEPPAR=1
    sed -i 's/--keeppar//g' .setuprc
fi

# Run the gold test if needed (in the background)
if [ ! -e ${RUN_GOLD_DIR}${OUTPUT} ] ; then
    echo "(Re)creating the gold"

    # Clean up after previous runs and create directories
    # Do not try to reuse previous version of object or run directory to avoid weird problems
    # This is supposed to be rare operation unless you call "make gold-clean" often
    rm -rf $GOLD_DIR $GOLD_SHA_FILE
    mkdir -p $GOLD_DIR

    # Prepare object from gold SHA1 ($GOLD_COMMIT)
    # Use current setup script and compilers configuration
    git clone -q $PIERNIK_REPO $GOLD_CLONE
    [ -e .setuprc ] && cp .setuprc $GOLD_CLONE
    cp python/piernik_setup.py ${GOLD_CLONE}/python/piernik_setup_today.py
    (
	cd $GOLD_CLONE
	git fetch -q $PIERNIK_REPO +refs/pull/*:refs/remotes/origin/pr/*
	git checkout -q $GOLD_COMMIT
	rsync -avxq --delete "${BASE_DIR}"/compilers/ ./compilers
	python3 python/piernik_setup_today.py $PROBLEM_NAME $SETUP_PARAMS -o $FLAT_PROBLEM_NAME  # prepares $OBJ in $GOLD_DIR
    )
    cp -a ${GOLD_CLONE}/$OBJ $GOLD_DIR

    # Copy whole directories from appropriate runs directories in case there is something more than piernik and problem.par
    mkdir -p $RUN_GOLD_DIR
    cp -a ${GOLD_CLONE}/runs/$(basename ${PROBLEM_NAME})_${FLAT_PROBLEM_NAME}/* $RUN_GOLD_DIR

    # Compile the gold version of Piernik
    sed -i 's/-fcheck=all/& -fcheck=no-array-temps/' ${GOLD_DIR}$OBJ"/Makefile"
    make -j -C ${GOLD_DIR}$OBJ > ${GOLD_DIR}make.stdout 2>&1

    # Run the gold version of Piernik problem
    (
	cd $RUN_GOLD_DIR
	eval $RUN_COMMAND ./${PIERNIK} $GOLD_PARAMS > gold_stdout
    )

    # Leave a mark to avoid unnecessary recalculation of the gold problem
    echo $GOLD_COMMIT > $GOLD_SHA_FILE
fi &

# Clean up after previous runs and create directories
rm -rf ${TEST_DIR}$OBJ ${TEST_DIR}runs
mkdir -p $TEST_DIR

# Prepare object from current files (including any uncommited or untracked changes)
python3 setup $PROBLEM_NAME $SETUP_PARAMS -o $FLAT_PROBLEM_NAME  # prepares $OBJ in $BASE_DIR
mv ${OBJ} ${TEST_DIR}
# OPT: By careful updating ${TEST_DIR}${OBJ} one can achieve some speedups on compilation

# Copy whole directories from appropriate runs directories in case there is something more than piernik and problem.par
mkdir -p $RUN_TEST_DIR
cp -a runs/$(basename ${PROBLEM_NAME})_${FLAT_PROBLEM_NAME}/* $RUN_TEST_DIR

# Compile current versions of piernik
sed -i 's/-fcheck=all/& -fcheck=no-array-temps/' ${TEST_DIR}$OBJ"/Makefile"
make -j -C ${TEST_DIR}$OBJ > ${TEST_DIR}/make.stdout 2>&1

# Detect whether we have genuine riemann setup (which is unlikely to run in RTVD)
grep -qi '^ *solver_str *= *"riemann"' ${RUN_GOLD_DIR}/problem.par && RIEMANN=1 || RIEMANN=0
# ToDo set this up in the config file to make it easier to implement dual RTVD/Riemann gold tests

# Run the tests on current version of Piernik
if [ $RIEMANN == 1 ] ; then
    (
	cd $RUN_TEST_DIR
	eval $RUN_COMMAND ./${PIERNIK} > test_Riemann_stdout
    )
else
    RUN_TEST_DIR2=${TEST_DIR}runs/$(basename ${PROBLEM_NAME})"_Riemann"/
    cp -a $RUN_TEST_DIR $RUN_TEST_DIR2
    (
	cd $RUN_TEST_DIR
	eval $RUN_COMMAND ./${PIERNIK} "-n '&NUMERICAL_SETUP solver_str = \"RTVD\" /'"  > test_RTVD_stdout
    ) &
    (
	RIEM_ERR=test_Riemann_stderr
	cd $RUN_TEST_DIR2
	eval $RUN_COMMAND ./${PIERNIK} "-n '&NUMERICAL_SETUP solver_str = \"Riemann\" /'"  > test_Riemann_stdout 2> $RIEM_ERR
	[ -s $RIEM_ERR ] && ( echo ${PROBLEM_NAME}": Riemann failed " ; grep "Error"  $RIEM_ERR | sort | grep -vE '(meaningful|Following)' ) 1>&2
    )
fi

wait
# Here background jobs should be finished

./bin/gdf_distance ${RUN_GOLD_DIR}${OUTPUT} ${RUN_TEST_DIR}${OUTPUT} 2>&1 | tee $GOLD_LOG
grep 'Difference of datafield `' $GOLD_LOG |\
    sed 's/.*`\([^ ]*\)[^ ] *: \(.*\)/\1 \2/' |\
    awk '{a[$1]=$2} END {for (i in a) printf("log10(%s),", i); print ""; for (i in a) printf("%s,", (a[i]>0.)?(log(a[i])/log(10.)):1.) ; print ""}' |\
    sed 's/,$//' |\
    tee $GOLD_CSV

if [ $RIEMANN == 0 ] ; then
    # The tool gdf_distance distance is supposed to return values in [0..1] range
    # Map log10(0.) to 1. and failed Riemann to 2. (both impossible as a results of log10(gdf_distance)
    if [ -e ${RUN_TEST_DIR2}${OUTPUT} ] ; then
	./bin/gdf_distance ${RUN_TEST_DIR}${OUTPUT} ${RUN_TEST_DIR2}${OUTPUT} 2>&1 | tee $RIEM_LOG
	grep 'Difference of datafield `' $RIEM_LOG |\
	    sed 's/.*`\([^ ]*\)[^ ] *: \(.*\)/\1 \2/' |\
	    awk '{a[$1]=$2} END {for (i in a) printf("log10(%s),", i); print ""; for (i in a) printf("%s,", (a[i]>0.)?(log(a[i])/log(10.)):1.) ; print ""}' |\
	    sed 's/,$//' |\
	    tee $RIEM_CSV
    else
	(
	    echo "Riemann failed"
	    echo 2.
	) | tee $RIEM_CSV
    fi
fi

[ $HAS_KEEPPAR == 1 ] && echo " --keeppar" >> .setuprc

# Fail if gold distance is not 0.
[ $( ( grep "^Total difference between" $GOLD_LOG || echo 1 ) | awk '{print $NF}' ) == 0 ] || exit 1

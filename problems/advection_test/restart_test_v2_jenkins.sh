#!/bin/bash

RUNDIR=runs/advection_test

rm -rf ${RUNDIR} compare.log ${RUNDIR}/compare.log 2> /dev/null

./setup advection_test -p problem.par.restart_test_v2 | grep -v '[FC]'

(
    cd $RUNDIR

    echo "run_id = ts1"
    echo "Start:   t = 0.0, nproc = 1"
    mpiexec -n 1 ./piernik > ts1.out
    echo "Restart: t = 1.0, nproc = 1"
    mpiexec -n 1 ./piernik -n '$END_CONTROL tend = 2.0 /' >> ts1.out
    echo "Finish:  t = 2.0"
    echo
    echo "run_id = ts2"
    echo "Start:   t = 0.0, nproc = 5"
    mpiexec -n 5 ./piernik -n '$OUTPUT_CONTROL run_id = "ts2" /' > ts2.out
    echo "Restart: t = 1.0, nproc = 3"
    mpiexec -n 3 ./piernik -n '$END_CONTROL tend = 2.0 /' -n '$OUTPUT_CONTROL run_id = "ts2" /' >> ts2.out
    echo "Finish:  t = 2.0"
) 2> ${RUNDIR}/ts.err

[ ! -z $YT ] && source $YT
./bin/gdf_distance ${RUNDIR}/moving_pulse_ts{1,2}_0002.h5 | tee ${RUNDIR}/compare.log
[ $( grep "^Total difference between" ${RUNDIR}/compare.log | awk '{print $NF}' ) == 0 ] || exit 1

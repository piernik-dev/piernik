#!/bin/bash
set -e

export PYTHONPATH=/mnt/xarth/local/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/mnt/xarth/local/bin:$PATH

rm -rf runs/advection_test

echo "Compilation started"
./setup advection_test -p problem.par.restart_test_v2 -c gnu47 &> /dev/null
echo "Compilation finished"

cp problems/advection_test/*py runs/advection_test

cd runs/advection_test
echo "Starting run for tend = 1.0 and nproc = 1 run_id = ts1"
mpiexec -n 1 ./piernik > /dev/null
echo "Finished run for tend = 1.0 and nproc = 1 run_id = ts1"

echo "Starting run for tend = 2.0 and nproc = 1 run_id = ts1"
mpiexec -n 1 ./piernik -n '$END_CONTROL tend = 2.0 /' > /dev/null
echo "Finished run for tend = 2.0 and nproc = 1 run_id = ts1"

echo "Starting run for tend = 1.0 and nproc = 5 run_id = ts2"
mpiexec -n 5 ./piernik -n '$OUTPUT_CONTROL run_id = "ts2" /' > /dev/null
echo "Finished run for tend = 1.0 and nproc = 5 run_id = ts2"

echo "Starting run for tend = 2.0 and nproc = 3 run_id = ts2"
mpiexec -n 3 ./piernik -n '$END_CONTROL tend = 2.0 /' \
   -n '$OUTPUT_CONTROL run_id = "ts2" /' > /dev/null
echo "Finished run for tend = 2.0 and nproc = 3 run_id = ts2"

echo "Starting testsuite"
python advection_test.py
echo "Testsuite finished"

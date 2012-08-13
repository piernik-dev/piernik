#!/usr/bin/python
from mpi4py import MPI
from array import array
#import numpy as np
import os
import sys
try:
    import argparse
    have_argparse = True
except ImportError:
    from optparse import OptionParser
    have_argparse = False


def get_np():
    temp = os.getenv('PBS_NP')
    if temp is None:
        print "It seems we're not running in PBS shell"
        temp = -1
    return int(temp)

if have_argparse:
    parser = argparse.ArgumentParser(description='Run PIERNIK via MPI_Spawn')
    parser.add_argument('-e', '--exe', type=str, default='obj/piernik',
                       help='Path to executable')
    parser.add_argument('-a', '--args', type=str,
                       default='-p obj/ -w obj/',
                       help='arguments passed to PIERNIK')
    run_args = parser.parse_args()
else:
    parser = OptionParser()
    parser.add_option('-e', '--exe', type=str, default='obj/piernik',
                     help='Path to executable')
    parser.add_option('-a', '--args', type=str,
                     default='-p obj/ -w obj/',
                     help='arguments passed to PIERNIK')
    run_args = parser.parse_args()
    (run_args, add_args) = parser.parse_args()

num_procs = get_np()
if num_procs < 2:
    print "I need more than two processors to run"
    sys.exit(2)

signal = array('i', [0])

piernik = MPI.COMM_SELF.Spawn(run_args.exe, maxprocs=num_procs - 1,
                             args=run_args.args.split())
status = MPI.Status()
while True:
    piernik.Recv([signal, MPI.INT],
         MPI.ANY_SOURCE, MPI.ANY_TAG,
         status)
    print signal
    if signal[0] == 1:
        print "Python script: TSL was updated"
        print "GRID magic can happen here \o/"
#  if status.Get_tag() == 1:
    if signal[0] <= 0:
        break

piernik.Disconnect()
print "done in python too"

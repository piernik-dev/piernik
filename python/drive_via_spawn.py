#!/usr/bin/python
from mpi4py import MPI
from array import array
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


def parse_mpisignals(fn='../src/base/mpi_signals.F90'):
    import re
    enum = re.compile("enumerator", re.IGNORECASE)
    equal = re.compile("=")

    f = file(fn, 'r')
    lines = [line.split('::')[-1].strip()
            for line in  filter(enum.search, f.readlines())]
    f.close()

    mpi_signals = {}
    i = 0
    for line in lines:
        if equal.search(line):
            temp = line.split("=")
            name, i = temp[0], int(temp[-1])
        else:
            name = line
        mpi_signals[name[4:].lower().strip()] = i
        i += 1

    return mpi_signals


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
mpi_signals = parse_mpisignals('src/base/mpi_signals.F90')

piernik = MPI.COMM_SELF.Spawn(run_args.exe, maxprocs=num_procs - 1,
                             args=run_args.args.split())
status = MPI.Status()
while True:
    piernik.Recv([signal, MPI.INT],
         MPI.ANY_SOURCE, MPI.ANY_TAG,
         status)
    if not signal[0] in mpi_signals.values():
        print("I've got signal %i but I have no clue what to do with it :/" %
        signal[0])

    if signal[0] == mpi_signals['tsl_updated']:
        print "Python script: TSL was updated"
        print "GRID magic can happen here \o/"
    elif signal[0] == mpi_signals['clean_exit']:
        print "It seems that Piernik finished execution gracefully \o/"
#  if status.Get_tag() == 1:
    if signal[0] <= 0:
        break

piernik.Disconnect()

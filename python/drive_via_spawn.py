#!/usr/bin/python
from mpi4py import MPI
from array import array
import os
import sys
#import numpy as np

def get_np():
   temp = os.getenv('PBS_NP')
   if temp is None:
      print "It seems we're not running in PBS shell"
      temp = -1
   return int(temp)

num_procs = get_np()
if num_procs < 2:
   print "I need more than two processors to run"
   sys.exit(2)

executable = 'obj/piernik'
args = ['-p', 'runs/sedov']
signal  = array('i', [0])
#n = np.array([0])

piernik = MPI.COMM_SELF.Spawn(executable, maxprocs=num_procs-1, args=args)
status = MPI.Status()
while True:
#  piernik.Reduce(sendbuf=None,
#                 recvbuf=[n, MPI.INT],
#                 op=MPI.MAX, root=MPI.ROOT)
   piernik.Recv([signal, MPI.INT],
       MPI.ANY_SOURCE, MPI.ANY_TAG,
       status)
   print signal
   if signal[0] == 10:
      print "Python script: TSL was updated"
      print "GRID magic can happen here \o/"
#  if status.Get_tag() == 1:
   if signal[0] < 0:
      break

piernik.Disconnect()
print "done in python too"

#!/usr/bin/python
from mpi4py import MPI
from array import array

executable = 'obj/piernik'
args = ['-p', 'runs/sedov']

n  = array('i', [0])

piernik = MPI.COMM_SELF.Spawn(executable, maxprocs=2, args=args)
while n[0] >= 0:
   piernik.Reduce(sendbuf=None,
                  recvbuf=[n, MPI.INT],
                  op=MPI.MAX, root=MPI.ROOT)
   if n[0] == 10
      print "Python script: TSL was updated"
      print "GRID magic can happen here \o/"

piernik.Disconnect()
print "done in python too"

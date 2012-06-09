#!/usr/bin/python
import sys
import getopt
import os
import subprocess as sp
import shutil as sh
import xmlrpclib


class MakeTest(object):
   def __init__(self,test):

      self.initpath = os.getcwd()
      self.runpath  = os.path.join(self.initpath,'runs',test)
      self.test     = test
      os.chdir(self.runpath)
      retcode = sp.call(["mpiexec","./piernik"])
      if retcode != 0:
         sys.exit(retcode)
      self.runtest(test)
      os.chdir(self.initpath)

   def put_png(self):
      server = xmlrpclib.ServerProxy("http://piernik:p1ern1k@hum/piernik/login/xmlrpc")
      for file in os.listdir(self.runpath):
         if file.find('png') != -1:
            server.wiki.putAttachment(self.test+'/'+file, xmlrpclib.Binary(open(self.runpath+'/'+file).read()))

   def testJeans (self):
      sp.call(["gnuplot","verify.gpl"])
      self.put_png()

   def testMaclaurin (self):
      from maclaurin import Maclaurin_test
      Maclaurin_test(self.runpath+'/maclaurin_sph_0001.h5')
      self.put_png()

   def testSedov (self):
      print "test not implemented"

   def output(self):
      print self.initpath
      print self.runpath

   def runtest(self,test):
      tests = { "jeans": self.testJeans,
                "maclaurin": self.testMaclaurin,
                "sedov": self.testSedov}[test]()
      #tests.get(test)

def usage():
   print __doc__

def main(argv):
   try:
      opts, args = getopt.getopt(argv, "ht", ["help", "test="])
   except getopt.GetoptError:
      usage()
      sys.exit(2)
   for opt, arg in opts:
      if opt in ("-h", "--help"):
         usage()
         sys.exit()
      elif opt in ("-t", "--test"):
         test=arg

   # add piernik modules
   sys.path.append(sys.path[0]+'/python')
   t = MakeTest(test)
   t.output()

if __name__ == "__main__":
    main(sys.argv[1:])


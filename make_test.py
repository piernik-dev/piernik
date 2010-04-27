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
      os.chdir(self.runpath)
      sp.call(["./piernik"])
      self.runtest(test)
      os.chdir(self.initpath)

   def testJeans (self):
      sp.call(["gnuplot","verify.gpl"])
      server = xmlrpclib.ServerProxy("http://piernik:p1ern1k@hum/piernik/login/xmlrpc")
      for file in os.listdir(self.runpath):
         if file.find('png') != -1:
            server.wiki.putAttachment('jeans/'+file, xmlrpclib.Binary(open(file).read()))
   def testSedov (self):
      print "test not implemented"

   def output(self):
      print self.initpath
      print self.runpath

   def runtest(self,test):
      tests = { "jeans": self.testJeans,
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

   t = MakeTest(test)
   t.output()

if __name__ == "__main__":
    main(sys.argv[1:])


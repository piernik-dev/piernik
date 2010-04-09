#!/usr/bin/python

import subprocess as sp

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

b = bcolors()

def remove_binaries(files):
   list = []
   for file in files:
      if(not sp.Popen('file -bi '+file, stdout=sp.PIPE, shell=True, executable="/bin/bash").communicate()[0].startswith('text') ):
         print b.WARNING + "QA:  " + b.ENDC + file + " is a binary. I will not test it."
      else:
         list.append(file)
   return list

def qa_checks(files,options):
   print b.OKBLUE + '"I am the purifier, the light that clears all shadows." - seal of cleansing inscription' + b.ENDC
   files = remove_binaries(files)
   qa_trailing_spaces(files,options)
   qa_nonconforming_tabs(files,options)

def qa_nonconforming_tabs(files,options):
   print b.OKGREEN + "QA: " + b.ENDC + "Checking for tabs"
   wrong_files = []
   for file in files:
      if (len(sp.Popen('grep "	" '+file, shell=True, executable="/bin/bash",stdout=sp.PIPE).communicate()[0])):
         print b.FAIL + "QA:  " + b.ENDC + "non conforming tab detected in " + file
         wrong_files.append(file)
   if (len(wrong_files)):
      print b.FAIL + "Error detected! " + b.ENDC + "I will not let you commit unless you force me!!!"
      if options.force:
         print "Damn! You are determined to make low quality commit :-("
      else:
         exit()

def qa_trailing_spaces(files,options):
   print b.OKGREEN + "QA: " + b.ENDC + "Removing trailing whitespaces"
   for file in files:
      sp.Popen('sed -i -e "s/\s\+$//" '+file, shell=True, executable="/bin/bash")
      if options.verbose:
         print b.OKGREEN + "QA: " + b.ENDC + " done cleansing in " + file.rstrip('\n')

if __name__ == "__main__":
   from optparse import OptionParser
   usage = "usage: %prog [options] FILES"
   parser = OptionParser(usage=usage)
   parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="make lots of noise [default]")
   parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose",
                  help="be vewwy quiet (I'm hunting wabbits)")
   parser.add_option("-f", "--force",
                  action="store_true", dest="force",
                  help="commit despite errors (It will be logged)")
   (options, args) = parser.parse_args()
   if len(args) < 1:
      parser.error("incorrect number of arguments")
   qa_checks(args,options)

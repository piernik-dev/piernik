#!/usr/bin/python

import re, sys
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
         print b.WARNING + "QA:  " + b.ENDC + file + " is not a text file. I will not test it."
      else:
         list.append(file)
   return list

def select_sources(files):
   list = []
   for file in files:
      if( re.search("F90$",file)):
         list.append(file)
   return list

def qa_checks(files,options):
   print b.OKBLUE + '"I am the purifier, the light that clears all shadows." - seal of cleansing inscription' + b.ENDC
   runpath = sys.argv[0].split("qa.py")[0]
   files = remove_binaries(files)
   f90files = select_sources(files)
   qa_trailing_spaces(files,options)
   qa_nasty_spaces(f90files,options)
   nt = qa_labels(f90files,options)
   nt += qa_nonconforming_tabs(files,options)
   nt += qa_implicit_saves(f90files,options,runpath)
   wc = qa_depreciated_syntax(f90files,options)
   wc += qa_crude_write(f90files,options)
   wc += qa_magic_integers(f90files,options)
   if (wc):
      s =  b.WARNING + "%i warnings detected. " % wc + b.ENDC + "Do you wish to proceed? (y/N) "
      if( raw_input(s) != 'y' ):
         print "See you later!"
         exit()
   else:
      print b.OKGREEN + "No warnings detected. " + b.ENDC + "If everyone were like you, I'd be out of business!"
   if (len(nt)):
      print b.FAIL + "%i error(s) detected! " % len(nt) + b.ENDC + "I will not let you commit unless you force me!!!"
      if options.force:
         print "Damn! You are determined to make low quality commit :-("
      else:
         exit()
   else:
      print b.OKGREEN + "Yay! No errors!!! " + b.ENDC

def qa_implicit_saves(files,options,runpath):
   print b.OKGREEN + "QA: " + b.ENDC + "Checking for implicit saves"
   wrong_files = []
   dirname = runpath.split("svnci.py")[0]
   for file in files:
      p = sp.Popen(dirname+'implicit_save.sh '+file, shell=True, executable="/bin/bash", stdout=sp.PIPE).communicate()[0]
      if (len(p)):
         print p.rstrip()
         wrong_files.append(file)
   return wrong_files

def qa_labels(files,options):
   print b.OKGREEN + "QA: " + b.ENDC + "Checking for labels"
   wrong_files = []
   for file in files:
      if (len(sp.Popen('grep "^[0-9]" '+file, shell=True, executable="/bin/bash",stdout=sp.PIPE).communicate()[0])):
         print b.FAIL + "QA:  " + b.ENDC + "label detected in " + file
         wrong_files.append(file)
   if(len(wrong_files) > 0):
      print b.FAIL + "Seriously!? Don't tell me your using \"goto\" too. It's XXI century..." + b.ENDC
   return wrong_files

def qa_crude_write(files,options):
   print b.OKGREEN + "QA: " + b.ENDC + "Checking for crude writes to stdout (suppress this by appending QA_WARN at the eol)"
   warning = 0
   for name in files:
      file = open(name)
      n = 1
      for line in file.readlines():
         if( re.search("write *\( *\*", line, flags=re.IGNORECASE) and not re.search("QA_WARN",line)):
            print b.WARNING + "!! crude write  " + b.ENDC + name + " @ L%i => " % n + line.strip()
            warning += 1
         n += 1
      file.close()
   return warning

#ToDo: look for dimension([1-9]
#ToDo: also warn about strings shorter than 10 characters. These also can make funny errors.
def qa_magic_integers(files,options):
   print b.OKGREEN + "QA: " + b.ENDC + "Checking for hardcoded integers"
   warning = 0
   for name in files:
      file = open(name)
      n = 1
      for line in file.readlines():
         if( re.search("\(len=[1-9][0-9]",line, flags=re.IGNORECASE) and not re.search("QA_WARN",line)):
            print b.WARNING + "!! magic integer  " + b.ENDC + name + " @ L%i => " % n + line.strip()
            warning += 1
         n += 1
      file.close()
   return warning

def qa_nonconforming_tabs(files,options):
   print b.OKGREEN + "QA: " + b.ENDC + "Checking for tabs"
   wrong_files = []
   for file in files:
      if (len(sp.Popen('grep "\t" '+file, shell=True, executable="/bin/bash",stdout=sp.PIPE).communicate()[0])):
         print b.FAIL + "QA:  " + b.ENDC + "non conforming tab detected in " + file
         wrong_files.append(file)
   return wrong_files

def qa_nasty_spaces(files,options):
   print b.OKGREEN + "QA: " + b.ENDC + "Removing spurious whitespaces"
   for file in files:
      sp.Popen('sed -i -e "s/end *if/endif/;s/end *do/enddo/;s/end *where/endwhere/;s/only *:/only:/;s/if(/if (/;s/where(/where (/;s/while(/while (/;s/forall(/forall (/;s/ case(/ case (/" '+file, shell=True, executable="/bin/bash")
      if options.verbose:
         print b.OKGREEN + "QA: " + b.ENDC + " done cleansing in " + file.rstrip('\n')

def qa_trailing_spaces(files,options):
   print b.OKGREEN + "QA: " + b.ENDC + "Removing trailing whitespaces"
   for file in files:
      sp.Popen('sed -i -e "s/\s\+$//" '+file, shell=True, executable="/bin/bash")
      if options.verbose:
         print b.OKGREEN + "QA: " + b.ENDC + " done cleansing in " + file.rstrip('\n')

def qa_depreciated_syntax(files,options):
   print b.OKGREEN + "QA: " + b.ENDC + "Checking for depreciated syntax"
   warning = 0
   for name in files:
      file = open(name)
      n = 1
      for line in file.readlines():
         if( re.search("^\s{1,12}(?i)(?:real(?:\s|,)|integer(?:\s|,)|logical(?:\s|,|\()|character(?:\s|,))(?!.*::)",line, flags=re.IGNORECASE) and
            not re.search("(?i)\sfunction\s",line, flags=re.IGNORECASE)):
               print b.WARNING + "!! lacking ::   " + b.ENDC + name + " @ L%i => " % n + line.strip()
               warning=warning+1
         if( re.search("^\s{1,12}(?i)use[\s](?!.*only)",line, flags=re.IGNORECASE)):
               print b.WARNING + "!! greedy use   " + b.ENDC + name + " @ L%i => " % n + line.strip()
               warning=warning+1
         if( re.search("^\s{1,12}(?i)character(?![(])",line, flags=re.IGNORECASE)):
               print b.WARNING + "!! wrong syntax " + b.ENDC + name + " @ L%i => " % n + line.strip()
               warning=warning+1
         n = n + 1
      file.close()
   return warning

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

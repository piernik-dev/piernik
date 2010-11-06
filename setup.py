#!/usr/bin/python

import os, re
import subprocess as sp
import numpy as np

columns = 120

is_f90       = re.compile("f90$", re.IGNORECASE)
not_svn_junk = re.compile(".*(?!pro).*", re.IGNORECASE)
test     = re.compile(r'pulled by',re.IGNORECASE).search
test2    = re.compile(r"\$Id",re.IGNORECASE).search
have_use = re.compile(r"^\s{0,9}use\s",re.IGNORECASE).search
have_inc = re.compile(r"^#include\s",re.IGNORECASE).search
cpp_junk = re.compile("(?!#define\s_)", re.IGNORECASE)

def striplist(l):
   return([x.strip() for x in l])

def strip_leading_path(l):
   return([x.rpartition('/')[2] for x in l])

def remove_suf(l):
   return([x.partition('.')[0] for x in l])

def pretty_format(fname,list,col):
   str = fname
   for item in list:
      if(len(str) + len(item) + 2 > int(col)):
         print str + "\\"
         str = "\t"
      str= str + item + " "
   print str

def pretty_format_suf(fname,list,suf,col):
   print fname
   str = "\t"
   for item in list:
      if(len(str) + len(item) + len(suf) + 2 > int(col)):
         print str + "\\"
         str = "\t"
      str= str + item + suf + " "
   print str
   print ""

class DirectoryWalker:
    # a forward iterator that traverses a directory tree

    def __init__(self, directory):
        self.stack = [directory]
        self.files = []
        self.index = 0

    def __getitem__(self, index):
        while 1:
            try:
                file = self.files[self.index]
                self.index = self.index + 1
            except IndexError:
                # pop next directory from stack
                self.directory = self.stack.pop()
                self.files = os.listdir(self.directory)
                self.index = 0
            else:
                # got a filename
                fullname = os.path.join(self.directory, file)
                if os.path.isdir(fullname) and not os.path.islink(fullname):
                    self.stack.append(fullname)
                return fullname

f90files = [f for f in DirectoryWalker('src') if(is_f90.search(f))]
for f in DirectoryWalker('problems/mcrtest'):         # BEWARE: testing on mcrtest
   if(is_f90.search(f)): f90files.append(f)

defines  = sp.Popen(["echo '#include \"piernik.def\"' > foo.f90 && cpp $cppflags -dM foo.f90 && rm foo*"], stdout=sp.PIPE, shell="/bin/bash").communicate()[0].rstrip().split("\n")
our_defs = [f.split(" ")[1] for f in filter(cpp_junk.match,defines)]
our_defs.append("ANY")

files = ['src/base/defines.c']
tags  = ['']   # BEWARE missing tag for defines.c
uses  = [[]]
incl  = ['']

for f in f90files:
   tag  = ""
   keys = []
   luse = []
   linc = []
   for line in file(f):
      if test2(line):
         tag  = line.strip()
      if test(line):
         keys = striplist(line.split(" ")[3:])
      if have_use(line):
         luse.append(line.split()[1].rstrip(","))
      if have_inc(line):
         linc.append( line.split('"')[1] )
# Logic here should be improved...
   if(len(keys) == 0  or (len(keys) == 1 and keys[0] in our_defs)):
      files.append(f)
      tags.append(tag)
      uses.append( np.unique(luse) )
      incl.append( np.unique(linc) )
   if(len(keys) == 3):
      if((keys[1] == "&&" and (keys[0] in our_defs and keys[2] in our_defs)) or
         (keys[1] == "||" and (keys[0] in our_defs or  keys[2] in our_defs))   ):
         files.append(f)
         tags.append(tag)
         uses.append( np.unique(luse) )
         incl.append( np.unique(linc) )

stripped_files  = strip_leading_path(files)

stripped_files.append("version.F90")   # adding version
incl.append('')
uses.append([])

files_to_build = remove_suf(stripped_files)

pretty_format_suf("SRCS = \\", stripped_files, '', columns)
pretty_format_suf("OBJS = \\", files_to_build, '.o', columns)

for i in range(0,len(files_to_build)):
   deps = files_to_build[i]+".o: "+stripped_files[i]+" "
   d = ""
   if(np.size(incl[i]) > 0): d += ' '.join(incl[i])+' '
   if(np.size(uses[i]) > 0): d += '.o '.join(set(uses[i]).intersection(files_to_build))+'.o'
   pretty_format(deps, d.split(), columns)

addons = []

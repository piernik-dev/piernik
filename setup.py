#!/usr/bin/python

import os, re, shutil, sys
import subprocess as sp
from optparse import OptionParser
try:
   import multiprocessing
   mp = True
except ImportError:
   print "No multiprocessing: The make command will be run in single thread unless you specified MAKEFLAGS += -j<n> in compiler setup."
   mp = False

columns = 90

is_f90       = re.compile("f90$", re.IGNORECASE)
is_header    = re.compile("h$", re.IGNORECASE)
not_svn_junk = re.compile(".*(?!pro).*", re.IGNORECASE)
test     = re.compile(r'pulled by',re.IGNORECASE).search
test2    = re.compile(r"\$Id",re.IGNORECASE).search
have_use = re.compile(r"^\s{0,9}use\s",re.IGNORECASE).search
have_inc = re.compile(r"^#include\s",re.IGNORECASE).search
have_mod = re.compile(r"^\s*module\s+(?!procedure)",re.IGNORECASE).search
cpp_junk = re.compile("(?!#define\s_)", re.IGNORECASE)

desc='''
EXAMPLE:
> cd obj
> ./newcompiler <settingsname>
> make
then copy files to your run directory (optional), e.g.
> cp {piernik,problem.par} ../runs/<problem>
> cd ../runs/<problem>

to run PIERNIK:
edit problem.par as appropriate, e.g.
* add var names for visualisation => var(<number>)='<name>'
* change domain dimensions/resolution => DOMAIN_SIZES
* change domain divisions for parallel processing => MPI_BLOCKS
* change frequency of data dumps => dt_* entries
* etc.
execute
> ./piernik
or for <np> parallel processes
> mpirun -n <np> ./piernik

HEALTH WARNINGS:
* the contents of \'./obj\' and \'./runs/<problem>\' are overwritten
  each time setup <problem>\' is run, unless you specify -obj <postfix>
  in which case the contents of runs/<problem>_<postfix> will be only updated
* the def file \'piernik.def\' is copied only for reference, to change flags
  with which the source is compiled edit \'./problems/<problem>/piernik.def\'
* by default PIERNIK will read the configuration file \'problem.par\' from the
working directory, to use alternative configurations execute
\'./piernik <directory with an alternative problem.par>\'
Enjoy your work with the Piernik Code!
'''

head_block1='''LIBS +=${STATIC} -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz ${DYNAMIC}

RM ?= /bin/rm
MV ?= /bin/mv
ECHO ?= /bin/echo
ifdef CHECK_MAGIC
\tRM  = /bin/true
\tMV  = /bin/true
\tF90 = /bin/true
endif
ifndef PRECOMP
\tPRECOMP=cpp
endif
ifeq ("$(SILENT)","1")
MAKEFLAGS += -s
define ECHO_FC
@$(ECHO) [FC] $<
endef
define ECHO_CC
@$(ECHO) [CC] $<
endef
else
define ECHO_FC
endef
define ECHO_CC
endef
endif

.PHONY: print_fc

all: date print_fc $(PROG)

print_fc:
ifeq ("$(SILENT)","1")
\t@$(ECHO) FC = $(F90) $(CPPFLAGS) $(F90FLAGS) -c
\t@$(ECHO) CC = $(CC) $(CPPFLAGS) $(CFLAGS) -c
endif

$(PROG): $(OBJS)
\t@$(ECHO) $(F90) $(LDFLAGS) -o $@ '*.o' $(LIBS)
\t@$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
\t@AO1=`mktemp _ao_XXXXXX`;\\
\tAO2=`mktemp _ao_XXXXXX`;\\
\t$(ECHO) $(OBJS) | tr ' ' '\\n' | sort > $AO1;\\
\t$(ECHO) *.o     | tr ' ' '\\n' | sort > $AO2;\\
\tif [ `join -v 2 $AO1 $AO2 | wc -l` -gt 0 ] ; then\\
\t\t$(ECHO) -n "WARNING: unused object files: ";\\
\t\tjoin -v 2 $AO1 $AO2 | tr '\\n' ' ';\\
\t\t$(ECHO);\\
\tfi;\\
\t$(RM) $AO1 $AO2

date:
\t@if [ ! -f env.dat ]; then\\
'''

head_block2='''
\t\tsed -n '/Id:/p' *.h *.c *.F90 | column -t >> env.dat; \\
\t\tsed -e '/^$$/ d' -e "/^\// d" piernik.def >> env.dat; \\
\tfi;

version.F90: date
\t@if [ -e version.F90 ]; then unlink version.F90; fi;
\t@( $(ECHO) -e "module version\\n   implicit none\\n   public\\n"; \\
\twc -l env.dat | awk '{print "   integer, parameter :: nenv = "$$1"+0"}'; \\
\t$(ECHO) -e "   character(len=128), dimension(nenv) :: env\\ncontains\\n   subroutine init_version\\n\t\timplicit none"; \\
\tawk '{printf("\\t\\t env(%i) = \\"%s\\"\\n",NR,$$0)}' env.dat; \\
\t$(ECHO) -e "    end subroutine init_version\\nend module version" ) > version.F90; \\
\t$(ECHO) 'generated version.F90';

clean:
\t$(RM) $(PROG) $(OBJS) *.mod

clean-run:
\t$(RM) *.bck *~ *.h5 *.res *.log *.tsl *.out *.tmp core*

clean-all:
\t$(RM) $(PROG) $(OBJS) *.mod *.bck *~ *.h5 *.res *.log *.tsl *.out *.tmp core* *.f *.dbg

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
\t$(ECHO_FC)
ifdef USE_GNUCPP
\t@$(PRECOMP) $(CPPFLAGS) $< $(patsubst %.F90,%_.f90,$<) 2>&1 | sed -n "/!warning: missing terminating \' character/p" 
\t$(F90) $(F90FLAGS) -c $(patsubst %.F90,%_.f90,$<) && ( $(RM) $(patsubst %.F90,%_.f90,$<); $(MV) $(patsubst %.F90,%_.o,$<) $(patsubst %.F90,%.o,$<) )
else
\t$(F90) $(CPPFLAGS) $(F90FLAGS) -c $<
endif

.c.o:
\t$(ECHO_CC)
\t$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

%.o : %.mod

'''

def striplist(l):
   return([x.strip() for x in l])

def strip_leading_path(l):
   return([x.rpartition('/')[2] for x in l])

def remove_suf(l):
   return([x.partition('.')[0] for x in l])

def pretty_format(fname,list,col):
   out = ""
   sl = True
   str = fname
   for item in list:
      if(len(str) + len(item) + 2 > int(col)):
         out += str + "\\\n"
         str = "\t"
         sl = False
      str= str + item + " "
   if(str != "\t"): out += str
   if(sl):
      return str.rstrip("\\\n")+"\n"
   else:
      return out.rstrip("\\\n")+"\n"

def pretty_format_suf(fname,list,suf,col):
   out = fname+'\n'
   str = "\t"
   for item in list:
      if(len(str) + len(item) + len(suf) + 2 > int(col)):
         out += str + "\\\n"
         str = "\t"
      str= str + item + suf + " "
   if(str != "\t"): out += str
   return out.rstrip("\\\n")+"\n"

def get_stdout(cmd):
   nul_f = open(os.devnull, 'w')
   process = sp.Popen([cmd], stdout=sp.PIPE, shell="/bin/bash", stderr = nul_f)
   nul_f.close()
   return process.communicate()[0]

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

#print desc
usage = "usage: %prog [options] FILES"
parser = OptionParser(usage=usage)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
   help="try to confuse the user with some diagnostics ;-)")
parser.add_option("-q", "--laconic", action="store_true", dest="laconic", default=False,
   help="compress stdout from make process")
parser.add_option("-n", "--nocompile", action="store_true", dest="nocompile", default=False,
   help='''Create object directory, check for circular dependencies, but do not compile or prepare run directory.
In combination with --copy will prepare portable object directory.''')
parser.add_option("--problems", dest="show_problems", action="store_true",
   help="print available problems and exit")
parser.add_option("-u","--units", dest="show_units", action="store_true",
   help="print available units and exit")
parser.add_option("--copy", dest="hard_copy", action="store_true",
   help="hard-copy source files instead of linking them")
parser.add_option("-l","--linkexe", dest="link_exe", action="store_true",
   help="do not copy obj/piernik to runs/<problem> but link it instead")
parser.add_option("-p", "--param", dest="param",
   help="use FILE instead problem.par", metavar="FILE", default="problem.par")
parser.add_option("-d", "--define", dest="cppflags",
   help="add precompiler directives, use comma-separated list", metavar="CPPFLAGS")
parser.add_option("-c", "--compiler", dest="compiler", default="default.in",
   help="choose specified config from compilers directory", metavar="FILE")
parser.add_option("-o", "--obj", dest="objdir", metavar="POSTFIX", default='',
   help="use obj_POSTFIX directory instead of obj/ and runs/<problem>_POSTFIX rather than runs/<problem>")

all_args = []
try:
   for line in file(".setuprc"):
      all_args += line.split()
except IOError:
   pass
all_args += sys.argv[1:]
(options, args) = parser.parse_args(all_args)

if(options.verbose):
   print "Setup options:"
   print options
   print "Setup arguments:"
   print args

if(options.show_problems):
   print get_stdout("cat problems/*/info")
   exit()

if(options.show_units):
   print get_stdout("grep uses ./src/base/constants.F90")
   exit()

if (len(args) < 1):
   parser.error("incorrect number of arguments")

# set problem dir
probdir = 'problems/'+args[0]+'/'

# parse cppflags
if(options.cppflags):
   cppflags = '-D' + ' -D'.join(options.cppflags.split(","))
else:
   cppflags = ""

# parse compiler
if(not re.search('\.in$',options.compiler)):
   compiler = options.compiler + '.in'
else:
   compiler = options.compiler


objdir = 'obj'
rundir = 'runs/'+args[0]
if(len(options.objdir)>0):
   objdir += '_'+options.objdir
   rundir += '_'+options.objdir+'/'
else:
   rundir += '/'

if(os.path.isdir(objdir)): shutil.rmtree(objdir)
os.mkdir(objdir)

f90files = []
allfiles = []
for f in DirectoryWalker('src'):
   if(is_f90.search(f)): f90files.append(f)
   if(is_header.search(f)): allfiles.append(f)

for f in DirectoryWalker(probdir):
   if(is_f90.search(f)): f90files.append(f)

allfiles.append(probdir+"piernik.def")
allfiles.append(probdir+options.param)

cmd = "echo '#include \"%spiernik.h\"' > foo.f90 && cpp %s -dM -I%s foo.f90 && rm foo*" % ('src/base/', cppflags, probdir)
defines  = sp.Popen([cmd], stdout=sp.PIPE, shell=True).communicate()[0].rstrip().split("\n")
if(options.verbose):
    print cmd
    print "Defined symbols:"
    for defin in defines: print defin

our_defs = [f.split(" ")[1] for f in filter(cpp_junk.match,defines)]
our_defs.append("ANY")
if(options.verbose):
   print "our_defs:"
   print our_defs

files = ['src/base/defines.c']
tags  = ['']   # BEWARE missing tag for defines.c
uses  = [[]]
incl  = ['']
module = dict()

for f in f90files:
   keys_logic1 = False
   keys_logic2 = False
   tag  = ""
   keys = []
   luse = []
   linc = []
   for line in open(f):       # Scan original files
      if test(line):
         keys = striplist(line.split(" ")[3:])
      if have_inc(line):
         linc.append(line.split('"')[1])

   keys_logic1 = len(keys) == 0  or (len(keys) == 1 and keys[0] in our_defs)
   if(len(keys) == 3):
      keys_logic2 = (keys[1] == "&&" and (keys[0] in our_defs and keys[2] in our_defs)) or (keys[1] == "||" and (keys[0] in our_defs or  keys[2] in our_defs))

   if(keys_logic1 or keys_logic2):
      cmd =  "cpp %s -I%s -I%s %s" % ( cppflags, probdir, 'src/base', f)
      for line in get_stdout(cmd).split('\n'):    # Scan preprocessed files
         if test2(line):
            tag  = line.strip()
         if have_use(line):
            luse.append(line.split()[1].rstrip(","))
         if have_mod(line):
            module.setdefault(line.split()[1], remove_suf(strip_leading_path([f]))[0])
      files.append(f)
      tags.append(tag)
      uses.append( list(set(luse) ) )
      incl.append( list(set(linc) ) )

#for i in iter(module):
#   if module[i] != i:
#      print "File",module[i]+".F90 contains an alien module", i

allfiles.extend(files)

for f in allfiles:
   if(options.hard_copy):
      shutil.copy(f,objdir)
   else:
      os.symlink('../'+f,objdir+'/'+strip_leading_path([f])[0])

makefile_head = open('compilers/'+compiler,'r').readlines()
m = open(objdir+'/Makefile', 'w')

stripped_files  = strip_leading_path(files)

stripped_files.append("version.F90")   # adding version
incl.append('')
uses.append([])
module.setdefault('version', 'version')

files_to_build = remove_suf(stripped_files)

for f in makefile_head: m.write(f)
m.write( pretty_format_suf("SRCS = \\", stripped_files, '', columns)   )
m.write( pretty_format_suf("OBJS = \\", files_to_build, '.o', columns) )
m.write( "\nCPPFLAGS += %s\n" % cppflags )
if( "PGPLOT" in our_defs ): m.write("LIBS += -lpgplot\n")
if( "SHEAR" in our_defs or "MULTIGRID" in our_defs ): m.write("LIBS += `pkg-config --libs fftw3`\n")
if( "POISSON_FFT" in our_defs): m.write("LIBS += `pkg-config --libs fftw3` `pkg-config --libs lapack`\n")
if( options.laconic ):
   m.write("SILENT = 1\n\n")
else:
   m.write("SILENT = 0\n\n")
m.write(head_block1)
m.write("\t\t$(ECHO) \"%s\" > env.dat; \\" % (sys.argv[0]+ " " + " ".join(all_args)))
m.write(head_block2)

for i in range(0,len(files_to_build)):
   deps = files_to_build[i]+".o: "+stripped_files[i]+" "
   d = ""
   if(len(incl[i]) > 0): d += ' '.join(incl[i])+' '
   for j in uses[i]:
      if j in module:
         if module[j]+'.F90' in stripped_files:
            d += module[j]+'.o '
         else:
            print "Warning: ",module[j]+'.F90',"referenced in",stripped_files[i],"not found!"
      else:
         print "Warning: module",j," from file ",stripped_files[i],"not found!"
   m.write( pretty_format(deps, d.split(), columns) )

m.close()

fatal_problem = False

if (not options.nocompile):
   makejobs = ""
   if (mp):
      makejobs = "-j%i" % multiprocessing.cpu_count()
   makecmd = "make %s -C %s" % ( makejobs, objdir)
   if( sp.call([makecmd], shell=True) != 0):
      print '\033[91m' + "It appears that '%s' crashed. Cannot continue." % makecmd + '\033[0m'
      exit()

try: os.makedirs(rundir)
except OSError:
   print '\033[93m' + "Found old run." + '\033[0m' + " Making copy of old 'problem.par'."
   try:
      if(os.path.isfile(rundir+'problem.par')): shutil.move(rundir+'problem.par',rundir+'problem.par.old')
   except (IOError, OSError ): print '\033[91m' + "Problem with copying 'problem.par' to 'problem.par.old'." + '\033[0m'
   try:
      if(os.path.isfile(rundir+'piernik')): os.remove(rundir+'piernik')
   except (IOError, OSError ): print '\033[91m' + "Problem with removing old 'piernik' executable from '%s'." % rundir.rstrip('/') + '\033[0m'
   try:
      if(os.path.isfile(rundir+'piernik.def')): os.remove(rundir+'piernik.def')
   except (IOError, OSError ): print '\033[91m' + "Problem with removing old 'piernik.def' from '%s'." % rundir.rstrip('/') + '\033[0m'

if (not options.nocompile):
   if(options.link_exe):
      try: os.symlink("../../"+objdir+"/piernik", rundir+'piernik')
      except (IOError, OSError ): print '\033[91m' + "Symlinking 'piernik' failed." + '\033[0m'; fatal_problem = True
   else:
      try: shutil.copy(objdir+"/piernik", rundir+'piernik')
      except (IOError, OSError, shutil.Error ): print '\033[91m' + "Copying 'piernik' failed." + '\033[0m'; fatal_problem = True

try: shutil.copy(objdir+"/"+options.param, rundir+'problem.par')
except IOError: print '\033[91m' + "Failed to copy 'problem.par' to '%s'." % rundir.rstrip('/') + '\033[0m';  fatal_problem = True
try: shutil.copy(objdir+"/piernik.def", rundir+'piernik.def')
except IOError: print '\033[91m' + "Failed to copy 'piernik.def' to '%s'." % rundir.rstrip('/') + '\033[0m'

if (options.nocompile):
   print '\033[93m' + "Compilation of '%s' skipped on request." % args[0] + '\033[0m' + " You may want to run 'make -C %s' before running the Piernik code." % objdir
   makejobs = ""
   if (mp):
      makejobs = "-j%i" % multiprocessing.cpu_count()
   makecmd = "LC_ALL=C make %s -C %s CHECK_MAGIC=yes piernik" % ( makejobs, objdir)
   output = sp.Popen(makecmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE).communicate()
   if re.search(r"Circular", output[1]):
      print '\033[91m' + "Circular dependencies foud in '%s'." % objdir + '\033[0m'
else:
   if (fatal_problem): print '\033[93m' + "'%s' compiled, but '%s' may not be ready to run." % (args[0], rundir.rstrip('/') ) + '\033[0m'
   else: print '\033[92m' + "'%s' ready in '%s'." % (args[0], rundir.rstrip('/') ) + '\033[0m'


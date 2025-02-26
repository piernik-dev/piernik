#!/usr/bin/env python3

import os
import re
import shutil
import sys
from DirWalk import DirectoryWalker
import subprocess as sp
import tempfile
from optparse import OptionParser

try:
    import multiprocessing
    mp = True
except ImportError:
    print("No multiprocessing: The make command will be run in single thread \
unless you specified MAKEFLAGS += -j<n> in compiler setup.")
    mp = False
try:
    import pickle
    pickle_avail = True
except ImportError:
    pickle_avail = False

columns = 90

is_f90 = re.compile(r"\.f90$", re.IGNORECASE)
is_header = re.compile(r"\.h$", re.IGNORECASE)
test = re.compile(r'pulled by', re.IGNORECASE).search
overriding = re.compile(r'overrides', re.IGNORECASE).search
have_use = re.compile(r"^\s*use\s", re.IGNORECASE).search
have_inc = re.compile(r"^#include\s", re.IGNORECASE).search
have_mod = re.compile(r"^\s*module\s+(?!procedure)", re.IGNORECASE).search
cpp_junk = re.compile(r"(?!#define\s_)", re.IGNORECASE)


head_block1 = '''ifneq (,$(findstring h5pfc, $(F90)))
LIBS += ${STATIC} -lz ${DYNAMIC}
else
LIBS += ${STATIC} -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz ${DYNAMIC}
endif

RM ?= /bin/rm
MV ?= /bin/mv
ECHO ?= /bin/echo
ifdef CHECK_MAGIC
\tMV  = /bin/true
\tF90 = /bin/true
endif
ifndef PRECOMP
\tPRECOMP=cpp
endif
ifeq ("$(SILENT)","1")
MAKEFLAGS += -s
ifdef PNAME
override PNAME+=""
endif
define ECHO_FC
@$(ECHO) [$(PNAME)FC] $<
endef
define ECHO_CC
@$(ECHO) [$(PNAME)CC] $<
endef
else
define ECHO_FC
endef
define ECHO_CC
endef
endif

CPPFLAGS := $(CPPFLAGS) $(shell $(F90) -c $(CPPFLAGS) $(F90FLAGS) ../compilers/tests/mpi_allgatherv_bug.F90 && $(F90) $(LDFLAGS) -o mpi_allgatherv_bug mpi_allgatherv_bug.o $(LIBS) && mpirun -np 1 ./mpi_allgatherv_bug || echo -DFORBID_F08)
CPPFLAGS := $(CPPFLAGS) $(shell $(F90) $(CPPFLAGS) $(F90FLAGS) ../compilers/tests/mpi_f08.F90 2> /dev/null && echo -DMPIF08 || echo -DNO_MPIF08_AVAILABLE)
CPPFLAGS := $(CPPFLAGS) $(shell $(F90) $(CPPFLAGS) $(F90FLAGS) ../compilers/tests/mpi.F90 2> /dev/null || echo -DNO_ALL_MPI_FUNCTIONS_AVAILABLE)
CPPFLAGS := $(CPPFLAGS) $(shell $(F90) $(CPPFLAGS) $(F90FLAGS) ../compilers/tests/F2018.F90 2> /dev/null || echo -DNO_F2018)

all: env.dat print_setup $(PROG)

check_mpi:
\t@$(F90) $(CPPFLAGS) $(F90FLAGS) ../compilers/tests/mpi_f08.F90
\t@$(F90) $(CPPFLAGS) $(F90FLAGS) ../compilers/tests/mpi.F90

$(PROG): $(OBJS) check_mpi
ifeq ("$(SILENT)","1")
\t@$(ECHO) $(PNAME)FC = $(F90) $(CPPFLAGS) $(F90FLAGS) -c
\t@$(ECHO) $(PNAME)CC = $(CC) $(CPPFLAGS) $(CFLAGS) -c
endif
\t@$(ECHO) $(F90) $(LDFLAGS) -o $@ '*.o' $(LIBS)
\t@$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
\t@touch mpi_f08.o mpi.o a.out mpi_allgatherv_bug mpi_allgatherv_bug.o F2018.o
\t@$(RM) mpi_f08.o mpi.o a.out mpi_allgatherv_bug mpi_allgatherv_bug.o F2018.o
\t@AO1=`mktemp _ao_XXXXXX`;\\
\tAO2=`mktemp _ao_XXXXXX`;\\
\t$(ECHO) $(OBJS) | tr ' ' '\\n' | sort > $$AO1;\\
\t$(ECHO) *.o     | tr ' ' '\\n' | sort > $$AO2;\\
\tif [ `join -v 2 $$AO1 $$AO2 | wc -l` -gt 0 ] ; then\\
\t\t$(ECHO) -n "WARNING: unused object files: ";\\
\t\tjoin -v 2 $$AO1 $$AO2 | tr '\\n' ' ';\\
\t\t$(ECHO);\\
\tfi;\\
\t$(RM) $$AO1 $$AO2

print_setup:
\t$(shell tail -n 1 .setup.call)

env.dat: piernik.def *.h $(SRCS_V)
'''

head_block2 = '''
\tawk '{print}' piernik.def | sed -e '/^$$/ d' -e "/^\\// d" ) > env.dat
\t@$(ECHO) "Recent history:" >> env.dat
\t@git log -5 --decorate --graph 2> /dev/null | sed -e 's/"//g' >> env.dat

os_detect.F90:
\t@( OS=INVALID; \\
\tgrep -q LINUX constants.F90 && gcc -dM -E - </dev/null | grep -q  __linux__ && OS=LINUX; \\
\tgrep -q APPLE constants.F90 && gcc -dM -E - </dev/null | grep -q  __APPLE__ && OS=APPLE; \\
\t$(ECHO) "module os_detect"; \\
\t$(ECHO) "   implicit none"; \\
\t$(ECHO) "   public"; \\
\t$(ECHO) "contains"; \\
\t$(ECHO) "   integer function which_os()"; \\
\t$(ECHO) "      use constants, only: "$$OS; \\
\t$(ECHO) "      implicit none"; \\
\t$(ECHO) "      which_os = "$$OS; \\
\t$(ECHO) "   end function which_os"; \\
\t$(ECHO) "end module os_detect" ) > os_detect.F90; \\
\t$(ECHO) 'generated os_detect.F90'

version.F90: env.dat
\t@( $(ECHO) "module version"; \\
\t$(ECHO) "   implicit none"; \\
\t$(ECHO) "   public"; \\
\twc -l env.dat | awk '{print "   integer, parameter :: nenv = "$$1"+0"}'; \\
\tawk '{if (length($0)>s) s=length($0)} END {print "   character(len="s+10"), dimension(nenv) :: env\\ncontains\\n   subroutine init_version\\n      implicit none"}' env.dat; \\
\tawk '{printf("      env(%i) = \\"%s\\"\\n",NR,$$0)}' env.dat; \\
\t$(ECHO) "   end subroutine init_version"; \\
\t$(ECHO) "end module version" ) > version_.F90; \\
\tif [ -e version.F90 ]; then diff version_.F90 version.F90 > /dev/null || unlink version.F90 ; fi; \\
\tif [ ! -e version.F90 ]; then mv version_.F90 version.F90 ; $(ECHO) 'generated version.F90'; fi

clean:
\t$(RM) $(PROG) $(OBJS) *.mod

clean-run:
\t$(RM) *.bck *~ *.h5 *.res *.log *.tsl *.out *.tmp core*

clean-all:
\t$(RM) $(PROG) $(OBJS) *.mod *.bck *~ *.h5 *.res *.log *.tsl *.out *.tmp core* *.f *.dbg

dep.png: dep.dot
\tdot dep.dot -T png -o dep.png -Gnodesep=0.05 -Granksep=0.3 -Gnslimit=10 -Nshape=box -Nfontsize=10 -Nheight=0 -Nwidth=0

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
\t$(ECHO_FC)
\t$(F90) $(CPPFLAGS) $(F90FLAGS) -c $<

.c.o:
\t$(ECHO_CC)
\t$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

%.o : %.mod

'''


def striplist(ll):
    return ([x.strip() for x in ll])


def strip_leading_path(ll):
    return ([x.split('/')[-1] for x in ll])
#  return ([x.rpartition('/')[2] for x in ll])


def remove_suf(ll):
    return ([x.split('.')[0] for x in ll])
#  return ([x.partition('.')[0] for x in ll])


def pretty_format(fname, list, col):
    out = ""
    sl = True
    str = fname
    for item in list:
        if (len(str) + len(item) + 2 > int(col)):
            out += str + "\\\n"
            str = "\t"
            sl = False
        str = str + item + " "
    if (str != "\t"):
        out += str
    if (sl):
        return str.rstrip("\\\n") + "\n"
    else:
        return out.rstrip("\\\n") + "\n"


def pretty_format_suf(fname, list, suf, col):
    out = fname + '\n'
    str = "\t"
    for item in list:
        if (len(str) + len(item) + len(suf) + 2 > int(col)):
            out += str + "\\\n"
            str = "\t"
        str = str + item + suf + " "
    if (str != "\t"):
        out += str
    return out.rstrip("\\\n") + "\n"


def get_stdout(cmd):
    nul_f = open(os.devnull, 'w')
    process = sp.Popen(
        [cmd], stdout=sp.PIPE, shell="/bin/bash", stderr=nul_f)
    nul_f.close()
    if sys.version_info >= (3, 0, 0):
        return process.communicate()[0].decode('utf-8')
    else:
        return process.communicate()[0]


def list_info(dir, indent):
    tab = []
    name = " " * (2 * indent) + os.path.basename(dir)
    try:
        file = open(dir + "/info", "r")
        il = 0
        for ll in file:
            tab.append([name if (il == 0) else "", ll.strip()])
            il += 1
        if (il == 0):
            tab.append([name, '\033[93m' + "empty info" + '\033[0m'])
        file.close()
    except IOError:
        if (dir != "problems"):
            tab.append([name, '\033[91m' + "no info" + '\033[0m'])
    for f in sorted(os.listdir(dir)):
        if (os.path.isdir(dir + "/" + f)):
            tab += list_info(dir + "/" + f, indent + 1)
    return tab


def setup_piernik(data=None):
    options, args, all_args, sys_args = piernik_parse_args(data)
    if (options.recycle_cmd):
        if (pickle_avail):
            if (os.path.isfile('.lastsetup')):
                output = open('.lastsetup', 'rb')
                options = pickle.load(output)
                args = pickle.load(output)
                output.close()
            else:
                print("No .lastsetup found")
        else:
            print("Can't read last setup. No Pickle available!")
    else:
        if (pickle_avail):
            output = open('.lastsetup', 'wb')
            pickle.dump(options, output)
            pickle.dump(args, output, -1)
            output.close()

    if (options.verbose):
        print("Setup options:")
        print(options)
        print("Setup arguments:")
        print(args)

    if (options.show_problems):
        tp = list_info("problems", -1)
        maxlen = 0
        for p in tp:
            maxlen = max(maxlen, len(p[0]))
        for p in tp:
            print("%-*s : %s" % (maxlen, p[0], p[1]))

    if (options.show_units):
        print(get_stdout("grep uses ./src/base/units.F90"))

    if (options.show_units or options.show_problems):
        sys.exit()

    if (len(args) < 1):
        sys.stderr.write('\033[91m' + "\nNo problem_name has been provided" + '\033[0m' + "\n")
        exit()

    if (len(args) > 1):
        sys.stderr.write('\033[93m' + "Ignored spurious arguments: " + '\033[0m' + "%s\n" % args[1:])

    # set problem dir
    probdir = 'problems/' + args[0] + '/'
    if (not os.path.isdir(probdir)):
        print("\033[91mCannot find problem directory '%s'." % probdir +
              '\033[0m')
        sys.exit()

    # parse cppflags
    cppflags = ""
    if (options.cppflags):
        for flag_grp in options.cppflags:
            for flag in flag_grp.split(","):
                if (len(flag) > 0):
                    cppflags += ' -D' + flag

    # parse compiler
    if (not re.search(r'\.in$', options.compiler)):
        compiler = options.compiler + '.in'
    else:
        compiler = options.compiler

    objdir = 'obj'
    rundir = 'runs/' + os.path.basename(args[0])
    if (len(options.objdir) > 0):
        objdir += '_' + options.objdir
        rundir += '_' + options.objdir + '/'
    else:
        rundir += '/'

    if (os.path.isdir(objdir)):
        shutil.rmtree(objdir)
    os.mkdir(objdir)

    print("Using compiler settings from \033[93m" + compiler + "\033[0m")
    sc = open(objdir + "/.setup.call", "w")
    sc.write(
        " ".join(sys_args) +
        "\n#effective call (after evaluation of .setuprc*):\n#" + "./setup " +
        " ".join(all_args) + "\n#Using compiler settings from " + compiler + "\n")
    sc.close()

    f90files = []
    allfiles = []
    for f in DirectoryWalker('src'):
        if (is_f90.search(f)):
            f90files.append(f)
        if (is_header.search(f)):
            allfiles.append(f)

    '''Take subproblem files, ignore subdirectories,
    append files from parent directory if their names are new'''

    probfiles = {}
    pdir = os.path.dirname(probdir)  # strip trailing '/'
    while (os.path.basename(pdir) not in ("..", "problems")):
        for f in os.listdir(pdir):
            if (os.path.basename(f) not in probfiles):
                probfiles[os.path.basename(f)] = pdir + '/' + f
        pdir = os.path.dirname(pdir)
    for ff in probfiles:
        f = probfiles[ff]
        if (is_f90.search(f)):
            f90files.append(f)
        if (is_header.search(f) or ff == "piernik.def" or ff == options.param):
            allfiles.append(f)

    piernikdef = ""
    problempar = ""
    pf = probdir + "info"  # intended lack of inheritance from parent directory
    if (not os.path.isfile(pf)):
        print("\033[93mCannot find optional file " + pf + "\033[0m")
    req_prob = ["piernik.def", "initproblem.F90", options.param]
    req_missing = False
    for pf in req_prob:
        found = False
        for lf in allfiles + f90files:
            if (os.path.basename(lf) == pf):
                if (not os.path.isfile(lf)):
                    print("\033[91mRequired file " + lf + " is not a regular file\033[0m")
                    req_missing = True
                else:
                    found = True
                    if (pf == "piernik.def"):
                        piernikdef = lf
                    if (pf == options.param):
                        problempar = lf
        if not found:
            req_missing = True
            print("\033[91mCannot find required file " + pf + "\033[0m")
    if (req_missing):
        sys.exit()

    foo_fd, foo_path = tempfile.mkstemp(suffix=".f90", dir='.')
    cmd = "echo '#include \"%spiernik.h\"' > \"%s\"" % ('src/base/', foo_path)
    cmd += " && cpp %s -dM -I%s \"%s\" && rm \"%s\"" % (
        cppflags, os.path.dirname(piernikdef), foo_path, foo_path)
    defines = get_stdout(cmd).rstrip().split("\n")
    if (options.verbose):
        print(cmd)
        print("Defined symbols:")
        for defin in defines:
            print(defin)

    our_defs = [f.split(" ")[1] for f in filter(cpp_junk.match, defines)]
    our_defs.append("ANY")
    if (options.verbose):
        print("our_defs:")
        print(our_defs)

    files = ['src/base/defines.c']

    uses = [[]]
    incl = ['']
    module = dict()

    for f in f90files:  # exclude links that symbolise file locks
        if (not os.access(f, os.F_OK)):
            print('\033[93m' + "Warning: Cannot access file:" + '\033[0m', f)
            f90files.remove(f)

    for f in f90files:
        for line in open(f):
            if overriding(line):
                o_cnt = 0
                for word in line.split(" "):
                    w = word.rstrip()
                    if f90files.count(w) > 0:
                        o_cnt += 1
                        if (os.path.dirname(f).split("/").count("problems") < 1):
                            print('\033[93m' + "Warning:" + '\033[0m' +
                                  " Only user problems should use the " +
                                  "override feature, not " + f)
                        if os.path.basename(w) == os.path.basename(f):
                            f90files.remove(w)
                            print("Overriding " + w + " by " + f)
                        else:
                            print('\033[93m' + "Warning:" + '\033[0m' +
                                  " Refused overriding " + w + " by " + f +
                                  " due to basename mismatch. Expect errors.")
                if (o_cnt == 0):
                    print('\033[93m' + "Warning:" + '\033[0m' +
                          " No overridable target found for directive '" +
                          line.rstrip() + "'. Expect errors.")

    for f in f90files:
        keys_logic1 = False
        keys_logic2 = False
        keys = []
        luse = []
        linc = []
        for line in open(f):       # Scan original files
            if test(line):
                keys = striplist(line.split(" ")[3:])
            if have_inc(line):
                linc.append(line.split('"')[1])
        # We don't check for recurrsive includes, se we add this one manually
        if ("piernik.h" in linc):
            linc.append("piernik.def")

        keys_logic1 = len(
            keys) == 0 or (len(keys) == 1 and keys[0] in our_defs)
        if (len(keys) == 3):
            keys_logic2 = (
                (keys[1] == "&&" and
                 (keys[0] in our_defs and keys[2] in our_defs)) or
                (keys[1] == "||" and (keys[0] in our_defs or
                                      keys[2] in our_defs)))

        if (keys_logic1 or keys_logic2):
            # workaround the fact that we're using cpp and some of use clauses may depend on __INTEL_COMPILER
            cmd = "cpp %s -I%s -I%s -I%s %s" % (cppflags + " -D__INTEL_COMPILER", probdir, 'src/base', os.path.dirname(piernikdef), f)
            # Scan preprocessed files
            lines = get_stdout(cmd).split('\n')
            for line in lines:
                if have_use(line):
                    luse.append(line.split()[1].rstrip(","))
                if have_mod(line):
                    module.setdefault(line.split()[1],
                                      remove_suf(strip_leading_path([f]))[0])
            files.append(f)
            uses.append(list(set(luse)))
            incl.append(list(set(linc)))

    # for i in iter(module):
    #   if module[i] != i:
    #      print "File",module[i]+".F90 contains an alien module", i

    allfiles.extend(files)

    for f in allfiles:
        if (options.hard_copy):
            shutil.copy(f, objdir)
            # Perhaps we should check for overwriting duplicates here too
        else:
            try:
                os.symlink('../' + f,
                           objdir + '/' + strip_leading_path([f])[0])
            except FileExistsError:
                print("Possible duplicate link or a name clash :", f)
                raise

    if (options.hard_copy):
        otdir = objdir + "/tests/"
        os.mkdir(otdir)
        ctdir = "compilers/tests/"
        for f in os.listdir(ctdir):
            shutil.copy(ctdir + f, otdir)

    if (options.param != 'problem.par'):
        os.symlink(options.param, objdir + '/' + 'problem.par')

    makefile_head = open('compilers/' + compiler, 'r').readlines()
    try:
        makefile_problem = open(probdir + "Makefile.in", 'r').readlines()  # BEWARE: Makefile.in is not inherited from the parent problem yet
    except (IOError):
        makefile_problem = ""
    m = open(objdir + '/Makefile', 'w')

    stripped_files = strip_leading_path(files)
    stripped_files_v = strip_leading_path(files)

    for f in ("version", "os_detect"):  # adding version.F90 and os_detect.F90, which are generated by the Makefile
        stripped_files.append(f + ".F90")
        incl.append('')
        if f == "os_detect":
            uses.append(["constants"])
        else:
            uses.append([])
        module.setdefault(f, f)

    known_external_modules = ["hdf5", "h5lt", "iso_c_binding", "iso_fortran_env", "fgsl",
                              "ifposix", "ifport", "ifcore"]  # Ugly trick: these modules are detected by -D__INTEL_COMPILER
    known_external_modules.append("mpi_f08" if "-DMPIF08" in cppflags.split() else "mpi")

    files_to_build = remove_suf(stripped_files)

    if options.piernik_debug:
        m.write("PIERNIK_DEBUG=1\n")
    for f in makefile_head:
        m.write(f)
    for f in makefile_problem:
        m.write(f)
    m.write(pretty_format_suf("SRCS_V = \\", stripped_files_v, '', columns))
    m.write("SRCS = $(SRCS_V) version.F90\n")
    m.write(pretty_format_suf("OBJS = \\", files_to_build, '.o', columns))
    m.write("\nCPPFLAGS += %s\n" % cppflags)
    if (isinstance(options.f90flags, str)):
        m.write("\nF90FLAGS += %s\n" % options.f90flags)
    if ("SHEAR" in our_defs or "MULTIGRID" in our_defs):
        if ("NO_FFT" not in our_defs):
            m.write("LIBS += $(shell pkg-config --libs fftw3)\n")
        m.write("CPPFLAGS += $(shell pkg-config --libs fftw3 2> /dev/null || echo '-DNO_FFT')\n")
    if ("PIERNIK_OPENCL" in our_defs):
        m.write("LIBS += $(shell pkg-config --libs fortrancl)\n")
        m.write("F90FLAGS += $(shell pkg-config --cflags fortrancl)\n")

    m.write("SILENT = %d\n\n" % (1 if options.laconic else 0))
    m.write(head_block1.replace("./compilers", "") if options.hard_copy else head_block1)
    m.write("\t@( $(ECHO) \"%s\"; \\" % ("./setup " + " ".join(all_args)))
    m.write(head_block2)

    for i in range(0, len(files_to_build)):
        deps = files_to_build[i] + ".o: " + stripped_files[i] + " "
        if (files_to_build[i] == "defines"):
            deps += "piernik.def "
        d = ""
        if (len(incl[i]) > 0):
            d += ' '.join(incl[i]) + ' '
        for j in uses[i]:
            if j in module:
                if module[j] + '.F90' in stripped_files:
                    d += module[j] + '.o '
                else:
                    print("Warning: " + module[j] + '.F90' + "referenced in" +
                          stripped_files[i] + " not found!")
            else:
                if j not in known_external_modules:
                    print("Warning: module " + j + " from file " +
                          stripped_files[i] + " not found!")
        m.write(pretty_format(deps, d.split(), columns))
    m.close()

    # The following code generates simplified module dependency graph data.
    # The simplification is done by removing dependencies of a given element
    # that are common with dependencies of its parent.
    # This is the minimum set of dependencies that are necessary for the make
    # process.
    # This implementation probably may fail on circular dependencies.

    # Collect the data in dictionary of dependence sets
    dep = dict()
    for i in range(0, len(files_to_build)):
        dep[files_to_build[i]] = set()
        for u in uses[i]:
            if u in module:
                if module[u] + '.F90' in stripped_files:
                    dep[files_to_build[i]].add(module[u])

    # save a copy
    dep_s = dict()
    for m in dep:
        dep_s[m] = set(dep[m])

    # construct set of parent dependences
    for m in dep_s:
        dep[m] = set()
        for d in dep_s[m]:
            dep[m] = dep[m].union(dep_s[d])

    # iterate, until all parent dependences are propagated
    # note that cheaper algorithms should exist. We will look for them when
    # we reach few hundreds of modules :-P
    dep_i = dict()
    while True:
        cnt = 0
        for m in dep:
            dep_i[m] = set(dep[m])
        for m in dep:
            for d in dep[m].union(dep_s[m]):
                for pd in dep[m].union(dep_s[d]):
                    if pd not in dep_i[m]:
                        dep_i[m].add(pd)
                        cnt = cnt + 1
        for m in dep_i:
            dep[m] = dep[m].union(dep_i[m])
        if (cnt == 0):
            break

    # assign color to boxes according to their location in source tree
    # some more colors that can be used: tan, lavender, beige, mauve
    colors = {'src': 'red', 'src/base': 'green', 'src/base/mpi': 'seagreen', 'src/grid': 'purple', 'src/particles': 'salmon',
              'src/multigrid': 'orange', 'src/gravity': 'gold3', 'src/scheme': 'brown', 'src/grid/refinement': 'magenta',
              'src/fluids': 'blue', 'src/fluids/cosmicrays': 'cyan', 'src/IO': 'yellow', 'problems': 'pink', 'src/magnetic': 'gray'}

    # create dictionary that maps file to directories
    dirs = {}
    for f in files:
        temp = f.rstrip('.F90').split('/')
        key = temp[-1]
        dirs[key] = [temp[0]]
        if len(temp) > 2:
            for directory in temp[1:-1]:
                dirs[key].append("%s/%s" % (dirs[key][-1], directory))
        dirs[key] = set(dirs[key])

    # write the connectivity file
    dd = open(objdir + '/dep.dot', "w")
    dd.write("digraph piernik {\n")
    dd.write("\t label=\"Dependency graph for the \'" + args[0] + "\' problem\"\n")
    for m in dep:
        for mod in dep_s[m].difference(dep[m]):
            try:
                longest_key = max(dirs[mod].intersection(set(colors.keys())))
                dd.write('\t "%s" [color="%s"];\n' % (mod, colors[longest_key]))
                del dirs[mod]  # prevent duplicates
            except KeyError:
                pass
            dd.write('\t "' + mod + '" -> "' + m + '"\n')
    dd.write('\t subgraph legend {\n\t\t label = "Legend"\n')
    roots = []
    for k in colors:
        dd.write('\t\t "%s" [color="%s"]\n' % (k, colors[k]))
        llen = 0
        ldir = ""
        for kk in colors:
            if k != kk:
                if kk in k:
                    if len(kk) > llen:
                        llen = len(kk)
                        ldir = kk
        if llen > 0:
            dd.write('\t\t "%s" -> "%s"\n' % (ldir, k))
        else:
            roots.append(k)
    dd.write('\t }\n')
    for r in roots:
        dd.write('\t {"piernik" -> "%s" [style = "invisible"; dir=none]}\n' % r)
        # Formally we should calculate the bottom node, not assume it'll be just "piernik"
    dd.write("}\n")
    dd.close()

    fatal_problem = False

    if (not options.nocompile):
        makejobs = ""
        if (mp):
            makejobs = "-j%i" % multiprocessing.cpu_count()
        makecmd = "make %s -C %s" % (makejobs, objdir)
        if (sp.call([makecmd], shell=True) != 0):
            sys.exit('\033[91m' + "It appears that '%s' crashed. Cannot continue." % makecmd + '\033[0m')

    try:
        os.makedirs(rundir)
    except OSError:
        if (not options.keep_par):
            print('\033[93m' + "Found old run." + '\033[0m' +
                  " Making copy of old 'problem.par'.")
            try:
                if (os.path.isfile(rundir + 'problem.par')):
                    shutil.move(rundir + 'problem.par', rundir + 'problem.par.old')
            except (IOError, OSError):
                print('\033[91m' +
                      "Problem with copying 'problem.par' to 'problem.par.old'." +
                      '\033[0m')
        else:
            if (os.path.isfile(rundir + 'problem.par')):
                print('\033[93m' + "Found old run." + '\033[0m' +
                      " Preserving copy of old 'problem.par'.")
        try:
            if (os.path.isfile(rundir + 'piernik')):
                os.remove(rundir + 'piernik')
        except (IOError, OSError):
            print('\033[91m' +
                  "Problem with removing old 'piernik' executable from '%s'." %
                  rundir.rstrip('/') + '\033[0m')
        try:
            if (os.path.isfile(rundir + 'piernik.def')):
                os.remove(rundir + 'piernik.def')
        except (IOError, OSError):
            print('\033[91m' +
                  "Problem with removing old 'piernik.def' from '%s'." %
                  rundir.rstrip('/') + '\033[0m')

    if (options.link_exe):
        try:
            if (os.path.islink(rundir + 'piernik')):
                os.remove(rundir + 'piernik')
        except (IOError, OSError):
            print('\033[91m' +
                  "Problem with removing old 'piernik' executable from '%s'." %
                  rundir.rstrip('/') + '\033[0m')
        try:
            os.symlink("../../" + objdir + "/piernik", rundir + 'piernik')
        except (IOError, OSError):
            print('\033[91m' + "Symlinking 'piernik' failed." + '\033[0m')
            fatal_problem = True
    else:
        if (not options.nocompile):
            try:
                shutil.copy(objdir + "/piernik", rundir + 'piernik')
            except (IOError, OSError, shutil.Error):
                print('\033[91m' + "Copying 'piernik' failed." + '\033[0m')
                fatal_problem = True

    try:
        if (not (options.keep_par and os.path.isfile(rundir + 'problem.par'))):
            shutil.copy(objdir + "/" + options.param, rundir + 'problem.par')
    except IOError:
        print('\033[91m' + "Failed to copy 'problem.par' to '%s'." %
              rundir.rstrip('/') + '\033[0m')
        fatal_problem = True

    # Add a link to original problem directory to make it easier to access auxiliary files stored there
    if options.link_problem:
        orig_prob_link = rundir + "problem"
        if (os.path.islink(orig_prob_link)):
            os.remove(orig_prob_link)
        try:
            os.symlink("../../" + probdir, orig_prob_link)
        except FileExistsError:
            print('\033[91m' + "Cannot create the link to original problem because '" + orig_prob_link + "' exists (and is not a symbolic link)" + '\033[0m')
            fatal_problem = True

    if (options.nocompile):
        print("\033[93mCompilation of '%s' skipped on request." % args[0] +
              """
\033[0mYou may want to run 'make -C %s' before running the Piernik code."""
              % objdir)
        makejobs = ""
        if (mp):
            makejobs = "-j%i" % multiprocessing.cpu_count()
        makecmd = "LC_ALL=C make %s -C %s CHECK_MAGIC=yes piernik" % (
            makejobs, objdir)
        output = sp.Popen(
            makecmd, shell=True, stderr=sp.PIPE, stdout=sp.PIPE).communicate()
        if re.search(r"Circular", output[1].decode("ascii")):
            print('\033[91m' +
                  "Circular dependencies found in '%s'." % objdir + '\033[0m')
    else:
        if (fatal_problem):
            print('\033[93m' +
                  "'%s' compiled, but '%s' may not be ready to run." % (
                      args[0], rundir.rstrip('/')) + '\033[0m')
        else:
            print('\033[92m' + "'%s' ready in '%s'." % (
                args[0], rundir.rstrip('/')) + '\033[0m')


def piernik_parse_args(data=None):
    epilog_help = """Call 'source bin/bash_completion.sh' for some
autocompletion. Frequently used options (like --linkexe, --laconic or
-c <configuration>) can be stored in .setuprc and .setuprc.${HOSTNAME} files.
"""
    usage = "Usage: %prog [problem_name|--last|--problems|--units|--help] [options ...]"
    try:
        parser = OptionParser(usage=usage, epilog=epilog_help)
    except TypeError:
        parser = OptionParser(usage=usage)

    parser.add_option("-v", "--verbose", action="store_true", dest="verbose",
                      default=False, help="""Try to confuse the user with some
diagnostics ;-)""")

    parser.add_option("-q", "--laconic", action="store_true", dest="laconic",
                      default=False, help="Compress the stdout from make.")

    parser.add_option("--debug", action="store_true", dest="piernik_debug",
                      default=False, help="Use debug set of compiler flags.")

    parser.add_option("-n", "--nocompile", action="store_true",
                      dest="nocompile", default=False, help='''Create object
directory, check for circular dependencies, but do not compile.
In combination with --copy will prepare portable object directory.
''')

    parser.add_option("--problems", dest="show_problems", action="store_true",
                      help="Print available problems and exit.")

    parser.add_option("-u", "--units", dest="show_units", action="store_true",
                      help="Print available measurement units and exit.")

    parser.add_option("--last", dest="recycle_cmd", action="store_true",
                      default=False, help="""Call setup with last used options
and arguments (require Pickles).""")

    parser.add_option("--copy", dest="hard_copy", action="store_true",
                      help="Hard-copy source files instead of linking them.")

    parser.add_option("-l", "--linkexe", dest="link_exe", action="store_true",
                      help="""Do not copy obj/piernik to runs/<problem> but link
it instead.""")

    parser.add_option("-p", "--param", dest="param", metavar="FILE",
                      help="Use FILE instead of problem.par .",
                      default="problem.par")

    parser.add_option("-d", "--define", dest="cppflags", metavar="CPPFLAGS",
                      action="append", help="""Add precompiler directives, '-d
DEF1,DEF2' is equivalent to '-d DEF1 -d DEF2' or '--define DEF1 -d DEF2'.""")

    parser.add_option("--f90flags", dest="f90flags", metavar="F90FLAGS",
                      help="Pass additional compiler flags to F90FLAGS .")

    parser.add_option("-c", "--compiler", dest="compiler",
                      default="default.in", help="""Choose specified config from
compilers directory.""", metavar="FILE")

    parser.add_option("-o", "--obj", dest="objdir", metavar="POSTFIX",
                      default='', help="""Use obj_POSTFIX directory instead of obj/ and
runs/<problem>_POSTFIX rather than runs/<problem> .""")

    parser.add_option("-k", "--keeppar", action="store_true", dest="keep_par",
                      help="Do not override existing problem.par file with the default one.")

    parser.add_option("--linkproblem", action="store_true", dest="link_problem",
                      help="Make a symbolic link to the original problem directory in the run directory.")

    if data is None:
        all_args = []
        try:
            for frc in (".setuprc", ".setuprc." + os.uname()[1]):
                for line in open(frc):
                    all_args += line.split()
        except IOError:
            pass
        all_args += sys.argv[1:]
    else:
        all_args = data.split()

    (options, args) = parser.parse_args(all_args)
    if (len(args) < 1 and not (options.recycle_cmd or
                               options.show_problems or
                               options.show_units)):
        parser.print_help()

    return options, args, all_args, sys.argv if data is None else []


if __name__ == "__main__":
    setup_piernik()

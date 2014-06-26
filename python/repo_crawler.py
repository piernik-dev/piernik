#!/usr/bin/env python

import os
import pysvn
import tempfile
import piernik_setup as psetup
import subprocess as sp
import shutil
from drive_via_spawn import parse_mpisignals


PIERNIK_REPO=r"svn+ssh://ladon/piernik/piernik/trunk"
REV_START=7700
REV_END=7733
PROBLEM="mcrwind"
SETUP_CMD="%s --param problem.par.build --compiler gnu47 --debug" % PROBLEM
SETUP_CMD_OLD="%s --param problem.par.build --compiler gnudbg" % PROBLEM
HDF5_FILE="crwind_final_tst_0001.h5"

def get_rev(rev):
    return pysvn.Revision(pysvn.opt_revision_kind.number, rev)


def compare_with_previous(tdir, fname, pref):
    gold_file = os.path.join(tdir, 'previous.h5')
    if os.path.isfile(gold_file):
        if sp.call([os.path.join(test_dir, 'piernik_problem.py'), gold_file, fname]) is not 0:
            shutil.move('diff.png', os.path.join(test_dir, "diff_%s" % pref))
            shutil.move('diff_bare.png',
                        os.path.join(test_dir, "bare_%s" % pref))
    shutil.move(fname, gold_file)
    return

#test_dir = tempfile.mkdtemp()
test_dir = '/dev/shm/dupa'
os.mkdir(test_dir)

client = pysvn.Client()
client.checkout(PIERNIK_REPO, test_dir, revision=get_rev(REV_START))
cwd = os.getcwd()
shutil.copy('piernik_problem.py', test_dir)
os.chdir(test_dir)

prev_rev = REV_START
for rev in range(REV_START, REV_END):
    client.update('.', revision=get_rev(rev))
    try:
        psetup.setup_piernik(SETUP_CMD)
    except IOError:
        psetup.setup_piernik(SETUP_CMD_OLD)

    run_dir = os.path.join(test_dir, "runs", PROBLEM)
    piernik_exe = os.path.join(run_dir, "piernik")
    hdf_file = os.path.join(run_dir, HDF5_FILE)
    if not os.path.isfile(piernik_exe):
        print "Failed to compile piernik, revision will be skipped."
    else:
        sp.call([piernik_exe, "-w", run_dir, "-p", run_dir])

    if not os.path.isfile(hdf_file):
        print "Failed to execute piernik, revision will be skipped."
    else:
        pref = "%i_%i.png" % (prev_rev, rev)
        compare_with_previous(test_dir, hdf_file, pref)
        prev_rev = rev
    shutil.rmtree(run_dir)

print test_dir

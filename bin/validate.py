#!/usr/bin/python

# @todo
#  get rid of paths
#  sanitize
#  problem_name is de facto var passed directly to setup, so it can be
#    anything, make it configurable

import sys
sys.path.append("python/")
import os
import json
import tempfile
import hashlib
import argparse
from DirWalk import DirectoryWalker
try:
    import requests
except ImportError:
    print "You need request: http://pypi.python.org/pypi/requests"
    sys.exit(-1)
try:
    import pysvn
except ImportError:
    print "You need pysvn"
    sys.exit(-1)


def run_jenkins_job(branch, setup_args):
    JENKINS = "http://ladon:8080/job/compile.%s.setup/build" % branch

    file_param = {'name': 'patch.diff', 'file': 'file0'}
    prob_param = {'name': 'problem_name', 'value': setup_args}

    files = {'file0': open(f.name, 'rb')}
    payload = {
        'json': json.dumps({'parameter': [prob_param, file_param]}),
        'Submit': 'Build'
    }
    r = requests.post(JENKINS, data=payload, files=files,
                      auth=requests.auth.HTTPBasicAuth(USER, PASSWORD))


parser = argparse.ArgumentParser(
    description='Validate repo changes using jenkins'
)
parser.add_argument('setup_args', type=str, nargs='?', default='',
                    help="string that is directly passed to setup")
parser.add_argument('-p', '--pretend', action='store_true')
parser.add_argument('-a', '--all', action='store_true')
args = parser.parse_args()

if os.path.exists('.jenkins'):
    f = open('.jenkins', 'rb')
    USER = f.readline().strip()
    PASSWORD = f.readline().strip()
    f.close()
else:
    print("You need to create .jenkins file with:")
    print("user")
    print("pass")
    sys.exit(-1)


svn = pysvn.Client()
f = tempfile.NamedTemporaryFile(delete=False)
diff_data = svn.diff('.', '.')
f.write(diff_data)
f.close()

diff_hash = hashlib.md5()
tab = diff_data.split()
tab.sort()
diff_hash = hashlib.md5()
diff_hash.update("".join(tab))
hash_file = open('.diff_hash', 'w')
hash_file.write(diff_hash.hexdigest())
hash_file.close()

entry = svn.info('.')
url = entry.url.replace(entry.repos + '/', '')
if url.find('branches') == -1:
    branch = url.partition("piernik/")[2].partition("/")[0]
else:
    branch = url.partition("branches/")[2].partition("/")[0]


if args.pretend:
    print("Validation was not performed due to --pretend")
else:
    if args.all:
        for file in DirectoryWalker(os.path.abspath('.')):
            if "initproblem.F90" in file and not os.path.islink(file):
                directory, problem = os.path.split(file)
                if not os.path.exists(os.path.join(directory, 'OBSOLETE')):
                    arg = "%s %s" % (os.path.split(
                        directory)[-1], args.setup_args)
                    run_jenkins_job(branch, arg)
    else:
        if args.setup_args == '':
            args.setup_args = "sedov"
        run_jenkins_job(branch, args.setup_args)
os.unlink(f.name)

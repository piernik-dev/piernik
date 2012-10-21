#!/usr/bin/python

# @todo
#  get rid of paths
#  sanitize
#  problem_name is de facto var passed directly to setup, so it can be
#    anything, make it configurable

import sys
import os
import json
import tempfile
import hashlib
import argparse
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

parser = argparse.ArgumentParser(
   description='Validate repo changes using jenkins'
)
parser.add_argument('setup_args', type=str, nargs=1,
                    help="string that is directly passed to setup")
parser.add_argument('-p', '--pretend', action='store_true')
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
url = entry.url.replace(entry.repos+'/','')
if url.find('branches') == -1:
   branch = url.partition("piernik/")[2].partition("/")[0]
else:
   branch = url.partition("branches/")[2].partition("/")[0]

JENKINS="http://ladon:8080/job/compile.%s.setup/build" % branch

file_param={'name': 'patch.diff', 'file': 'file0'}
prob_param={'name': 'problem_name', 'value': "".join(args.setup_args)}

files={'file0': open(f.name, 'rb')}
payload = {
   'json': json.dumps({'parameter' : [prob_param, file_param]}),
   'Submit': 'Build'
}

if args.pretend:
   print("Validation was not performed due to --pretend")
else:
   r = requests.post(JENKINS, data=payload, allow_redirects=True, files=files,
                     auth=requests.auth.HTTPBasicAuth(USER, PASSWORD))
os.unlink(f.name)

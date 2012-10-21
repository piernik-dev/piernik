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

if os.path.exists('.jenkins'):
   f = open('.jenkins', 'rb')
   USER = f.readline().strip()
   PASSWORD = f.readline().strip()
   f.close()
else:
   print("You need to create .jenkins file with format:")
   print("user")
   print("pass")

f = tempfile.NamedTemporaryFile(delete=False)
svn = pysvn.Client()
f.write(svn.diff('.', '.'))
f.close()
JENKINS="http://ladon:8080/job/compile.trunk.setup/build"

file_param={'name': 'patch.diff', 'file': 'file0'}
prob_param={'name': 'problem_name', 'value': 'kepler'}

files={'file0': open(f.name, 'rb')}
payload = {
   'json': json.dumps({'parameter' : [prob_param, file_param]}),
   'Submit': 'Build'
}

r = requests.post(JENKINS, data=payload, allow_redirects=True, files=files,
                  auth=requests.auth.HTTPBasicAuth(USER, PASSWORD))
os.unlink(f.name)

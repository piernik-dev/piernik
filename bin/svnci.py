#!/usr/bin/python
"""Piernik svnci.py

Usage: svnci.py PATH comment

where comment is "" limited string.
Examples:
   svnci.py src/base/dataio.F90 "changes in dataio.F90"
   svnci.py . "doxygen comment"
"""
import sys
from optparse import OptionParser
import pysvn
import qa

def print_files(list,action):
   print action
   if(not len(list)):
      print "  None"
   else:
      for file in list:
         print "  "+file

def main():
   usage = "usage: %prog [options] PATH COMMENT"
   parser = OptionParser(usage=usage)
   parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="make lots of noise [default]")
   parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose",
                  help="be vewwy quiet (I'm hunting wabbits)")
   parser.add_option("-f", "--force",
                  action="store_true", dest="force",
                  help="commit despite errors (It will be logged)")
   (options, args) = parser.parse_args()
   if len(args) != 2:
      parser.error("incorrect number of arguments")

   client = pysvn.Client()
   entry = client.info('.')
   url = entry.url.replace(entry.repos+'/','')

   if(url.find('branches') == -1):
      comment =  "[source:public/trunk public]: " + sys.argv[2]
   else:
      branch = url.partition("branches/")[2].partition("/")[0]
      comment = "[source:public/branches/"+branch+" "+branch+"]: " + sys.argv[2]

   changes = client.status(sys.argv[1])
   fadd = [f.path for f in changes if f.text_status == pysvn.wc_status_kind.added]
   fmod = [f.path for f in changes if f.text_status == pysvn.wc_status_kind.modified]
   print_files(fadd,      'files to be added')
   print_files(fmod,   'files to be modified')
   print_files([f.path for f in changes if f.text_status == pysvn.wc_status_kind.deleted],    'files to be removed')
   print_files([f.path for f in changes if f.text_status == pysvn.wc_status_kind.conflicted], 'files that are conflicted')
   print_files([f.path for f in changes if f.text_status == pysvn.wc_status_kind.unversioned],'files that are unversioned')
   for item in fmod:
      fadd.append(item)
   qa.qa_checks(fadd,options)
   print "Commiting..."
   client.checkin(sys.argv[1],comment)

if __name__ == "__main__":
   main()

#!/usr/bin/python
"""Piernik svnci.py

Usage: svnci.py PATH comment

where comment is "" limited string.
Examples:
   svnci.py src/base/dataio.F90 "changes in dataio.F90"
   svnci.py . "doxygen comment"
"""
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
   usage = "usage: %prog [-vqpf] PATH -m COMMENT"
   parser = OptionParser(usage=usage)
   parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="make lots of noise [default]")
   parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose",
                  help="be vewwy quiet (I'm hunting wabbits)")
   parser.add_option("-p", "--pretend",
                  action="store_true", dest="pretend",
                  help="I will commit... Haha just kiddin'")
   parser.add_option("-m", "--message",
                  action="store", dest="message", type="string",
                  help="commit message")
   parser.add_option("-f", "--force",
                  action="store_true", dest="force",
                  help="commit despite errors (It will be logged)")
   (options, args) = parser.parse_args()
   if len(args) < 1:
      parser.error("I need at least one file to commit")

   if(not options.message):
      parser.error("Empty message is not allowed")

   client = pysvn.Client()
   entry = client.info('.')
   url = entry.url.replace(entry.repos+'/','')

   if(url.find('branches') == -1):
      rel_path = "public/trunk"
      comment =  "[source:"+rel_path+" public]: " + options.message
   else:
      branch = url.partition("branches/")[2].partition("/")[0]
      rel_path = "public/branches/"+branch
      comment = "[source:"+rel_path+" "+branch+"]: " + options.message

   all_changes = [client.status(f) for f in args]
   changes     = [change for dir_change in all_changes for change in dir_change]
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
   if (options.pretend):
      print "Just kiddin', you used pretend option!"
   else:
      client.checkin(args,comment)

if __name__ == "__main__":
   main()

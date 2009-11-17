#!/usr/bin/python
"""Piernik svnci.py

Usage: svnci.py PATH comment

where comment is "" limited string.
Examples:
	svnci.py src/base/dataio.F90 "changes in dataio.F90"
	svnci.py . "doxygen comments"
"""
import sys
import getopt
import pysvn

def main():
	# parse command line options
	if (len(sys.argv) != 3):
		print __doc__
		sys.exit(2)

	client = pysvn.Client()
	entry = client.info('.')
	url = entry.url.replace(entry.repos+'/','')

	if(url.find('branches') == -1):
		comment =  "[source:public/trunk public]: " + sys.argv[2]
	else:
		branch = url.partition("branches/")[2].partition("/")[0]
		comment = "[source:public/branches/"+branch+" "+branch+"]: " + sys.argv[2]

	client.checkin(sys.argv[1],comment)

if __name__ == "__main__":
	main()

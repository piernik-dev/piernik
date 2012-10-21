'''http://effbot.org/librarybook/os-path-walk-example-3.py'''

import os
import sys


class DirectoryWalker:
    '''Forward iterator that traverses a directory tree'''

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
                try:
                    self.files = os.listdir(self.directory)
                    self.index = 0
                except OSError:
                    print "\033[91mCannot open problem directory '%s'." % \
                        self.directory + '\033[0m'
                    sys.exit(-1)
            else:
                # got a filename
                fullname = os.path.join(self.directory, file)
                if os.path.isdir(fullname) and not os.path.islink(fullname):
                    self.stack.append(fullname)
                return fullname

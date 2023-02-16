#!/usr/bin/env python

import re
import sys

if len(sys.argv) != 2:
    print("Usage: %s LOGFILE" % sys.argv[0])
    sys.exit(1)

f = open(sys.argv[1])
data = f.readlines()
f.close()

namelist = re.compile('\&[A-Z]')
change = re.compile('\*\s[A-Z]')
slash = re.compile('0:\s{1,4}/$')

for line in data:
    line = line.strip()
    if re.search(namelist, line):
        print(' $%s' % (line.split('&')[-1]))
        params = []
        n = 0
    if re.search(change, line):
        ln = line.split('* ')[-1].strip(',')
        ln = re.sub(' *"', '"', ln)
        l1, l2 = ln.split('=')
        l1 = l1.lower().strip()
        l2s = l2.split(',')
        l2n = ''
        for s in l2s:
            l2ns = s.strip()
            if (len(l2ns.split('.')) == 2):
                eend = ''
                l2espl = l2ns.split('E')
                if (len(l2espl) == 2):
                    l2ns, eend = l2espl[0], 'e' + str(int(l2espl[1]))
                while (l2ns[-1] == '0'):
                    l2ns = l2ns[:-1]
                l2ns = l2ns + eend
            if (len(l2n) > 0):
                l2n = l2n + ', '
            l2n += l2ns
        n = max(len(l1), n)
        params.append([l1, l2n])
    if re.search(slash, line.strip()):
        for p in params:
            print('    %s = %s' % (p[0].ljust(n), p[1]))
        print(' /\n')

#!/usr/bin/env python

'''Stupid h5 content comparison script'''

import sys
from yt.mods import load

THRESHOLD = 1e-9

if len(sys.argv) != 3:
    print("Wrong number of arguments!")
    sys.exit(-1)

PF1 = load(sys.argv[1])
PF2 = load(sys.argv[2])

DATA1 = PF1.h.all_data()
DATA2 = PF2.h.all_data()

if not PF1.h.field_list == PF2.h.field_list:
    print("Fields in files differ!")
    sys.exit(-1)

for field in PF1.h.field_list:
    if abs(DATA1[field] - DATA2[field]).max() >= THRESHOLD:
        print("Field %s differs" % field)
        sys.exit(-1)

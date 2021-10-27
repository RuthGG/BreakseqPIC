#!/bin/env python

import sys

from sam import *

handle = open(sys.argv[1], "r") if len(sys.argv)>1 else sys.stdin

for sam in SamParser(handle):
	print "@%s\n%s\n+\n%s" % (sam.qname, sam.seq, sam.qual)

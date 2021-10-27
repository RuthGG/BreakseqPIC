#!/bin/env python

import sys

from sam import *

refs = set()

handle = open(sys.argv[2], "r") if len(sys.argv)>2 else sys.stdin
for sam in  SamParser(handle):
	refs.add(sam.qname)

for sam in SamParser(open(sys.argv[1], "r")):
	if sam.qname in refs:
		print >> sys.stderr, sam
	else:
		print sam

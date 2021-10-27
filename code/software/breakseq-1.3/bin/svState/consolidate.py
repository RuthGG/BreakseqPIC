#!/usr/bin/env python

import sys
from biopy.io import PSL

if len(sys.argv) < 3:
	print "Usage: %s <best hits> <PSL file>..." % sys.argv[0]
	exit()

best=int(sys.argv[1])

for hitList in PSL.hash(PSL.parse(sys.argv[2:], sort=True), max=best).values():
	for hit in hitList: print hit

#!/usr/bin/env python

import sys
from biopy.io import PSL

if len(sys.argv) <= 3:
	print "Usage: %s <Percent ID> <Percent Coverage> <PSL file>..." % sys.argv[0]
	exit()

id=int(sys.argv[1])
cover=int(sys.argv[2])

for file in sys.argv[3:]:
	out = open(file+".i%sc%s.filter"%(id,cover), "w")
	for hit in PSL.parse(file):
		scaledID = hit.percentIdentity *(hit.percentOverlay/100.0)
		if scaledID>=id and hit.percentCoverage>=cover:
			out.write(str(hit)+"\n")
	out.close()	

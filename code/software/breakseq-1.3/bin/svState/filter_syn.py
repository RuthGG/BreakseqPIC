#!/usr/bin/env python

import sys
from biopy.io import SV
from biopy.io import PSL
from biopy.io import Net
from biopy.io.CList import Interval

if len(sys.argv) <= 3:
	print "Usage: %s <SV file> <Net file> <PSL file>..." % sys.argv[0]
	exit()

svs=SV.hash(SV.parse(sys.argv[1]))
nets=Net.hash(Net.parse(sys.argv[2]))

def is_syntenic(sv, hit):
	syntenic=False
	if sv.name in nets:
		fills=nets[sv.name]
		for fill in fills:
			target=Interval(fill.tName, fill.tStart, fill.tEnd)
			if sv.is_overlapping(target):
				query=Interval(fill.qName, fill.qStart, fill.qEnd)
				syntenic=hit.is_overlapping(query)
				if syntenic: break
	return syntenic

for file in sys.argv[3:]:
	out = open(file+".syn", "w")
	for hit in PSL.parse(file):
		sv=svs[hit.qName]
		sv=Interval(sv.name, sv.start, sv.end)
		if is_syntenic(sv, Interval(hit.tName, hit.tStart, hit.tEnd)):
			out.write(str(hit)+"\n")
	out.close()

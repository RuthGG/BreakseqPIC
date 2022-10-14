#!/usr/bin/env python

import sys

from biopy.io import SV

if len(sys.argv) < 2:
	print "Usage: %s <txt file>" % sys.argv[0]
	exit()

def set_attribute(sv, key, value):
	if value is not None and value != "" and value != "N.A":
		sv.attributes[key]=value

txt=open(sys.argv[1])
txt.readline()

for line in txt:
	cols=line.rstrip().split("\t")
	sv=SV.Call()
	sv.name=cols[4]
	sv.source=cols[0]
	sv.feature=cols[3]
	sv.start=int(cols[12])
	sv.end=int(cols[13])
	set_attribute(sv, "Phase", cols[1])
	set_attribute(sv, "Id", cols[2])
	set_attribute(sv, "Iseq", cols[14])
	set_attribute(sv, "NonTempltIseq", cols[11])
	set_attribute(sv, "Microhomology", cols[16])
	print sv

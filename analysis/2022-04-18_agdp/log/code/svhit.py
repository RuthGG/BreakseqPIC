#!/bin/env python

import sys, os

from sam import *

files=sys.argv[1:]

calls={}
cats={"A":0,"B":1,"C":2}

for file in files:
	handle = open(file, "r")
	for sam in SamParser(handle):
		key = sam.rname 
		cat = 0
		if key[-2]==":":
			cat=cats[key[-1]]
			key=key[0:-2]
		if key not in calls:
			calls[key]=[0,0,0]
		calls[key][cat]+=1
	handle.close()

for call in sorted(calls.items()):
	print "%s\t%s\t%s\t%s\t%s" % (call[0], call[1][0], call[1][1], call[1][2], sum(call[1]))

#!/usr/bin/env python

import sys

from biopy.io import PSL
from biopy.io import Fasta
from biopy.io import SV

if len(sys.argv) < 6:
	print "Usage: %s <SV file> <PSL A> <PSL B> <PSL C> <Genome Dir>" % sys.argv[0]
	exit()

setA=PSL.hash(PSL.parse(sys.argv[2], sort=True), max=1)
setB=PSL.hash(PSL.parse(sys.argv[3], sort=True), max=1)
setC=PSL.hash(PSL.parse(sys.argv[4], sort=True), max=1)

seqs=Fasta.Seqs(sys.argv[5])

def get_comment(setX, label, id):
	return "# %s\t"%label + ("N/A" if id not in setX else "\t".join(setX[id][0].summary().split("\t")[3:]))

def get_hit(id, func):
	if id in setA:
		if id not in setB: return setA[id]
		else: return func(setA[id], setB[id])
	elif id in setB:
		return setB[id]
	else: 
		return None

for sv in SV.parse(sys.argv[1], base=seqs):
	rectified=0
	rectifiable=sv.is_indel() and ((sv.id in setA and sv.id in setB) or (sv.id in setC))
	if rectifiable:
		rectified=1
		if sv.is_insertion():
			worst=get_hit(sv.id, min)
			if sv.id not in setC or (sv.id in setA and sv.id in setB and setC[sv.id]<worst):
				rectified=2
		elif sv.is_deletion():
			best=get_hit(sv.id, max)
			if sv.id in setC and (best is None or setC[sv.id]>best):
				rectified=2
	if rectified==2:
		sv.event(flip=True)
	sv.rect(rectified)
	print sv

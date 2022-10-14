#!/usr/bin/env python

import sys
from biopy.io import SV

if len(sys.argv) <= 3:
	print "Usage: %s <Original SV GFF> <Rectified SV GFF 1> <Rectified SV GFF 2>..." % sys.argv[0]
	exit()

sv_calls=SV.parse(sys.argv[1])
gffs=[]
for file in sys.argv[2:]:
	gffs.append(SV.parse(file))

for i in range(len(sv_calls)):
	sv=sv_calls[i]
	statuses=[]
	for gff in gffs:
		statuses.append(gff[i].rect())
	for s in statuses:
		status=int(s)
		if status>0:
			event=sv.event()
			if status==2: event=sv.event_flipped()
			sv.ancestral_state(event)
			break
	sv.rect(":".join(statuses))
	print sv

#!/usr/bin/env python

import sys

from biopy.io import SV

if len(sys.argv) < 2:
	print "Usage: %s <SV file>..." % sys.argv[0]
	exit()

class counter:

	def __init__(self):
		self.counts={}

	def increment(self, key=None, value=1):
		if not isinstance(key, tuple):
			key=(key,"Count")
		if key[0] not in self.counts:
			self.counts[key[0]]={}
		if key[1] not in self.counts[key[0]]:
			self.counts[key[0]][key[1]]=value
		else:
			self.counts[key[0]][key[1]]+=value

	def __str__(self):
		s="Type"
		total=0
		cols=set()
		for row in self.counts:
			for col in self.counts[row]:
				cols.add(col)
		
		for col in sorted(cols):
			s+="\t"+col

		for rowindex in sorted(self.counts.keys()):
			s += "\n"+rowindex
			for col in sorted(cols):
				row=self.counts[rowindex]
				count = 0 if col not in row else row[col]
				total+=count
				s += "\t%s" % count
		s += "\nTOTAL\t%s\n" % total
		return s

def count_rectifiable(rectification, sv, element):
	s = sv.is_rectified()
	r = sv.is_rectifiable()
	if r:
		rectification.increment(("Rectifiable%s"%r,element))
		rectification.increment(("Rectifiable",element))
		c = sv.is_consistently_rectifiable()
		if c:
			rectification.increment(("Consistently Rectifiable%s"%c,element))
			rectification.increment(("Consistently Rectifiable",element))
	else:
		rectification.increment(("Unrectifiable",element))
	rectification.increment(("Rectified%s"%s,element))


bounds=[i for i in range(10000, 100000, 10000)]
bounds.append(sys.maxint)

def count_bounds(boundaries, sv, element):
	bp = int(sv.trace().split(":")[-1])
	for bound in bounds:
		if bp<=bound:
			boundaries.increment((">" if bound==sys.maxint else str(bound), element))
			break

def bcounter():
	c=counter()
	for bound in bounds: c.increment((">" if bound==sys.maxint else str(bound), "-"), 0)
	return c

for f in sys.argv[1:]:
	sources=counter()
	mechs=counter()
	events=counter()
	both=counter()
	rectification=counter()
	rectification_breakdown=counter()
	rearrange=counter()
	rectifiables=counter()
	rect_rearrange=counter()
	prox_count=bcounter()
	prox_rect_count=bcounter()
	prox_conr_count=bcounter()
	for sv in SV.parse(f):
		sources.increment(sv.source)
		mechs.increment(sv.mech())
		events.increment(sv.event())
		both.increment((sv.event(),sv.mech()))
		count_rectifiable(rectification, sv, "N.A." if "State" not in sv.attributes else ("R."+sv.attributes["State"]))
		count_rectifiable(rectification_breakdown, sv, sv.mech())
		if sv.is_rectifiable():
			rectifiables.increment(("R."+sv.attributes["State"],sv.mech()))
		if sv.is_consistently_rectifiable():
			if sv.is_intra(): rect_rearrange.increment(("CR.Intra",sv.mech()))
			elif sv.is_inter(): rect_rearrange.increment(("CR.Inter",sv.mech()))
		if sv.is_intra():
			count_bounds(prox_count, sv, sv.mech())
			if sv.is_rectifiable(): count_bounds(prox_rect_count, sv, sv.mech())
			if sv.is_consistently_rectifiable(): count_bounds(prox_conr_count, sv, sv.mech())		
		if sv.is_intra(): rearrange.increment(("Intra",sv.mech()))
		elif sv.is_inter(): rearrange.increment(("Inter",sv.mech()))

	print f
	print sources
	print events
	print mechs
	print both
	print rectification
	print rectification_breakdown
	print rearrange
	print rectifiables
	print rect_rearrange
	print prox_count
	print prox_rect_count
	print prox_conr_count

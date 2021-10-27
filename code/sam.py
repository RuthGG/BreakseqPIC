#!/bin/env python

import sys

class Sam:
	
	def __init__(self, aline):
		fields=aline.split("\t")
		self.aline=aline
		self.qname=fields[0]
		self.flag=int(fields[1])
		self.rname=fields[2]
		self.pos=int(fields[3])
		self.mapq=int(fields[4])
		self.ciagr=fields[5]
		self.mrnm=fields[6]
		self.mpos=int(fields[7])
		self.isize=int(fields[8])
		self.seq=fields[9]
		self.qual = fields[10]
		self.opt = fields[11] if len(fields)>11 else None

	def __str__(self):
		return self.aline

class SamParser:

	def __init__(self, handle):
		self.handle=handle

	def __iter__(self):
		return self

	def next(self):
		line=self.handle.readline()
		while line is not None and line.startswith("@"):
			line=self.handle.readline()
		line=line.strip()
		if line is None or len(line)==0:
			raise StopIteration 
		else:
			try:
				return Sam(line)
			except Exception as e:
				print >> sys.stderr, "ALINE: "+line
				raise e

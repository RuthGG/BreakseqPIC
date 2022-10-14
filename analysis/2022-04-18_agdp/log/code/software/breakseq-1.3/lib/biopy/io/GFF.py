#!/usr/bin/env python

import sys

from biopy.io.CList import Interval

class Entry:

	def __init__(self, entry=None):
		fields = entry.rstrip().split("\t") if entry is not None else ["","","","0","0",".",".","."]
		has_attr = len(fields)>8
		for i in 3,4: fields[i]=int(fields[i])
		self.name, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame=fields[0:-1] if has_attr else fields[0:]
		self.attributes={}
		if len(fields)>8:
			for attr in fields[-1].split(";"):
				attr=attr.strip()
				if len(attr)==0: continue
				sp=attr.find(" ")
				attr=(attr[0:sp],attr[sp+1:].strip("\""))
				self.attributes[attr[0]]=attr[-1]

	def __str__(self):
		return "\t".join([self.name, self.source, self.feature, str(self.start), str(self.end), self.score, self.strand, self.frame])+"\t"+"; ".join(["%s \"%s\"" % i for i in sorted(self.attributes.items())])

	def interval(self):
		return Interval(self.name, self.start, self.end)

def parse(file):
	entries=[]
	for line in open(file, "r"):
		if line.startswith("#"): continue
		try:
			entry=Entry(line)
			entries.append(entry)
		except:
			print >> sys.stderr, "Unable to parse line: %s" % line
			raise
	return entries 

def hash(entries):
	hash={}
	for entry in entries:
		if entry.feature not in hash:
			hash[entry.feature]=[]
		hash[entry.feature].append(entry)
	return hash

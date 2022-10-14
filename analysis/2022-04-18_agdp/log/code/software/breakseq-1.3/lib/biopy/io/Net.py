#!/usr/bin/env python

import sys

class Entry:

	def __init__(self, entry):
		self.fields=entry.strip().split(" ")
		self.type=self.fields[0]
		self.level=self.__level(entry)
		if self.type=="net":
			self.tName=self.fields[1]
			self.tSize=int(self.fields[2])
		else:
			self.tName=""
			for i in (1,2,5,6): self.fields[i]=int(self.fields[i])
			self.tStart,self.tSize,self.qName,self.qStrand,self.qStart,self.qSize=self.fields[1:7]
			self.tEnd=self.tStart+self.tSize
			self.qEnd=self.qStart+self.qSize
		self.children=[]
		
	def __str__(self):
		return ("%s\t%s"%(self.tName,self.tSize)) if self.level==0 else "\t".join([str(s) for s in (self.fields[0:7])])

	def __cmp__(self, other):
		return cmp(self.level, other.level)

	def __level(self, entry):
		i=0
		for c in entry: 
			if c==" ": i+=1 
			else: break
		adj=i%2==1
		if adj: i=i-1
		return i/2 + (1 if adj else 0)

	def isFill(self):
		return self.type=="fill"

	def isGap(self):
		return self.type=="gap"

	def add(self, child):
		child.tName=self.tName
		self.children.append(child)

	def accept(self, child):
		if self.level==0:
			return child.isFill() and child.level==1
		else:
			return self.type!=child.type and (self.level == child.level - (1 if child.isFill() else 0))

	def dump(self, max=-1, entry=None, indent=0):
		if entry is None:
			entry=self
		if max<0 or entry.level<=max:
			tabs=""
			for i in range(indent):
				tabs+="\t"
			print tabs+entry.type+"\t"+str(entry.coords())
			for child in entry.children:
				self.dump(max, child, indent+1)

	def coords(self):
		return "%s:%s-%s" % (self.tName, 1, self.tSize) if self.level==0 else "%s:%s-%s\t%s:%s-%s:%s" % \
			(self.tName, self.tStart+1, self.tEnd, self.qName, self.qStart+1, self.qEnd, self.qStrand)

def __build(entries, pos=0, parent=None):
	while pos<len(entries):
		curr=entries[pos]
		if parent is None or curr.level==0:
			parent=curr
			pos+=1
			continue
		if parent.accept(curr):
			parent.add(curr)
			entries[pos]=None
			pos+=1
			if pos<len(entries):
				next = entries[pos]
				if not parent.accept(next):
					if parent <= next:
						pos=__build(entries, pos, curr)
					else:
						break
		else:
			break
	return pos

def parse(file, max=-1):
	entries=[]
	for line in open(file, "r"):
		entry=Entry(line)
		if max<0 or entry.level<=max:
			entries.append(entry)
	__build(entries)
	nets=[]
	for entry in entries:
		if entry is not None:
			nets.append(entry)
	return nets

def hash(nets):
	hash={}
	for net in nets:
		hash[net.tName]=net.children
	return hash

#for i in range(-1,-2,-1):
#	print >> sys.stderr, "LEVEL %s"%i
#	for entry in parse(sys.argv[1], i):
#		entry.dump(i)

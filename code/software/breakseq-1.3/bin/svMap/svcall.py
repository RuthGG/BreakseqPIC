#!/bin/env python

import sys, os
import math

from sam import *

uni=open(sys.argv[1], "r") if len(sys.argv)>1 else sys.stdin
xun=open(sys.argv[2], "r") if len(sys.argv)>2 else None

class Call:

	def __init__(self, cline):
		fields=cline.rstrip().split("\t")
		cols=len(fields)
		self.id=fields[0]
		self.hits_A =int(fields[1]) if cols>1 else 0
		self.hits_B =int(fields[2]) if cols>2 else 0
		self.hits_C =int(fields[3]) if cols>3 else 0
		self.hits=int(fields[4]) if cols>4 else 0
		self.hits_log=0 if self.hits==0 else math.log(self.hits, 2)
		self.score = -1
		self.refs_log=-1
		self.ref_hits=-1

	def add_ref(self, ref_call):
		self.ref_call = ref_call
		self.ref_hits = 0 if ref_call is None else self.ref_call.hits
		self.refs_log = 0 if ref_call is None else self.ref_call.hits_log
		self.score = max(0, self.hits_log-self.refs_log)
	
	def __cmp__(self, other):
		if self.score < other.score:
			return -1
		elif self.score > other.score:
			return 1
		elif self.ref_hits < other.ref_hits:
			return 1
		elif self.ref_hits > other.ref_hits:
			return -1
		elif self.hits < other.hits:
			return -1
		elif self.hits > other.hits:
			return 1
		else:
			return 0

	def __str__(self):
		return "%s\t%s\t%s\t%s\t%s\t%s\t%.2f" % (self.id, self.hits_A, self.hits_B, self.hits_C, self.hits, (0 if self.ref_call is None else self.ref_call.hits), self.score)

def read_calls(handle):
	calls={}
	if handle is not None:
		for cline in handle:
			call = Call(cline)
			calls[call.id]=call
	return calls

uni_calls=read_calls(uni)
xun_calls=read_calls(xun)

for x in xun_calls:
	if x not in uni_calls:
		c = Call(x)
		uni_calls[x]=c

counts=[0,0,0]
for call in uni_calls.values():
	ref_call = None	if call.id not in xun_calls else xun_calls[call.id]
	call.add_ref(ref_call)
	if call.score>2: counts[0]+=1
	elif call.score>0: counts[1]+=1
	else: counts[2]+=1

for call in sorted(uni_calls.values(), reverse=True):
	print call

print >> sys.stderr, "Total: %s; High-support: %s; Low-support: %s; No-support: %s" % (sum(counts), counts[0], counts[1], counts[2])

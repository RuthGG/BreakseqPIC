#!/usr/bin/env python

import sys
from math import log

class Hit:

	dataLabels = [	'matches','misMatches','repMatches','nCount',
			'qNumInsert','qBaseInsert','tNumInsert','tBaseInsert',
			'strand','qName','qSize','qStart','qEnd','tName','tSize','tStart','tEnd',
			'blockCount','blockSizes','qStarts','tStarts' ]

	stringFields = [ 8, 9, 13, 18, 19, 20 ]

	def __init__(self, entry, isProtein=False):
		self.entry=entry.rstrip()
		fields=self.entry.split("\t")
		for lIndex in xrange(len(self.dataLabels)): 
			exec ('self.%s = "%s"' if lIndex in self.stringFields else 'self.%s = int("%s")') % (self.dataLabels[lIndex], fields[lIndex])
		self.isProtein = isProtein
		self.percentIdentity = self.calcPercentIdentity()
		self.score = self.calcScore()
		self.qHitSize = self.qEnd-self.qStart
		self.tHitSize = self.tEnd-self.tStart
		self.percentCoverage = 100 * float(self.qHitSize)/self.qSize
		self.percentOverlay  = 100 * float(min(self.qHitSize, self.tHitSize))/max(self.qHitSize, self.tHitSize)

	def summary(self):
		return "\t".join([str(s) for s in (self.qName, self.qStart, self.qEnd, self.tName, self.tStart, self.tEnd, self.strand, self.percentIdentity, self.percentCoverage, self.score)])

	def __str__(self):
		hit="\t".join([eval("str(self.%s)"%l) for l in self.dataLabels])
		return "\t".join([str(s) for s in (hit, self.percentIdentity, self.percentCoverage, self.percentOverlay, self.score)])

	def __cmp__(self, other):
		c = cmp(self.score, other.score)
		if c==0: c = cmp(self.percentOverlay, other.percentOverlay)
		return c

	def __calcMilliBad(self, isMrna):
		sizeMul = 3 if self.isProtein else 1
		qAliSize = sizeMul * (self.qEnd - self.qStart)
		tAliSize = self.tEnd - self.tStart
		aliSize = min(qAliSize, tAliSize)
		milliBad = 0
		if (aliSize <= 0):
			return 0
		sizeDif = qAliSize - tAliSize
		if (sizeDif < 0):
			if (isMrna):
				sizeDif = 0
			else:
        			sizeDif = -sizeDif
		insertFactor = self.qNumInsert
		if (not isMrna):
			insertFactor += self.tNumInsert
		total = (sizeMul * (self.matches + self.repMatches + self.misMatches))
		if (total != 0):
			milliBad = (1000 * (self.misMatches*sizeMul + insertFactor + round(3*log(1+sizeDif)))) / total
		return milliBad

	def calcPercentIdentity(self, isMrna=True):
		return 100.0 - self.__calcMilliBad(isMrna) * 0.1

	def calcScore(self):
		sizeMul = 3 if self.isProtein else 1
		return sizeMul * (self.matches + ( self.repMatches>>1)) - sizeMul * self.misMatches - self.qNumInsert - self.tNumInsert


def parse(files, sort=False):
	if not isinstance(files, list):
		files=[files]
	hits=[]
	for file in files:
		handle=open(file, "r")
		for entry in handle:
			if entry.startswith("psLayout"):
				for header in handle:
					if header.startswith("-"):
						break
				continue
			try:
				hits.append(Hit(entry))
			except ValueError:
				print >> sys.stderr, entry,
				continue
	if sort: hits.sort(reverse=True)
	return hits

def hash(hits, max=0):
	hash={}
	for hit in hits:
		if hit.qName not in hash:
			hash[hit.qName]=[]
		if max<=0 or len(hash[hit.qName]) < max: hash[hit.qName].append(hit)
	return hash

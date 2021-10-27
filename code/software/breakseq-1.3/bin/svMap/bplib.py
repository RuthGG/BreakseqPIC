#!/bin/env python
import sys
from Bio import SeqIO
from biopy.io import Fasta
from biopy.io import SV

if len(sys.argv)<=3:
	print "Usage: %s <SV GFF file> <reference dir> <window>" % sys.argv[0]
	exit()

reference=Fasta.Seqs(sys.argv[2], int(sys.argv[3]))

def print_seq(id, type, seq):
	print ">%s:%s\n%s" % (id, type, seq)

for sv in SV.parse(sys.argv[1], reference):
	flanks=sv.get_flanks()
	if sv.is_insertion():
		print_seq(sv.id, "A", flanks[0])
		print_seq(sv.id, "B", flanks[1])
	elif sv.is_deletion():
		print_seq(sv.id, "C", flanks[2])

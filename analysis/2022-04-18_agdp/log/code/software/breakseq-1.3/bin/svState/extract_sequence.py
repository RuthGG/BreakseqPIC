#!/bin/env python
import sys, os
from Bio import SeqIO
from biopy.io import Fasta
from biopy.io import SV

if len(sys.argv)<4:
	print "Usage: %s <SV GFF file> <genome dir> <output dir> [window=1000]" % sys.argv[0]
	exit()

input_file=sys.argv[1]
genome_dir=sys.argv[2]

out_prefix=sys.argv[3]+"/"+os.path.basename(input_file).replace(".gff", "")

outs=[open(out_prefix+".A.fa", "w"), open(out_prefix+".B.fa", "w"), open(out_prefix+".C.fa", "w")]
seqs=Fasta.Seqs(genome_dir, window=1000 if len(sys.argv)<5 else int(sys.argv[4]))

for sv in SV.parse(input_file, seqs):
	flanks=sv.get_flanks()
	for i in range(0,len(outs)):
		outs[i].write(">%s\n" % sv.id)
		outs[i].write(str(flanks[i])+"\n")

for o in outs: o.close()

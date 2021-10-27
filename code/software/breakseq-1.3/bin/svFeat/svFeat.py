#!/usr/bin/env python

import sys, os

from biopy.io import Config

if len(sys.argv) <= 2:
	print "%s <SV input GFF> <output dir>" % sys.argv[0]
	exit()

outdir=sys.argv[2]
outpre=outdir+"/"+os.path.basename(sys.argv[1]).replace(".gff", "")

config=Config.load()

seqdir=config["Genome_Seq_Directory"]
gencode=config["Gencode_GFF"]

print "Calculating physical features and intersecting gene annotation"

os.system("%s/build_features.py %s %s %s %s 50" % (sys.path[0], sys.argv[1], seqdir, gencode, outdir))

print "Finished feature analysis"

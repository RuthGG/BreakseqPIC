#!/usr/bin/env python

import sys, os

from biopy.io import Config

if len(sys.argv) <= 2:
	print "%s <SV input GFF> <output dir> [window=1000]" % sys.argv[0]
	exit()

outdir=sys.argv[2]
outpre=outdir+"/"+os.path.basename(sys.argv[1]).replace(".gff", "")
window=1000 if len(sys.argv) <=3 else int(sys.argv[3])

config=Config.load()

seqdir=config["Genome_Seq_Directory"]
server=config["Blat_Server"]
client=config["Blat_Client"]
twobit=config["Blat_TwoBit"]
synnet=config["Net_Chain"]

print "> Extracting sequences for alignment"

os.system("%s/extract_sequence.py %s %s %s %i" % (sys.path[0], sys.argv[1], seqdir, outdir, window))

for g in server.keys():
	print "> Aligning to ancestral genome %s and performing state analysis"%g
	g_outpre="%s.%s" % (outpre, g)
	for t in ("A", "B", "C"):
		input="%s.%s.fa" % (outpre,t)
		output="%s.%s" % (g_outpre, t)
		os.system("%s %s 8080 -minScore=0 -minIdentity=90 %s %s %s.psl" % (client, server[g], twobit[g], input, output))
		os.system("%s/filter_idt.py 90 90 %s.psl" % (sys.path[0], output))
		os.system("%s/filter_syn.py %s %s %s.psl.i90c90.filter" % (sys.path[0], sys.argv[1], synnet[g], output))
		os.system("%s/consolidate.py 1 %s.psl.i90c90.filter.syn > %s.best.psl" % (sys.path[0], output, output))
	os.system("%s/rectify.py %s %s.*.best.psl %s > %s.rectified.gff" % (sys.path[0], sys.argv[1], g_outpre, seqdir, g_outpre))

print "> Integrating results from multiple ancestral analyses"

os.system("%s/merge.py %s %s.*.rectified.gff > %s.gff" % (sys.path[0], sys.argv[1], outpre, outpre))

print "> Finished ancestral state analysis"

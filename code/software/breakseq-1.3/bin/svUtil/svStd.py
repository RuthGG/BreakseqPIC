#!/usr/bin/env python

import sys, os

from biopy.io import SV

if len(sys.argv) <= 2:
	print "%s <SV input GFF> <output dir>" % sys.argv[0]
	exit()

outdir=sys.argv[2]
outpre=outdir+"/"+os.path.basename(sys.argv[1]).replace(".gff", "")

out_gff=open(outpre+".gff", "w")
out_fna=open(outpre+".fna", "w")
for sv in SV.parse(sys.argv[1]):
	# filtering can be done here
	out_gff.write(str(sv)+"\n")
	if sv.is_insertion():
		out_fna.write(">%s\n%s\n"%(sv.id,sv.get_sequence()))
out_gff.close()
out_fna.close()

print "> Finished library standardization"

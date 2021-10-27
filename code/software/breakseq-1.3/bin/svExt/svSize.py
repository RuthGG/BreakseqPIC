#!/usr/bin/env python

import sys, os

from biopy.io import SV

if len(sys.argv) <= 3:
	print "%s <SV input GFF> <output dir> <min size> [max size (exclusive)]" % sys.argv[0]
	exit()

minsize=int(sys.argv[3])
maxsize=sys.maxint if len(sys.argv)<5 else int(sys.argv[4])
outdir=sys.argv[2]
outpre=outdir+"/"+os.path.basename(sys.argv[1]).replace(".gff", ".%i-%sbp"%(minsize,maxsize if maxsize<sys.maxint else "max"))

print "> Selecting SVs by size: >=%i and <%s" % (minsize, maxsize if maxsize<sys.maxint else "max")

out_gff=open(outpre+".gff", "w")
out_fna=open(outpre+".fna", "w")
for sv in SV.parse(sys.argv[1]):
	size=sv.size()
	if size>=minsize and size<maxsize:
		out_gff.write(str(sv)+"\n")
		if sv.is_insertion():
			out_fna.write(">%s\n%s\n"%(sv.id,sv.get_sequence()))
out_gff.close()
out_fna.close()

print "> Finished SV selection"

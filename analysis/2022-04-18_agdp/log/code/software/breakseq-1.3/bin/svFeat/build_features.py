#!/usr/bin/env python

import sys, os
import itertools

from biopy.io import Fasta
from biopy.io import SV
import Feature
import Fragility
import GC
import Gene

if len(sys.argv) < 5:
	print "Usage %s <SV GFF file> <Sequence Dir> <gencode GFF file> <output dir> [window=50]" % sys.argv[0]
	exit()

window=50 if len(sys.argv)<6 else int(sys.argv[5])
gencode=Gene.load_genes(sys.argv[3])
outdir=sys.argv[4]
outpre=outdir+"/"+os.path.basename(sys.argv[1]).replace(".gff","")

sv_calls = SV.parse(sys.argv[1], base=Fasta.Seqs(sys.argv[2], window=window))


seq_feats={\
        "Flex":Fragility.calc_flexibility,
        "Helix":Fragility.calc_stability,
        "GC":GC.calc_content
}
freq_feats={}
other_feats={\
	"Gene":(Gene.find_overlap, gencode)
}

def write_outs(outs, results=None, trailing=None, leading=None):
	if isinstance(results, str):
		results=(results,)
	for i in range(len(outs)):
		if leading  is not None: outs[i].write(leading)
		if results  is not None: outs[i].write(str(results[0] if len(results)==1 else results[i]))
		if trailing is not None: outs[i].write(trailing)

output=open(outpre+".gff", "w")
outfile=outpre+".feats.win%s.tab"%(window)
outs=[open(outfile, "w"), open(outfile+".A", "w"), open(outfile+".B", "w")]
header="ID\t"+"\t".join([feature for feature in itertools.chain(sorted(seq_feats.keys()), sorted(freq_feats.keys()), sorted(other_feats.keys()))])
write_outs(outs, header, "\n")
for sv in sv_calls:
	flanks=sv.get_flanks()
	write_outs(outs, sv.id)
	for feat in sorted(seq_feats.keys()):
		results=Feature.calc_property(seq_feats[feat], \
			(flanks[0],), 
			(flanks[1],))
		write_outs(outs, results, leading="\t")
		sv.attributes[feat]=results[0]
	for feat in sorted(freq_feats.keys()):
		results=Feature.calc_property(Feature.calc_freq, \
			(sv.name, sv.start, sv.start, freq_feats[feat]), 
			(sv.name, sv.end, sv.end, freq_feats[feat]))
		write_outs(outs, results, leading="\t")
		sv.attributes[feat]=results[0]
	for feat in sorted(other_feats.keys()):
		results=other_feats[feat][0](\
			sv.name, sv.start, sv.end, other_feats[feat][1])
		write_outs(outs, results, leading="\t")
		sv.attributes[feat]=results
	write_outs(outs, "\n")
	output.write(str(sv)+"\n")

output.close()
for out in outs: out.close()

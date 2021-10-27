from biopy.io.CList import *
from biopy.io import GFF

def find_overlap(name, start, end, intervals):
	result=[]
        i = Interval(name, start, end)
        found=intervals.find(i)
	for f in found:
		a=":".join(f.annotations)
		result.append(a)
	result.sort()
	return ",".join(result)

def load_genes(gencode_gff):
	genes=Intervals()
	for entry in GFF.parse(gencode_gff):
		if entry.feature == "gene":
			if "gene_type" in entry.attributes and entry.attributes["gene_type"]=="pseudogene":
				continue
			i = Interval(entry.name, entry.start, entry.end, [entry.attributes["gene_id"]])
			genes.add(i)
	genes.sort()
	return genes

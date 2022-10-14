#!/bin/env python

import sys, os

from biopy.io import SV

if len(sys.argv) <= 2:
	print "Usage: %s <SV GFF file> <output dir> [minimum size (bp) of events: 200 (default), 50 (min)]" % sys.argv[0]
	exit()

outdir=sys.argv[2]
prefix=os.path.basename(sys.argv[1]).replace(".gff","")
outpre=outdir+"/"+prefix

print "> Preprosssing for mechanism classification"

del_file=outpre+".del.txt"
ins_file=outpre+".ins.txt"
inv_file=outpre+".inv.txt"

del_out=open(del_file, "w")
ins_out=open(ins_file, "w")
inv_out=open(inv_file, "w")

del_out.write("chr\tstart\tend\tid\n")
ins_out.write("chr\tstart\tend\tinserted\tid\n")
inv_out.write("chr\tstart\tend\tid\n")

svs = SV.parse(sys.argv[1])

for sv in svs:
	chr=sv.name.replace("chr","")
	if sv.is_deletion():
		del_out.write("\t".join([chr, str(sv.start-1), str(sv.end), sv.id])+"\n")
	elif sv.is_insertion():
		ins_out.write("\t".join([chr, str(sv.start), str(sv.end-1)])+"\t"+sv.get_sequence()+"\t%s\n"%sv.id)
	elif sv.is_inversion():
		inv_out.write("\t".join([chr, str(sv.start-1), str(sv.end), sv.id])+"\n")

del_out.close()
ins_out.close()
inv_out.close()


print "> Classifying SV for formation mechanism"

rtWin=200 
rtGap=150
if len(sys.argv) > 3:
	minsize = int(sys.argv[3])
	if minsize<200:
		if minsize >= 50:
			rtWin=50
			rtGap=50
		else:
			raise Error("Minimum size of events should be 50bp")

os.system("perl %s/Breakpoint_Classification_Pipeline %s %s %s 0.5 200 50 85 %i %i %s %s > %s.log" % (sys.path[0],del_file,ins_file,inv_file,rtWin,rtGap,prefix,outdir,outpre))

print "> Postprocessing for mechanism classification"

mechs={}

txt=open(outdir+"/Output/%s.All_Output.txt"%prefix, "r")
txt.readline()
for line in txt:
	fields=line.split("\t")
	chr, start, end, event, mech = fields[0:5]
	unsure = fields[21].strip()
	seq = fields[-1].strip()
	event = event.title().strip()
	mech = mech.upper().strip()
	if unsure == "3" or unsure == "4":
		mech="NAHR_EXT"
	elif unsure == "5":
		mech="STEI_EXT"
	elif unsure == "6":
		mech="MTEI_EXT"
	if not chr.startswith("chr"):
		chr = "chr"+chr
	start=int(start) + (0 if event=="Insertion" else 1)
	end=int(end) + (1 if event=="Insertion" else 0)
	coord = (chr,start,end)
	mechs[coord]=mech
txt.close()

output=open(outpre+".gff", "w")
for sv in svs:
	coord = (sv.name, sv.start, sv.end)
	mech = mechs[coord]
	if mech is not None:
		sv.mech(mech)
	output.write(str(sv)+"\n")

output.close()

print "> Finished mechanism classification"

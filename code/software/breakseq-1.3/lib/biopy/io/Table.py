import csv

def parse(file, delimiter="\t", hasheader=True, textcols=[0], trimlist=True):
	table = {}
	handle = csv.reader(open(file), delimiter=delimiter, quotechar='"')
	header = handle.next() if hasheader else None
	for row in handle:
		for cindex in range(len(row)):
			if header is not None and cindex>=len(header): break
			val = row[cindex]
			col = cindex if header is None else header[cindex]
			if col not in table: table[col]=[]
			if cindex in textcols:
				table[col].append(val)
			elif val=="":
				table[col].append(None)
			else:
				table[col].append(float(val))
	if trimlist:
		for col in table:
			datalist=table[col]
			end=0
			while end-1 >= -len(datalist):
				if datalist[end-1] is None: end-=1
				else: break
			if end < 0:
				table[col]=datalist[0:end]
	return table


def calc_content(sequence):
	gc = 0
	at = 0
	for nt in sequence:
		nt = nt.upper()
		if nt == 'G' or nt == 'C':
			gc=gc+1
		elif nt == 'A' or nt == 'T':
			at=at+1
	total=float(gc+at)
	return gc/total if total > 0 else 0

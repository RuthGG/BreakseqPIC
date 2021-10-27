from biopy.io.CList import *

def calc_property(func, *paramsets):
	result=[]
	for params in paramsets:
		if params is not None and params != (None,):
			result.append(func(*params))
	if len(result)>1:
		avg=sum(result)/float(len(result))
		result.insert(0, avg)
	return result

def calc_freq(name, start, end, annotes):
	return len(find_annote(name, start, end, annotes))

def find_annote(name, start, end, annotes):
	i = Interval(name, start, end)
	return annotes.find(i)

import itertools

class Interval:

	def __init__(self, name, start, end=None, annotations=None, window=0):
		flanking=int(window/2)
		if end is None: end=start
		self.name=name
		self.start=min(start, end)
		self.end=max(start, end)
		self.start=self.start-flanking	
		self.end=self.end+flanking-(1 if start==end and flanking>0 else 0)
		self.size=self.end-self.start+1
		self.mid=(self.start+self.end+1)/2
		self.freq=1
		self.sub_intervals=[]
		self.annotations=[] if annotations is None else annotations

	def __str__(self):
		return self.name+":"+str(self.start)+"-"+str(self.end)+\
			("" if self.annotations is None or len(self.annotations)==0 else str(self.annotations))

	def __cmp__(self, other):
		if self.name != other.name:
			return cmp(self.name, other.name)
		elif self.start < other.start:
			return -1
		elif self.start > other.start:
			return 1
		else:
			if (self.end > other.end):
				return -1
			elif (self.end < other.end):
				return 1
			else:
				return 0

	def is_overlapping(self, other):
		return not (self.name!=other.name or self.start > other.end or self.end < other.start)

	def is_enclosing(self, other):
		return (self.name==other.name and self.start <= other.start and self.end >= other.end)

	def bin(self, size):
		bins=Intervals()
		for i in range(self.start, self.end+1, size):
			bin=Interval(self.name, i, min(i+size-1,self.end))
			bins.add(bin)
		return bins


class Intervals:

	def __init__(self, hashing=True):
		self.sets = {}
		self.hashing=hashing

	def __iter__(self):
		return itertools.chain.from_iterable(self.sets.itervalues())

	def load(self, file, icols=[0,1,2], acols=None, header=0, delim="\t", name_prefix=None, window=0):
		nregex = isinstance(delim, str)
		for line in open(file, "r"):
			try:
				if header>0:
					header=header-1
					continue
				line=line.strip(' \n\r')
				cols=line.split(delim) if nregex else delim.split(line)
				name=None if icols[0] is None else cols[icols[0]]
				if name_prefix is not None:
					name=name_prefix if name is None else name_prefix+name
				start=int(cols[icols[1]])
				end=int(cols[icols[2]])
				anno=None if acols is None or len(acols)<1 else [cols[a] for a in acols]
				self.add(Interval(name, start, end, anno, window))
			except ValueError:
				continue
		return self

	def add(self, *intervals):
		for interval in intervals:
			key=interval.name if self.hashing else None
			if key not in self.sets:
				self.sets[key]=[]
			self.sets[key].append(interval)
		return self

	def sort(self, allow_duplicate=True):
		for a in self.sets.values():
			self.__sort_intervals(a, allow_duplicate)
		return self

	def find(self, x):
		key=x.name if self.hashing else None
		if key in self.sets:
			return self.__find_intervals(self.sets[key], x)
		else:
			return []

	def __sort_intervals(self, a, allow_duplicate=True):
		a.sort()
		i=0
		while i < len(a)-1:
			if a[i].is_enclosing(a[i+1]):
				if allow_duplicate or a[i]!=a[i+1]:
					a[i].sub_intervals.append(a[i+1])
				del a[i+1]
				a[i].freq=a[i].freq+1
			else:
				if len(a[i].sub_intervals)>1:
					self.__sort_intervals(a[i].sub_intervals)
				i=i+1

	def __find_intervals(self, a, x, result=None):
		if result is None:
			result=[]
		lo = 0
		hi = len(a)
		while lo < hi:
			i = int((lo+hi)/2)
			if (a[i].is_overlapping(x)):
				self.__append_result(result, a[i], x)
				j=i+1
				while j < len(a):
					if (a[j].is_overlapping(x)):
						self.__append_result(result, a[j], x)
						j=j+1
					else:
						break
				j=i-1
				while j >= 0:
        	                        if (a[j].is_overlapping(x)):
						self.__append_result(result, a[j], x)
						j=j-1
                                	else:
                                        	break
				break
			else:
				if a[i] < x:
					lo = i+1
				elif a[i] > x: 
					hi = i
				else:
					return -1
		return result

	def __append_result(self, result, interval, x):
		result.append(interval)
		self.__find_intervals(interval.sub_intervals, x, result)

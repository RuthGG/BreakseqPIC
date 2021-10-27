import os, sys

def load(file='Config.txt', dir=None):
	config={}
	if dir is None:
		dir=os.path.dirname(sys.argv[0])
	path=os.path.join(dir,file)
	if not os.path.lexists(path):
		path=os.path.join(dir, "..", file)
	if not os.path.lexists(path):
		raise IOError("Configuration file '%s' not found."%path)
	for l in open(path):
		if not l.startswith("#") and l.strip()!="":
			kv = l.split("=", 1)
			config[kv[0].strip()]=eval(kv[1].strip())
	return config

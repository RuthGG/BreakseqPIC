#!/bin/env python

import sys

from sam import *

lib_len=int(sys.argv[1])
seed_len=int(sys.argv[2])
min_cover=int(sys.argv[3])

handle = open(sys.argv[4], "r") if len(sys.argv)>4 else sys.stdin

min_pos=lib_len/2 - (seed_len-min_cover) + 1
max_pos=lib_len/2 - min_cover + 1

for sam in SamParser(handle):
	if sam.pos < min_pos or sam.pos > max_pos:
		continue
	print sam

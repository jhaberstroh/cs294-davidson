#!/Users/markstrother/canopy/bin/python
import csv
import numpy as np
from itertools import permutations

def make_h(filename,h,N):
    with open(filename) as csvfile:
	readerator = csv.reader(csvfile)
	i=0;j=0
	for item in readerator:
	    h[i,j] = float(item[0])
	    if j == N-1:
		i+=1;j=0
	    else: j+=1

# FINISH THIS...
# You should compare with parse.py
# 
def make_teis(filename,teis,N):
    with open(filename) as f:
	for line in f:
	    if line[0] == '>':
		x = re.finall('[\d\.]+',line)
		arrayofeverything.append(x)
    for indexAndElement in arrayofeverything:
	abcd = indexAndElement[0:3]
	val=indexAndElement[4]
	for perm in permutations(abcd):

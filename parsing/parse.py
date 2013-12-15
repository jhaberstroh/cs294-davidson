#!/usr/bin/python
import re
import numpy as np

import os, sys, argparse 

# Tell the parser which file you're using and how many basis functions.
# The number of basis functions could be read using regular expressions,
# but that's some voodoo shit.
parser = argparse.ArgumentParser(description='Parses output from Psi4 SCF computation into a 2-D square matrix')
parser.add_argument('-i',type=str,nargs=1,help='The file containing the printed output from Psi4')
parser.add_argument('-b',default=[25],type=int,nargs=1,help='Number of basis functions')

args = parser.parse_args()

numBas = args.b[0]
filename = args.i[0]

arrayofeverything=[]

with open(filename) as f:
	for line in f:
		if line[0] == '>':
			#print line,
			x = re.findall('[\d\.]+', line)
			arrayofeverything.append(x)

def indexStripper(indexAndMatrixEl):
	mu = int(indexAndMatrixEl[0])
	nu = int(indexAndMatrixEl[1])
	rho = int(indexAndMatrixEl[2])
	sigma = int(indexAndMatrixEl[3])
	return mu,nu,rho,sigma

def indexConverter(mu,nu,rho,sigma):
	i = mu*numBas + nu
	j = rho*numBas + sigma
	return i,j

def permutations(mu,nu,rho,sigma):
	# This can be hella optimized, but I figure we only
        # need to do it once so why bother
        # The thing only takes a few seconds anyway.
	symmIndices=[]
	symmIndices.append(indexConverter(mu,nu,rho,sigma))
	symmIndices.append(indexConverter(nu,mu,rho,sigma))
	symmIndices.append(indexConverter(mu,nu,sigma,rho))
	symmIndices.append(indexConverter(nu,mu,sigma,rho))
	symmIndices.append(indexConverter(rho,sigma,mu,nu))
	symmIndices.append(indexConverter(sigma,rho,mu,nu))
	symmIndices.append(indexConverter(rho,sigma,nu,mu))
	symmIndices.append(indexConverter(sigma,rho,nu,mu))
	return symmIndices

N = numBas*numBas
matrix=np.zeros((N,N))

for l in arrayofeverything:
	mu,nu,rho,sigma = indexStripper(l)
        greeklist = permutations(mu,nu,rho,sigma)
	for p in greeklist:
		i,j = p
		matrix[i,j] = float(l[4])

# Prints the upper 10x10 square of a giant ass matrix
# I have no idea how to store the matrix in memory.
for i in range(0,10):
	for j in range(0,10):
		print matrix[i,j],
		print '\t',
		if(matrix[i,j]==0.0): print '\t',
	print '\n'





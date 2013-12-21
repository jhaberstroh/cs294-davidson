#!/Users/markstrother/canopy/bin/python


import os, sys, argparse 
import struct

# Just a little helper routine, to be used later
def diffList(lst1,lst2):
    returnList=[]
    for elem in lst1:
	if elem not in lst2: returnList.append(elem)
    for elem in lst2:
	if elem not in lst1: returnList.append(elem)
    return returnList


# Tell the parser which file you're using and how many basis functions.
# The number of basis functions could be read using regular expressions,
# but that's some voodoo shit.
#parser = argparse.ArgumentParser(description='Parses output from Psi4 SCF computation into a 2-D square matrix')
#parser.add_argument('-i',required=True,type=str,nargs=1,help='The file containing the printed output from Psi4')
#parser.add_argument('-b',required=True,default=[25],type=int,nargs=1,help='Number of basis functions (could be read from output, Psi4 file, but it\'s easeir this way.')

#args = parser.parse_args()
#
#numBas = args.b[0]
#filename = args.i[0]



from itertools import combinations
import numpy as np
import scipy as sc
K=1
N=2


hSpat = np.zeros((K,K))
hSpin = np.zeros((2*K,2*K))
# Going from spatial to spin OEIs is easy; < a | b > = 0
for i in range(0,K):
    for j in range(0,K):
	hSpin[2*i,2*j] = hSpat[i,j]
	hSpin[2*i+1,2*j+1] = hSpat[i,j]

# Going from spation to spin TEIs
teiSpat = np.zeros((K,K,K,K))
teiSpin = np.zeros((2*K,2*K,2*K,2*K))
for i in range(0,K):
    for j in range(0,K):
	for k in range(0,K):
	    for l in range(0,K):
		teiSpin[2*i,2*j,2*k,2*l]=teiSpat[i,j,k,l]
		teiSpin[2*i,2*j,2*k+1,2*l+1]=teiSpat[i,j,k,l]
		teiSpin[2*i+1,2*j+1,2*k,2*l]=teiSpat[i,j,k,l]
		teiSpin[2*i+1,2*j+1,2*k+1,2*l+1]=teiSpat[i,j,k,l]

teiDblBar = np.zeros((2*K,2*K,2*K,2*K))
for i in range(0,2*K):
    for j in range(0,2*K):
	for k in range(0,2*K):
	    for l in range(0,2*K):
		teiDblBar[i,j,k,l] = teiSpin[i,j,k,l]-teiSpin[i,j,l,k]
lst=range(0,2*K)
lstOcc=lst[0:N]
lstVirt=lst[N:2*K]

print "Occupied: "+str(lstOcc)
print "Virtual: "+str(lstVirt)

if K < N: nmax = K
elif N <= K: nmax = N

cmbOcc = []
cmbVirt = []

determinants = []
determinants.append(lstOcc)
for i in range(0,nmax):
    for elem in combinations(lstOcc,i+1):
        for elemVirt in combinations(lstVirt,i+1):
	    lstCopy = lst[0:N]
            for j in range(0,len(elem)):
                lstCopy[elem[j]] = elemVirt[j] 
  	    determinants.append(lstCopy)

Ndet=0                   
for det in determinants:
    print det    
    Ndet+=1
# Is there no better way to get binom coefficient??
# Scipy fail...
FullCIMatrix = np.zeros((Ndet,Ndet))
print "Number of determinants: "+str(Ndet)

# SINGLE ELECTRON OPERATOR, diagonals
for i in range(0,Ndet):
    print "CImatrixElem["+str(i)+","+str(i)+"]: "
    for m in determinants[i]:
        print "h["+str(m)+","+str(m)+"] +"
        FullCIMatrix[i,i] += hSpin[m,m] 

# SINGLE ELECTRON OPERATOR, off diagonal
Noff=0
for i in range(0,Ndet):
    for j in range(0,Ndet):
        if i == j: continue
	ndiff=0
	mSave=[]
	pSave=[]
        for m in range(0,N):
	    if determinants[i][m] not in determinants[j]:
		ndiff+=1
		# if ndiff > 1 

        if ndiff==1:
	    print "CIMatrixElem["+str(i)+","+str(j)+"]: "
	    mSave = diffList(determinants[i],determinants[j])
            print determinants[i],determinants[j]
	    print mSave[0],mSave[1]
	    FullCIMatrix[i,j] = hSpin[mSave[0],mSave[1]]
	    Noff+=1
#	elif ndiff==2 and (mSave[0]==pSave[1]):
#	    print "CIMatrixElem["+str(i)+","+str(j)+"]: "
#	    print determinants[i],determinants[j]
#	    print mSave[1],pSave[0]
#	    FullCIMatrix[i,j] = hSpin[mSave[1],pSave[0]]
#	    Noff+=1
#	elif ndiff==2 and (mSave[1]==pSave[0]):
#	    print "CIMatrixElem["+str(i)+","+str(j)+"]: "
#	    print determinants[i],determinants[j]
#	    print mSave[0],pSave[1]
#	    FullCIMatrix[i,j] = hSpin[mSave[0],pSave[0]]
#	    Noff+=1

print "Number of off diagonals: "+str(Noff)
for i in range(0,Ndet):
    for j in range(0,Ndet):
	print FullCIMatrix[i,j]

# DOUBLE ELECTRON OPERATOR, diagonals
for i in range(0,Ndet):
    for m in determinants[i]:
	for n in determinants[i]:
	    FullCIMatrix[i,i] += .5*teiDblBar[m,n,m,n]
	    # aka, < K | O2 | K > = sum sum <m n || m n>

# DOUBLE ELECTRON, off-diagonals
for i in range(0,Ndet):
    for j in range(0,Ndet):
	if i==j: continue
	ndiff=0
	for m in range(0,N):
	    if determinants[i][m] not in determinants[j]:
		ndiff+=1
		
	if ndiff==1:
	    mSave = diffList(determinants[i],determinants[j])
	    assert( len(mSave) == 2)
	    m = mSave[0]
	    p = mSave[1]
	    for n in determinants[i]:
		if n == m: continue
		FullCIMatrix[i,j] += teiDblBar[m,n,p,n]

	elif ndiff==2:
	    mSave = diffList(determinants[i],determinants[j])
	    assert( len(mSave) == 4)
	    m = mSave[0]
	    n = mSave[1]
	    p = mSave[2]
	    q = mSave[3]   
	    FullCIMatrix[i,j] = teiDblBar[m,n,p,q]

### FullCIMatrix is now done.
##   Dump to file.#

print "Full CI Matrix dimensions:",FullCIMatrix.shape
file_out = "derpitydoo.dat"
with open(file_out,'wb') as f2:
	f2.write("MOformat_");

	f2.write(struct.pack('I',FullCIMatrix.shape[0]));
	f2.write(struct.pack('I',8));
	for i in xrange(0,Ndet):
		for j in xrange(0,Ndet):
			f2.write(struct.pack('d',FullCIMatrix[i,j]))


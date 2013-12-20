import numpy as np
import scipy.linalg as la
import scipy.linalg.lapack as lapack
import scipy.sparse.linalg as sla
import scipy as sp
import random
import copy


#	Returns approximate (err, eval, evec, perp) for the m x m submatrix of A
def JDLoop(A, m, v0, V, tol, verbose = False):
	print "m =",m
	assert(type(A) == type(np.zeros(1)))
	assert(type(v0) == type(np.zeros(1)))
	n = A.shape[0]
	assert(v0.shape[0]== n and v0.shape[1]== 1)
	assert(V.shape[0] == n and V.shape[1] == m-1)
	
	t = copy.copy(v0)
	k = 0
	while k < .25:
		v0 = copy.copy(t)
		#print "\tT0:",t
		for i in xrange(m-1):
			#print np.dot(v0[:,0],V[:,i])
			#print "\tVECT:", V[:,i]
			t[:,0] -= np.dot(v0[:,0].T,V[:,i]) * V[:,i]
		#print "\tTF:", t
		#print "Norm v0, iteration",m,":",la.norm(v0)
		#print "Norm t:",la.norm(t)
		k = la.norm(v0) / la.norm(t)
		#print "Orthogonalization rescale kappa: ", k
	v0 = t / la.norm(t)

	Vm = np.append(V,v0,1)

	#print Vm.T.shape, A.shape, Vm.shape
	M = np.dot(Vm.T, np.dot(A, Vm))

	th, s = la.eigh(M, eigvals=(m-1,m-1))

	#print "S shape:", s.shape
	u = np.dot(Vm,s)
	r = np.dot(A, np.dot(Vm, s)) - th * u
	if la.norm(r) < tol:
		return la.norm(r), th, u, t

	diag = A.diagonal()
	diag.shape = (n,1)
	#print "A diag: ",diag
	I = np.eye(n)
	P = I - np.dot(u, u.T)
	A2 = A - (th * I)
	MAT = np.dot(P, np.dot(A2,P))

	#print MAT.shape, r.shape
	x, info = sla.minres(MAT,-r)
	if verbose:
		print "A:"
		print A
		print "r:"
		print r
		print "x:"
		print x
	
	t = x
	t /= la.norm(t)
	t.shape = (n,1)

	if verbose:
		print "PAPt = ",np.dot(P, np.dot(A2, np.dot(P, t)))
		print "U: ", u
		print "T: ", t
		print "R: ", r
	print "\t\t\t residual", la.norm(r)
	print "\t\t\t t dot u: ", np.dot(t.T,u)
	#assert(np.dot(t.T,u) <= 1e-5)

	return la.norm(r), th, u, t


def JDRound(A,v0, maxM=0, verbose = False):
	n = A.shape[0]
	V = np.zeros((n,0))
	if verbose:
		print "INPUT :::"
		print v0
	err = 1
	tol = 1e-5
	m = 1

	lam = 0
	if maxM == 0:
		maxM = n
	while err > tol and m <= maxM:
		err, lam, v0, perp = JDLoop(A,m,v0,V, tol)
		err = err / m
		V = np.append(V,perp,1)
		if verbose:
			print "Guess",m,", error,",err,":::"
			print v0
		m += 1
	print "Largest eigenvalue: ", lam
	if verbose:
		print "\tCorresponding evec: ", v0
	print "Final error: ", err
	return err, lam, v0
	

def JDRoutine(A, maxM = 0, maxRound = 10, verbose = False):
	n = A.shape[0]
	v0 = np.zeros((n,1))
	for i in xrange(n):
		v0[i] = random.uniform(-1,1)
	v0 /= la.norm(v0)

	tol = 1e-5
	err = tol + 1
	lam = 0
	rounds = 0
	while (err > tol and rounds < maxRound):
		print "Round ", rounds
		err, lam, v0 = JDRound(A,v0, maxM, verbose)
		rounds += 1


def main():
	#random.seed(8675309)
	#MT = np.zeros((3,3))
	#MT[0,0] = 1
	#MT[1,0] = -1
	#MT[0,1] = -1
	#MT[1,1] = 1
	#MT[2,2] = 2
	#JDRoutine(MT)
	

	#print "\n\n\n\n\n\n\nRandom large matrix:\n\n\n\n\n\n\n"
	#random.seed(90210)
	#N = 5
	#MT = np.zeros((N,N))
	#for i in xrange(200):
	#	pos1 = random.randint(0,N-1)
	#	pos2 = random.randint(0,N-1)
	#	num = random.uniform(-1,1)
	#	MT[pos1,pos2] = num
	#	MT[pos2,pos1] = num
	#JDRoutine(MT, 3, 10, False)
	
	
	print "\n\n\n\n\n\n\nRandom large matrix:\n\n\n\n\n\n\n"
	random.seed(90210)
	N = 25
	MT = np.zeros((N,N))
	for i in xrange(200):
		pos1 = random.randint(0,N-1)
		pos2 = random.randint(0,N-1)
		num = random.uniform(-1,1)
		MT[pos1,pos2] = num
		MT[pos2,pos1] = num
	JDRoutine(MT, 3, 100, False)

	
if __name__ == "__main__":
	main()

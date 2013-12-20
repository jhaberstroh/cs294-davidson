import numpy as np
import scipy.linalg as la
import scipy.linalg.lapack as lapack
import scipy.sparse.linalg as sla
import scipy as sp
import random
import copy


#	Returns approximate (err, eval, evec, perp) for the m x m submatrix of A
def JDLoop(A, m, perp, V, tol, verbose = False):
	print "m =",m
	assert(type(A) == type(np.zeros(1)))
	assert(type(perp) == type(np.zeros(1)))
	n = A.shape[0]
	assert(perp.shape[0]== n and perp.shape[1]== 1)
	assert(V.shape[0] == n and V.shape[1] == m-1)
	
	t = copy.copy(perp)
	k = 0
	while k < .001:
		perp = copy.copy(t)
		#print "\tT0:",t
		for i in xrange(m-1):
			print "\tDot product vm * vi for i=",i,"is (before removal): ",np.dot(V[:,i],perp[:,0],)
			#print "\tVECT:", V[:,i]

			t[:,0] -= (np.dot(V[:,i].T,t[:,0]) * V[:,i])
			print "\tDot product t  * vi for i=",i,"is (after removal):  ",np.dot(V[:,i],t[:,0].T,)
		#print "\tTF:", t
		#print "Norm v0, iteration",m,":",la.norm(v0)
		#print "Norm t:",la.norm(t)
		k = la.norm(t) / la.norm(perp)
		print "\t\tOrthogonalization rescale kappa: ", k
	perp= copy.copy(t)
	perp /= la.norm(t)
	vm = perp
	for i in xrange(m-1):
		print "\t\tPOST Dot product t  * vi for i=",i,"is:  ",np.dot(t[:,0].T,V[:,i])

	V = np.append(V,perp,1)

	M = np.dot(V.T, np.dot(A, V))

	th, s = la.eigh(M, eigvals=(m-1,m-1))

	u = np.dot(V,s)
	r = np.dot(A, np.dot(V, s)) - th * u
	if la.norm(r) < tol:
		print "\t\t\t residual", la.norm(r)
		print "\t\t\t SUCCESS!"
		return la.norm(r), th, u, t, vm

	diag = A.diagonal()
	diag.shape = (n,1)
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

	return la.norm(r), th, u, t, vm


def JDRound(A,v0, maxM=0, verbose = False):
	n = A.shape[0]
	V = np.zeros((n,0))
	V = np.append(V,v0,1)
	if verbose:
		print "INPUT :::"
		print v0
	err = 1
	tol = 1e-5
	m = 2
	 
	perp = np.zeros((n,1))
	for i in xrange(n):
		perp[i] = random.uniform(-1,1)
	perp /= la.norm(perp)

	lam = 0
	if maxM == 0:
		maxM = n
	while err > tol and m <= maxM:
		x = JDLoop(A,m,perp,V, tol)
		print "LENGTH OF X:", len(x)
		err, lam, v0, perp, vm = JDLoop(A,m,perp,V, tol)
		err = err / m
		V = np.append(V,vm,1)
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
	JDRoutine(MT, 4, 5, False)

	
if __name__ == "__main__":
	main()

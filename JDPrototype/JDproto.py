import numpy as np
import scipy.linalg as la
import scipy as sp
import random
import copy


#	Returns approximate (err, eval, evec, perp) for the m x m submatrix of A
def JDLoop(A, m, v0, V):
	assert(type(A) == type(np.zeros(1)))
	assert(type(v0) == type(np.zeros(1)))
	n = A.shape[0]
	assert(v0.shape[0]== n and v0.shape[1]== 1)
	assert(V.shape[0] == n and V.shape[1] == m-1)
	
	t = copy.copy(v0)
	for i in xrange(m-1):
		t[:,0] -= np.dot(v0[:,0],V[:,i]) * V[:,i]
	#print "Norm v0, iteration",m,":",la.norm(v0)
	#print "Norm t:",la.norm(t)
	k = la.norm(t) / la.norm(v0)
	#print "Orthogonalization rescale kappa: ", k
	v0 = t / la.norm(t)

	Vm = np.append(V,v0,1)

	#print Vm.T.shape, A.shape, Vm.shape
	M = np.dot(Vm.T, np.dot(A, Vm))

	th, s = la.eigh(M, eigvals=(m-1,m-1))

	#print "S shape:", s.shape
	u = np.dot(Vm,s)
	r = np.dot(A, np.dot(Vm, s)) - th * u

	diag = A.diagonal()
	diag.shape = (n,1)
	#print "A diag: ",diag
	t = (1/(diag-th)) * r
	t /= la.norm(t)
	#print "R: ", r
	#print "T: ", t

	return la.norm(r), th, u, t


def JDRound(A,v0, verbose = False):
	n = A.shape[0]
	V = np.zeros((n,0))
	if verbose:
		print "INPUT :::"
		print v0
	err = 1
	tol = 1e-5
	m = 1

	lam = 0
	while err > tol and m <= n:
		err, lam, v0, perp = JDLoop(A,m,v0,V)
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
	

def JDRoutine(A, verbose = False):
	n = A.shape[0]
	v0 = np.zeros((n,1))
	for i in xrange(n):
		v0[i] = random.uniform(-1,1)
	v0 /= la.norm(v0)

	tol = 1e-5
	err = tol + 1
	lam = 0
	rounds = 0
	while (err > tol and rounds < 10):
		print "Round ", rounds
		err, lam, v0 = JDRound(A,v0,verbose)
		rounds += 1


def main():
	random.seed(8675309)
	MT = np.zeros((3,3))
	MT[0,0] = 1
	MT[1,0] = -1
	MT[0,1] = -1
	MT[1,1] = 1
	MT[2,2] = 2
	JDRoutine(MT)
	

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
	JDRoutine(MT, False)

	
if __name__ == "__main__":
	main()

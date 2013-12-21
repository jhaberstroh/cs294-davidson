#include <vector>
#include <cmath>
#include "LargeMatrix.H"
#include "Real.H"

//#include "gmres.H"
//

//template <class Vector>
//class EyePrecondition
//{
//	public:
//		EyePrecondition() {}
//
//		Vector solve(const Vector& in) const
//			{
//				return in;
//			}
//};


//class SimpleMatrix
//{
//	private:
//		Real** arr;
//		int m_size;
//	
//	public:
//		SimpleMatrix(int a_size)
//			{
//				m_size = a_size;
//				arr = new Real*[m_size];
//				for (int i = 0 ; i < m_size ; i++)
//					{	
//						arr[i] = new Real[m_size];
//					}
//			}
//
//		~SimpleMatrix()
//			{
//				for (int i = 0 ; i < m_size ; i++)
//					{
//						delete [] arr[i];
//					}
//				delete[] arr;
//			}
//
//		Real& operator ()(int row, int col)
//			{
//				return arr[row][col];
//			}
//};
//
//
//Real norm(std::vector<Real> v)
//{
//	Real sum = 0;
//	for (int i = 0 ; i < v.size() ; i++)
//		{
//			sum += v[i] * v[i];
//		}
//	return std::sqrt(sum);
//}
//
//std::vector<double> minus(const std::vector<double>& p, const std::vector<double>& m)
//{
//	assert(p.size() == m.size());
//	std::vector<double> c(p.size());
//
//	for (int i = 0 ; i < p.size() ; i++)
//		{
//			c[i] = p[i] - m[i];
//		}
//	return c;
//
//}

//
//int main()
//	{
//		EyePrecondition<std::vector<Real> > precond;
//
//		Mtx::Format fmt = Mtx::FMT_OPT_TEI;
//		int rows = 25 * 25;
//		std::string mtx_name = "../parsing/matx_out_test.mot";
//		Mtx::LargeMatrix mt(mtx_name, rows, fmt, Mtx::FRead|Mtx::FWrite, 100000);
//
//		SimpleMatrix H(rows);
//
//		std::vector<Real> b(rows);
//		std::vector<Real> x(rows);
//
//		int m = 0;
//		int maxiter = 10;
//		Real q = 1e-6;
//
//		int result = GMRES(mt, x, b, precond, H, m, maxiter, q);
//	}




#include <tminres.hpp>
#include "SimpleVector.hpp"
#include <cmath>
#include <fstream>


class Preconditioner
{
public:
	//! Y = M\X
	virtual void Apply(const SimpleVector & X, SimpleVector & Y) const = 0;
};

class Operator
{
public:

	Operator(Mtx::LargeMatrix a_mtx) :
		m_mtx(a_mtx)
	{
	}

	//! Y = A*X;
	void Apply(const SimpleVector & X, SimpleVector & Y) const
	{
		std::vector<double> x(m_mtx.size());
		for (int i = 0 ; i < m_mtx.size() ; i++)
			{
				x[i] = X[i];
			}
		std::vector<double> y = m_mtx * x;
		for (int i = 0 ; i < m_mtx.size() ; i++)
			{
				Y[i] = y[i];
			}
	}
	
	void Print(std::ostream & os)
	{
		os<< "LargeMatrix class" <<std::endl;
	}

private:
	Mtx::LargeMatrix m_mtx;
};


/*!
 * \example SerialExample/ex1.cpp
 *
 * A simple example of MINRES() usage without preconditioner.
 */

int main()
{

	Mtx::Format fmt = Mtx::FMT_OPT_TEI;
	int rows = 25 * 25;
	std::string mtx_name = "../parsing/matx_out_test.mot";
	Mtx::LargeMatrix mt(mtx_name, rows, fmt, Mtx::FRead|Mtx::FWrite, 100000);

	//(1) Define the size of the problem we want to solve
	int size(rows);
	//(2) Define the linear operator "op" we want to solve.
	Operator op(mt);
	//(3) Define the exact solution (at random)
	SimpleVector sol(size);
	sol.Randomize( 1 );

	//(4) Define the "rhs" as "rhs = op*sol"
	SimpleVector rhs(size);
	op.Apply(sol, rhs);
	double rhsNorm( sqrt(InnerProduct(rhs,rhs)) );
	std::cout << "|| rhs || = " << rhsNorm << "\n";

	//(5) We don't use any preconditioner. Let prec be a null pointer.
	Preconditioner * prec = NULL;

	//(6) Use an identically zero initial guess
	SimpleVector x(size);
	x = 0;

	//(7) Set the minres parameters
	double shift(0);
	int max_iter(1000);
	double tol(1e-6);
	bool show(true);

	//(8) Solve the problem with minres
	MINRES(op, x, rhs, prec, shift, max_iter, tol, show);

	//(9) Compute the error || x_ex - x_minres ||_2
	subtract(x, sol, x);
	double err2 = InnerProduct(x,x);
	std::cout<< "|| x_ex - x_n || = " << sqrt(err2) << "\n";

	std::ofstream fid("ex1.m");
	op.Print(fid);
	fid<< "rhs = [";
	for(int i(0); i<size-1; ++i)
		fid<<rhs[i] <<"; ";
	fid<<rhs[size-1] <<"]; \n";
	fid<<"[ x, istop, itn, rnorm, Arnorm, Anorm, Acond, ynorm ] = minres(Op, rhs, [], 0, true, false, 100, 1e-6);\n";


	return 0;
}

#ifndef __TMINRES_H_DEFINED__
#define __TMINRES_H_DEFINED__
#include <cmath>
#include "SimpleVector.hpp"

class Preconditioner
{
public:
	//! Y = M\X
	virtual void Apply(const SimpleVector & X, SimpleVector & Y) const = 0;
};

class LargeMtxOperator
{
	public:

		LargeMtxOperator(Mtx::LargeMatrix a_mtx) :
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

void RemoveProjection(std::vector<double>& x, std::vector<double> const& u)
	{
		double u_len = 0;
		double projection = 0;
		for (int i = 0 ; i < x.size() ; i++)
			{
				u_len += u[i] * u[i];
				projection += x[i] * u[i];
			}
		u_len = std::sqrt(u_len);
		projection /= u_len;

		for (int i = 0 ; i < x.size() ; i++)
			{
				x[i] -= u[i] * projection;
			}
	}

class LargeMtxProjectionOperator
{

	public:
		LargeMtxProjectionOperator(const Mtx::LargeMatrix& a_mtx, const Mtx::colVect& a_u):
			m_mtx(a_mtx), m_u(a_u)
		{
		}

		void Apply(const SimpleVector & X, SimpleVector & Y) const
			{
				std::vector<double> x(m_mtx.size());
				for (int i = 0 ; i < m_mtx.size() ; i++)
					{
						x[i] = X[i];
					}
				
				RemoveProjection(x,m_u);
				std::vector<double> y = m_mtx * x;
				RemoveProjection(y,m_u);

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
		Mtx::colVect m_u;


};


//class MtxOperator
//{
//	public:
//		MtxOperator(double** a_mtx, int size):
//			m_mtx(a_mtx)
//		{
//		}
//
//		void Apply(const SimpleVector & X, SimpleVector & Y) const
//			{
//				std::vector<double> x(size);
//				for (int i = 0 ; i < m_mtx.size() ; i++)
//					{
//						x[i] = X[i];
//					}
//
//				// NEED PROPER MULTIPLICATION
//				std::vector<double> y = x;
//
//				for (int i = 0 ; i < m_mtx.size() ; i++)
//					{
//						Y[i] = y[i];
//					}
//			}
//
//		void Print(std::ostream & os)
//			{
//				os<< "LargeMatrix class" <<std::endl;
//			}
//
//	private:
//		double** m_mtx;
//		int size;
//
//};

#endif//__TMINRES_H_DEFINED__

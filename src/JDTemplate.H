#ifndef __JDTEMPLATE_H_DEFINED__
#define __JDTEMPLATE_H_DEFINED__

#include "Real.H"
#include "LargeMatrix.H"
#include <gsl/gsl_rng.h>
#include <vector>
#include <cmath>
#include <lapacke.h>
#include "LargeTminres.H"
#include "tminres.hpp"
#include "SimpleVector.hpp"

#define RESTART 10
#define MAXROUND 10
#define VERBOSE false
#define TOL 1e-6


template <class Vector>
Real dot(int size, const Vector& a, const Vector& b)
	{
		Real sum = 0;
		for (int i = 0 ; i < size ; i++)
			{
				sum += a[i] * b[i];
			}
		return sum;
	}

template <class Vector>
Real norm(int size, const Vector& v)
	{
		Real sum = 0;
		for (int i = 0 ; i < size ; i++)
			{
				sum += v[i] * v[i];
			}
		return std::sqrt(sum);
	}


template <class Vector>
Real normalize(int size, Vector& v)
	{
		Real norm_len = norm(size, v);
		for (int i = 0 ; i < size ; i++)
			{
				v[i] /= norm_len;
			}
	
	}



template <class Vector>
void GrahamSchmidt(int n, int m, Vector& perp, std::vector<Vector>& V)
	{
		
		double kappa = 0;
		while (kappa < .01)
			{ 
				for (int i = 0 ; i < m-1 ; i++)
					{
						Real v_proj_i = 0;
						for (int j = 0 ; j < n ; j++)
							{
								v_proj_i += V[i][j] * perp[j];
							}
						for (int j = 0 ; j < n ; j++)
							{
								perp[j] -= v_proj_i * V[i][j];
							}
					}
				kappa = norm(n, perp);
				normalize(n, perp);
			}
		V.push_back(perp);
	}


// A			Operator to diagonalize
// m			Current round of diagonalizatoin
// perp		Appended to V and used as output
// V 			Array of vectors
template <class Matrix, class Vector>
void JDLoop(const Matrix& A, int m, Vector& perp, std::vector<Vector>& V, Vector& eig_vec, double& eig_val, double& residue)
	{
		int n = A.size();
		GrahamSchmidt(n, m, perp, V);

		Real M[m * m];
		for (int i = 0 ; i < m ; i++)
			{
				for (int j = 0 ; j < m ; j++)
					{
						M[i + m*j] = dot(m, V[i], (A * V[j]));
					}
			}

		lapack_int size_mtx = m;
		double eig_vals[m];
		int result = LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U', size_mtx, M, m, eig_vals);

		int max_index = m-1;

		//std::cout << "Eigenval"<<max_index<<": " << eig_vals[max_index]<<std::endl;
		//std::cout << std::endl;

		eig_val = eig_vals[max_index];
		Real* max_evec = M + (m * max_index);

		for (int j = 0 ; j < n ; j++)
			{
				eig_vec[j] = 0;
			}

		for (int i = 0 ; i < m ; i++)
			{
				for (int j = 0 ; j < n ; j++)
					{
						eig_vec[j] += max_evec[i] * V[i][j];
					}
			}

		Vector r = A * eig_vec;
		SimpleVector res_simple(n);
		SimpleVector perp_simple(n);
		for (int i = 0 ; i < n ; i++)
			{
				r[i] -= eig_vec[i] * eig_val;
				res_simple[i] = r[i];
				perp_simple[i] = 0;
			}

		std::cout << "Residue norm: " << norm(n,r) << std::endl;

		LargeMtxProjectionOperator op(A,eig_vec);
		Preconditioner * prec = NULL;

		double shift(0);
		int max_iter(1000);
		double tol(1e-6);
		bool show(false);
	
		MINRES(op, perp_simple, res_simple, prec, shift, max_iter, tol, show);

		for (int i = 0 ; i < n ; i++)
			{
				perp[i] = perp_simple[i];
			}

		std::cout << "\t\t\tt dot u: " << dot(n,perp,eig_vec) <<std::endl;
	}



template <class Matrix, class Vector>
void JDRound(const Matrix& A, Vector& eig_vec, double& eig_val, double& residue)
	{
		int n = A.size();
		std::vector<Vector> V(0);
		Vector perp(eig_vec);

		int m = 0;
		while (residue > TOL && m < RESTART)
			{

				m++;
				JDLoop(A,m,perp,V, eig_vec,eig_val,residue);

			}
	}



template <class Matrix>
void JDRoutine(const Matrix& A)
	{
		int n = A.size();
		
		Mtx::colVect eig_vec(n);
		gsl_rng * r;
		const gsl_rng_type * T;
		gsl_rng_env_setup();
		T = gsl_rng_default;
		r = gsl_rng_alloc(T);
		gsl_rng_set(r, 90210);

		for (int i = 0 ; i < n ; i++)
			{
				eig_vec[i] = gsl_rng_uniform(r);
			}

		gsl_rng_free(r);


		double residue = TOL + 1;
		Real eig_val = 0;
		int round_counter = 1;
		while (residue > TOL && round_counter <= MAXROUND)
			{
				std::cout << "Round : "<<round_counter <<std::endl;
				JDRound(A, eig_vec, eig_val, residue);
				round_counter ++;
			}

	}


#endif//__JDTEMPLATE_H_DEFINED__

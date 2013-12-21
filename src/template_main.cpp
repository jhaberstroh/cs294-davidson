#include "JDTemplate.H"
#include "LargeMatrix.H"
#include "lapacke.h"

#define N 5
#define LDA N

int main()
	{
		Mtx::Format fmt = Mtx::FMT_OPT_TEI;
		int rows = 25 * 25;
		std::string mtx_name = "../parsing/matx_out_test.mot";
		Mtx::LargeMatrix mt(mtx_name, rows, fmt, Mtx::FRead|Mtx::FWrite, 100000);

 //n = N, lda = LDA, info;	
 // 	double w[N];
 //       double a[LDA*N] = {
 //           1.96, -6.49, -0.47, -7.20, -0.65,
 //           0.00,  3.80, -6.39,  1.50, -6.34,
 //           0.00,  0.00, 4.17, -1.51, 2.67,
 //           0.00,  0.00, 0.00,  5.70, 1.80,
 //           0.00,  0.00, 0.00,  0.00, -7.10
 //       };


 //  printf( "LAPACKE_dsyev (row-major, high-level) Example Program Results\n" );
 //       /* Solve eigenproblem */
 //       info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, a, lda, w );
 //       /* Check for convergence */
 //       if( info > 0 ) {
 //               printf( "The algorithm failed to compute eigenvalues.\n" );
 //               exit( 1 );
 //       }
 //       /* Print eigenvalues */
 //       print_matrix( "Eigenvalues", 1, n, w, 1 );
 //       /* Print eigenvectors */
 //       print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
 //       exit( 0 );

		JDRoutine(mt);
	}

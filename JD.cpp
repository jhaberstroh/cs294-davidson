//g++ -framework Accelerate JD.cpp
//Follows the algorithm here: http://web.eecs.utk.edu/~dongarra/etemplates/node138.html
#include <iostream>
#include <vector> 
#include "Matrix/LargeMatrix.H"
#include <stdlib.h>
#include <Accelerate/Accelerate.h>
//#include "blas.h"

#define N 3

using namespace std;

void matrixVectorMultiply(double matrix[][N], double vector[], double result[]){
    int i,j;
    
    for(i=0; i<N; i++){
        result[i] = 0;
        for(j=0; j<N; j++){
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

void matrixVectorMultiply(std::vector <std::vector <double> > matrix, std::vector<double> vector, std::vector<double> result){
    int i,j;
    
    for(i=0; i<N; i++){
        result[i] = 0;
        for(j=0; j<N; j++){
            result[i] += matrix[i][j] * vector[j];
        }
    }
}

int delta(int i, int j){
    if(i == j) return 1;
    return 0;
}

int main(){
    double A[N][N];
    int LWORK = 2*4*N, INFO;
    double t[N], v[N][N], u[N], uA[N], r[N], M[N][N], M_lapack[N*N], vA[N][N], WR[N], WI[N], VL[N], VR[N*N], WORK[LWORK];
    double IminUU[N][N], IminUU_lapack[N*N], AminTheta[N][N], AminTheta_lapack[N*N];
    int i, j;
    double tmag;
    
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            M[i][j] = 0;
            v[i][j] = 0;
        }
    }

    
    srand (time(NULL));

    //initialize the random matrix to substitute for an electronic structure Hamiltonian
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            A[i][j] = rand() % 1000;
        }
    }
    
    //(1)
    //initialize the random eigenvector guess
    for(i=0; i<N; i++) {
        t[i] = rand() % 1000;
        //v[i][0] = t[i];
    }
    
    //(2)
    for(int m=1; m<N; m++){
        //(3)
        //orthogonalization
        for(i=1; i<m; i++){
            //(4)
            double vt = 0;
            for(j=0; j<N; j++) vt += v[j][i] * t[j];
            for(j=0; j<N; j++) t[j] -= vt * v[j][i];
        }//(5)
        //(6)
        //normalization
        tmag = 0;
        for(j=0; j<N; j++) tmag += t[j] * t[j];
        tmag = sqrt(tmag);
        for(j=0; j<N; j++) v[j][m-1] = t[j] / tmag;
        
        //check the new vector
        cout << "v's\n";
        for(i=0; i<N; i++){
            for(j=0; j<N; j++){
                cout << v[i][j] << " ";
            }
            cout << endl;
        }
        
        double vm[N];
        for(i=0;i<N;i++) vm[i] = v[i][m];
        matrixVectorMultiply(A, vm, vA[m-1]);
        //(7)
        for(i=1;i<=m;i++){
            //(8)
            M[i-1][m-1] = 0;
            for(j=0;j<N;j++){
                M[i-1][m-1] += v[j][i] * vA[m-1][j];
            }
        }//(9)
        
        //(10) calculate largest eigenpair of M with LAPACK
        char left = 'N';
        char right = 'V';
        int size = N;
        int LDVL = 1;
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                M_lapack[j+i*N] = M[j][i];
            }
        }
        //right eigenvectors are stored in the columns of VR in the same order as the corresponding eigevanlues in WR and WI (real,imaginary)
        dgeev_(&left,&right,&size,M_lapack,&size,WR,WI,VL,&LDVL,VR,&size,WORK,&LWORK,&INFO);
        cout << "eigenvalues\n";
        for(i=0;i<N;i++) cout << WR[i] << " ";
        cout << endl;
        cout << "eigenvectors\n";
        for(i=0;i<N;i++){
            for(j=0;j<N;j++) cout << VR[j*N + i] << " ";
            cout << endl;
        }
        
        //find the largest eigenvalue and its corresponding eigenvector
        double maxEigen = WR[0];
        int maxEigenLoc = 0;
        for(i=1;i<N;i++){
            if(WR[i] > maxEigen){
                maxEigen = WR[i];
                maxEigenLoc = i;
            }
        }
        
        //(11)
        for(i=0;i<N;i++){
            u[i] = 0;
            for(j=0;j<N;j++){
                u[i] += v[i][j] * VR[j*N + maxEigenLoc];
            }
        }
        
        //(12)
        for(i=0;i<N;i++){
            uA[i] = 0;
            for(j=0;j<N;j++){
                uA[i] += vA[i][j] * VR[j*N + maxEigenLoc];
            }
        }
        
        //(13)
        for(i=0;i<N;i++) r[i] = uA[i] - WR[maxEigenLoc] * u[i];
        //(14)
        double rmag = 0;
        for(i=0;i<N;i++) rmag += r[i] * r[i];
        rmag = sqrt(rmag);
        cout << "rmag: " << rmag << endl;
        //(15)
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                IminUU[i][j] = delta(i,j) - u[i] * u[j];
                AminTheta[i][j] = A[i][j] - WR[maxEigenLoc] * delta(i,j);
            }
        }
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                IminUU_lapack[j+i*N] = IminUU[j][i];
                AminTheta_lapack[j+i*N] = AminTheta[j][i];
            }
        }
        
        //(15)
        double C[N*N], D[N*N];
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,1,IminUU_lapack,N,AminTheta_lapack,N,0,C,N);
        cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,1,C,N,IminUU_lapack,N,0,D,N);
        
        int NRHS = 1;
        int IPIV[N];
        //r goes in as the residue, comes out as the new eigenvector guess t
        dgesv_(&size,&NRHS,D,&size,IPIV,r,&size,&INFO);
        for(i=0;i<N;i++) t[i] = r[i];
    }//(16)
}

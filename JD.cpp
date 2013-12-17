//g++ -framework Accelerate JD.cpp
//Follows the algorithm here: http://web.eecs.utk.edu/~dongarra/etemplates/node138.html
#include <iostream>
#include <vector> 
//#include "Matrix/LargeMatrix.H"
#include <stdlib.h>
#include <Accelerate/Accelerate.h>

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
    int i, j;
    int LWORK = 2*4*N, INFO;
    double t[N], v[N][N], u[N], uA[N], r[N], A[N][N], M[N][N], M_lapack[N*N], vA[N][N], WR[N], WI[N], VL[N], VR[N*N], WORK[LWORK];
    double IminUU[N][N], IminUU_lapack[N*N], AminTheta[N][N], AminTheta_lapack[N*N];
    double tmag;

    //initialize the random matrix to substitute for an electronic structure Hamiltonian
    //srand (time(NULL));
    srand(12345);
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            A[i][j] = rand() % 1000;
        }
    }
    
    cout << "Matrix we're solving:\n";
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    
    //(1)
    //initialize the random normalized eigenvector guess
    for(i=0; i<N; i++) {
        t[i] = rand() % 1000;
        //v[i][0] = t[i];
    }
    tmag = 0;
    for(j=0; j<N; j++) tmag += t[j] * t[j];
    tmag = sqrt(tmag);
    for(j=0;j<N;j++) t[j] = t[j] / tmag;
    
    cout << "Our random guess for the eigenvector: ";
    for(i=0;i<N;i++) cout << t[i] << " ";
    cout << endl;
    
    //iterate the whole process
    int maxIter = 1;
    for(int iter = 0; iter<maxIter; iter++){
        //reset the list of guesses, v, and the matrix M to zero at the beginning of each iteration
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                M[i][j] = 0;
                v[i][j] = 0;
            }
        }
        //cout << "iter: " << iter << endl;
        //(2)
        for(int m=1; m<=N; m++){
            cout << "t: ";
            for(i=0;i<N;i++) cout << t[i] << " ";
            cout << endl;
            for(j=0; j<N; j++) v[j][m-1] = t[j];
            //(3) - orthogonalization
            for(i=1; i<m; i++){
                //(4)
                double vt = 0;
                for(j=0; j<N; j++) vt += v[j][i-1] * t[j];
                cout << "vt: " << vt << endl;
                for(j=0; j<N; j++) v[j][m-1] -= vt * v[j][i-1];
            }//(5)
            //normalization
            tmag = 0;
            for(j=0; j<N; j++) tmag += v[j][m-1] * v[j][m-1];
            tmag = sqrt(tmag);
            for(j=0;j<N;j++) v[j][m-1] = v[j][m-1] / tmag;
           
            cout << "orthonormalized t: ";
            for(i=0;i<N;i++) cout << v[i][m-1] << " ";
            cout << endl;
            

            if(tmag == 0){
                cout << "Your eigenvector guess is the null vector.";
                abort();
            }
             //(6)
            
            
            //check the new vector
            cout << "v's\n";
            for(i=0; i<N; i++){
                for(j=0; j<N; j++){
                    cout << v[i][j] << " ";
                }
                cout << endl;
            }
            
            double vm[N], vAm[N];
            for(i=0;i<N;i++) vm[i] = v[i][m-1];
            matrixVectorMultiply(A, vm, vAm);
            for(i=0;i<N;i++) vA[i][m-1] = vAm[i];
            
            /*cout << "vector vA: ";
            for(i=0;i<3;i++) cout << vA[i][m-1] << " ";
            cout << endl;*/
            //(7)
            for(i=1;i<=m;i++){
                //(8)
                M[i-1][m-1] = 0;
                for(j=0;j<N;j++){
                    M[i-1][m-1] += v[j][i-1] * vA[m-1][j];
                }
            }//(9)
            
            /*cout << "matrix M\n";
            for(i=0;i<N;i++){
                for(j=0;j<N;j++){
                    cout << M[i][j] << " ";
                }
                cout << endl;
            }*/
            
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
            
            /*cout << "eigenvalues\n";
            for(i=0;i<N;i++) cout << WR[i] << " ";
            cout << endl;
            cout << "eigenvectors\n";
            for(i=0;i<N;i++){
                for(j=0;j<N;j++) cout << VR[j*N + i] << " ";
                cout << endl;
            }*/
            
            //find the largest eigenvalue and its corresponding eigenvector
            double maxEigen = WR[0] * WR[0];//want the maximum magnitude, don't care about the sign
            int maxEigenLoc = 0;
            //cout << "WR[0]*WR[0]: " << WR[0] * WR[0] << endl;
            //cout << "maxEigen: " << maxEigen << endl;
            for(i=1;i<N;i++){
                //cout << "WR[i] * WR[i]: " << WR[i] * WR[i] << endl;
                //cout << "maxEigen: " << maxEigen << endl;
                if(WR[i] * WR[i] > maxEigen){
                    maxEigen = WR[i] * WR[i];
                    maxEigenLoc = i;
                }
            }
            //cout << "maxEigenLoc: " << maxEigenLoc << endl;
            
            //(11)
            for(i=0;i<N;i++){
                u[i] = 0;
                for(j=0;j<N;j++){
                    u[i] += v[i][j] * VR[maxEigenLoc*N + j];
                }
            }
            
            /*cout << "vector u: ";
            for(i=0;i<N;i++) cout << u[i] << " ";
            cout << endl;*/
            
            //(12)
            for(i=0;i<N;i++){
                uA[i] = 0;
                for(j=0;j<N;j++){
                    uA[i] += vA[i][j] * VR[maxEigenLoc*N + j];
                }
            }
            /*cout << "vector uA: ";
            for(i=0;i<N;i++) cout << uA[i] << " ";
            cout << endl;*/
            
            //(13)
            for(i=0;i<N;i++) r[i] = uA[i] - WR[maxEigenLoc] * u[i];
            /*cout << "residual vector: ";
            for(i=0;i<N;i++) cout << r[i] << " ";
            cout << endl;*/
            //(14)
            double rmag = 0;
            for(i=0;i<N;i++) rmag += r[i] * r[i];
            rmag = sqrt(rmag);
            cout << "rmag: " << rmag << endl;
            if(rmag < 0.25){
                cout << "finished!\n";
                cout << "eigenvalue: " << WR[maxEigenLoc] << endl;
                cout << "eigenvector: " << endl;
                for(i=0;i<N;i++) cout << VR[maxEigenLoc*N + i] << " ";
                cout << endl;
                iter = maxIter;
                break;
            }
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
            /*cout << "matrix C: ";
            for(i=0;i<N*N;i++) cout << C[i] << " ";
            cout << endl;*/
            cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,N,N,N,1,C,N,IminUU_lapack,N,0,D,N);
            /*cout << "matrix D: ";
            for(i=0;i<N*N;i++) cout << D[i] << " ";
            cout << endl;*/
            
            int NRHS = 1;
            int IPIV[N];
            //r goes in as the residue, comes out as the new eigenvector guess t
            dgesv_(&size,&NRHS,D,&size,IPIV,r,&size,&INFO);
            /*cout << "new vector t: ";
            for(i=0;i<N;i++) cout << r[i] << " ";
            cout << endl;*/
            for(i=0;i<N;i++) t[i] = r[i];
        }//(16)
    }
}

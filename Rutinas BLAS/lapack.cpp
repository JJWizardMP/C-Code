//g++ lapack.cpp -o lapack -lf77blas && ./lapack

#include <iostream>

extern "C" {
	//lvl 1 Blas routine
	void saxpy_(int *N, float *SA, float *SX, int *INCX, float *SY, int *INCY);
	// lvl 2 Blas routine
	void sgemv_(char *TRANS, int *M, int *N, float *AlPHA, float *A, int *LDA, 
		float *X, int *INCX, float *BETA, float *Y, int *INCY);
	// lvl 3 Blas routine
	void sgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, float *ALPHA,
		float *A, int *LDA, float *B, int *LDB, float *BETA, float *C, int *LDC);
}

const int SIZE = 3;
void PrintVector(float *, int, char*);
void PrintMatrix(float [][SIZE], int, int, char*);

int main(){
	// Values needed for saxpy [ lvl 1 ]
	int n=3;
	float sa=2;
	float sx[SIZE] = {4,6,1}, sy[SIZE]= {9,7,6};
	int incx=1, incy=1;
	// Values needed for sgemv [ level 2]
	/* n=3, incx=1, incy=1*/
	int m=3;
	float alpha=3, beta=6;
	int lda=3;
	float A[][SIZE] = {{1,4,1},{1,6,4},{5,0,1}};
	float x[SIZE] = {1,3,2}, y[SIZE]= {2,4,1};
	// Values needed for sgemm [ lvl 3 ]
	/* m=3, n=3, lda=3, alpha=3, beta=6, matrix A*/
	int  k=3;
	int ldb=3, ldc=3;
	//float A[][SIZE] = {{1,4,1},{1,6,4},{5,0,1}};
	float B[][SIZE] = {{2,3,1},{1,5,1},{3,1,1}};
	float C[][SIZE] = {{5,2,2},{1,7,4},{1,4,1}};

	std::cout << "\nLevel 1 Blas saxpy routine:\n\n";
	std::cout << "SA: " << sa << "\n";
	PrintVector(sx, n, (char*)"SX: ");
	PrintVector(sy, n, (char*)"SY: ");
	// constant times a vector plus a vector.
 	saxpy_(&n, &sa, &sx[0], &incx, &sy[0], &incy);	
	PrintVector(sy, n, (char*)"\ny = a*x + y\nY: ");
	
	std::cout << "\nLevel 2 Blas segmv routine:\n\n";
	std::cout << "ALPHA: " << alpha << "\n";
	std::cout << "BETA: " << beta << "\n";
	PrintVector(x, m, (char*)"X: ");
	PrintVector(y, n, (char*)"Y: ");
	PrintMatrix(A, m, n, (char*)"A: ");
	// performs one of the matrix-vector operations
 	sgemv_((char*) "T", &m, &n, &alpha, &A[0][0], &lda, &x[0], &incx, &beta, &y[0], &incy);	
	PrintVector(y, n, (char*)"y = alpha*A*x+beta*y\nY: ");

	std::cout << "\nLevel 3 Blas sgemm routine:\n\n";
	std::cout << "ALPHA: " << alpha << "\n";
	std::cout << "BETA: " << beta << "\n";
	PrintMatrix(A, m, n, (char*)"A: ");
	PrintMatrix(B, m, n, (char*)"B: ");
	PrintMatrix(C, m, n, (char*)"C: ");
	// performs one of the matrix-matrix operations
 	sgemm_((char*)"T", (char*)"T", &m, &n, &k, &alpha, &A[0][0], &lda, 
		&B[0][0], &ldb, &beta, &C[0][0], &ldc);	
	std::cout << "C := alpha*op( A )*op( B ) + beta*C" << std::endl;
	PrintMatrix(C, m, n, (char*)"C: ");
}

void PrintVector(float *SA, int N, char *msg){
	std::cout << msg << "[ ";
	for(int i=0; i<N; i++){
		std::cout << SA[i] << " ";
	}
	std::cout << "]" << std::endl;
}
void PrintMatrix(float A[][SIZE], int M, int N, char *msg){
	std::cout << msg << std::endl;
	for(int i=0; i<M; i++){
		std::cout << "[ ";
		for(int j=0; j<N; j++){
			std::cout << A[i][j] << " ";
		}
		std::cout << "]\n";
	}
	std::cout << std::endl;
}

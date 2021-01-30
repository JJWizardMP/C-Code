//g++ dot.cpp -o dot -lf77blas

#include <iostream>
#include <stdio.h>
#include <iomanip>
using namespace std;


extern "C" {
	//lvl 1 Blas routine
	float snrm2_(int *N, float *X, int *INCX);
	void scopy_(int *N, float *SX, int *INCX, float *SY, int *INCY);
	void saxpy_(int *N, float *SA, float *SX, int *INCX, float *SY, int *INCY);
	// lvl 2 Blas routine
	void sgemv_(char *TRANS, int *M, int *N, float *AlPHA, float *A, int *LDA, 
		float *X, int *INCX, float *BETA, float *Y, int *INCY);
	// lvl 3 Blas routine
	void sgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K, float *ALPHA,
		float *A, int *LDA, float *B, int *LDB, float *BETA, float *C, int *LDC);
}

class Matrix{
	public:
	float **M;
	int m, n;
	// Constructors
	// 1: Create a zeros dinamic matrix
	Matrix(int, int);
	// 2: Read from MTX
	Matrix(char*);
	// Getters
	float **getMatrix(){ return M; };
	int getRows(){ return m; };
	int getColumns(){ return n; }
	// Methods
	float **NewMatrix(int, int);
	void Zeros(float ***, int, int);
	void Read_MatrixMTX(float ***, int *, int *, char*);	
};
//Functions
void PrintMatrix(float **, int, int, char*);

int main(){

	Matrix A((char*)"./A.mtx");
	Matrix B((char*)"./B.mtx");
	Matrix C(2,2);

	char *trans = (char*) "T";
	int k=1, lda=3;
	float sa=1;
	sgemm_(trans, trans, &A.m, &A.n, &k, &sa, 
			&A.M[0][0], &lda, &B.M[0][0], &lda, &sa, &C.M[0][0], &lda);

	PrintMatrix(A.getMatrix(), A.getRows(), A.getColumns(), (char*)"A: ");
	PrintMatrix(B.getMatrix(), B.getRows(), B.getColumns(), (char*)"B: ");
	PrintMatrix(C.getMatrix(), C.getRows(), C.getColumns(), (char*)"C: ");
}

Matrix::Matrix(int row, int col){
	M = NewMatrix(row, col);
	m = row;
	n = col;
	Zeros(&M, m, n);
}

Matrix::Matrix(char *handle){
	Read_MatrixMTX(&M, &m, &n, handle);
}

float **Matrix::NewMatrix(int dinamic_row, int dinamic_column){
	float **M = new float *[dinamic_row];
	for(int i=0; i<dinamic_row; i++){
		M[i] = new float [dinamic_column];
	}
	return M;
}

void Matrix::Zeros(float ***M, int m, int n){
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			(*M)[i][j] = 0;
		}
	}
}

void Matrix::Read_MatrixMTX(float ***M, int *m, int *n, char* handle){
	char str1[150];

	FILE * pFile;
	//int m=0, n=0;
	int nzv=0,fil,col;
	float data;

	pFile = fopen (handle, "rt");
	fgets(str1,150,pFile);

	fscanf(pFile, "%d %d %d\n", &*m,&*n,&nzv);

	(*M) = NewMatrix(*m, *n);
	for(int i=0;i<nzv;i++){
		fscanf(pFile, "%d %d %f\n", &fil,&col,&data);
		(*M)[fil-1][col-1] = data;
	}
	//PrintMatrix(*M, *m, *n, (char*)"M: ");
	fclose(pFile);
}

void PrintMatrix(float **M, int m, int n, char* msg){
	cout << msg << "\n";
	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			cout << setprecision(5) << M[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}
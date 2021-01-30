// g++ -std=c++11 jacobi.cpp -o jacobi -lf77blas

#include <algorithm> 
#include <iostream> 
#include <vector> 
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
using namespace std; 

typedef std::vector<float> Vector;
typedef vector<vector<float>> Matrix;
int const SIZE=3;

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

void ZerosFill_Matrix(Matrix&, int, int);
void PrintVector(const Vector&, char*);
void PrintMatrix(const Matrix&, char*);
void DLU(const Matrix&, Matrix&, Matrix&, Matrix&);

void getCofactor(Matrix&, Matrix&, int, int, int);
int determinant(Matrix& , int );
void adjoint(Matrix&, Matrix&);
bool inverse(Matrix&, Matrix&);

int main(){

	Matrix M = { {3,-1,-1},
		     {-1,3,1},
		     {2,1,4}};

	float A[3][3] =  { {3,-1,-1},
		     {-1,3,1},
		     {2,1,4}};

	float B[3] = {1,3,7};
	float R1[3] = {0,0,0};

	Matrix D, L , U, inv;

	Vector b = {1.0,3.0,7.0};
	int n=3,m=3;

	// auxiliar vectors
	float sa = 1, sb=1;
	int lda=3, incx=1, incy=1;
	Vector r1 = {0.0,0.0,0.0};


	DLU(M, D, L, U);

	if (inverse(D, inv)){
		PrintMatrix(inv, (char*)"inverse = ");
	} 


 	sgemv_((char*) "N", &m, &n, &sa, &A[0][0], &lda, &B[0], &incx, &sb, &R1[0], &incy);


	for(int i=0; i<n; i++){
		cout << R1[i] << "  ";
	}
	/*PrintMatrix(M, (char*)"M = ");
	PrintMatrix(D, (char*)"D = ");
	PrintMatrix(L, (char*)"L = ");
	PrintMatrix(U, (char*)"U = ");*/

	return 0;
}

void ZerosFill_Matrix(Matrix& M, int m, int n){
	for(int i=0; i < m; i++){
		Vector V;
		for(int j=0; j<n; j++){
			V.push_back(0.0);
		}
		M.push_back(V);
	}
}

void PrintVector(const Vector& V, char* msg){ 
	cout << msg << "[ "; 
	for_each(V.begin(), V.end(), [](int a) { 
		cout << setprecision(2) << a << " "; 
	}); 
	cout << "]" << endl; 
} 

void PrintMatrix(const Matrix& M, char* msg) { 
	int m = M.size();
	int n = M[0].size();
	cout << msg << "\n";
	for (int i = 0; i < m; i++) { 
		for (int j = 0; j < n; j++){
			cout << setprecision(2) << M[i][j] << " ";
		}        
		cout << endl; 
	}
	cout << endl;
} 

void DLU(const Matrix& M, Matrix& D, Matrix& L, Matrix& U){
	int m = M.size();
	int n = M[0].size();

	ZerosFill_Matrix(D, m, n);
	ZerosFill_Matrix(L, m, n);
	ZerosFill_Matrix(U, m, n);

	for(int i=0; i<m; i++){
		for(int j=0; j<n; j++){
			if(i==j){
				D[i][j] = M[i][j];
			} else if(j>i) {
				U[i][j] = -1*M[i][j];
			} else if(i>j) {
				L[i][j] = -1*M[i][j];
			}
		}
	}
}

// Function to get cofactor of A[p][q] in temp[][]. n is current 
// dimension of A[][] 
void getCofactor(Matrix& A, Matrix& temp, int p, int q, int n) 
{ 
    int i = 0, j = 0; 
  
    // Looping for each element of the matrix 
    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
            //  Copying into temporary matrix only those element 
            //  which are not in given row and column 
            if (row != p && col != q) 
            { 
                temp[i][j++] = A[row][col]; 
  
                // Row is filled, so increase row index and 
                // reset col index 
                if (j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 
} 
  
/* Recursive function for finding determinant of matrix. 
   n is current dimension of A[][]. */
int determinant(Matrix& A, int n) 
{ 
    int D = 0; // Initialize result 
  
    //  Base case : if matrix contains single element 
    if (n == 1) 
        return A[0][0]; 
  
    Matrix temp; // To store cofactors
    ZerosFill_Matrix(temp, n, n);
  
    int sign = 1;  // To store sign multiplier 
  
     // Iterate for each element of first row 
    for (int f = 0; f < n; f++) 
    { 
        // Getting Cofactor of A[0][f] 
        getCofactor(A, temp, 0, f, n); 
        D += sign * A[0][f] * determinant(temp, n - 1); 
  
        // terms are to be added with alternate sign 
        sign = -sign; 
    } 
  
    return D; 
} 
  
// Function to get adjoint of A[N][N] in adj[N][N]. 
void adjoint(Matrix& A, Matrix& adj) 
{ 
    int N = A.size();
    if (N == 1) 
    { 
        adj[0][0] = 1; 
        return; 
    } 
  
    // temp is used to store cofactors of A[][] 
    int sign = 1;
    Matrix temp; 
    ZerosFill_Matrix(temp, N, N);
  
    for (int i=0; i<N; i++) 
    { 
        for (int j=0; j<N; j++) 
        { 
            // Get cofactor of A[i][j] 
            getCofactor(A, temp, i, j, N); 
  
            // sign of adj[j][i] positive if sum of row 
            // and column indexes is even. 
            sign = ((i+j)%2==0)? 1: -1; 
  
            // Interchanging rows and columns to get the 
            // transpose of the cofactor matrix 
            adj[j][i] = (sign)*(determinant(temp, N-1)); 
        } 
    } 
} 
  
// Function to calculate and store inverse, returns false if 
// matrix is singular 
bool inverse(Matrix& A, Matrix& inverse) 
{ 
    // Find determinant of A[][] 
    int N = A.size();
    int det = determinant(A, N); 
    if (det == 0) 
    { 
        cout << "Singular matrix, can't find its inverse"; 
        return false; 
    } 
  
    // Find adjoint 
    Matrix adj;
    ZerosFill_Matrix(adj, N, N);
    ZerosFill_Matrix(inverse, N, N);
    adjoint(A, adj); 
  
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
    for (int i=0; i<N; i++) 
        for (int j=0; j<N; j++) 
            inverse[i][j] = adj[i][j]/float(det); 
  
    return true; 
} 

# Online Python compiler (interpreter) to run Python online.
# Write Python 3 code in this online editor and run it.
import numpy as np

A = np.array([[3,-1,-1], [-1,3,1], [2,1,4]])
b = np.array([1,3,7])
D = np.array( [[3,0,0], [0,3,0], [0,0,4]])
L = np.array( [[0,0,0], [1,0,0], [-2,-1,0]])
U = np.array( [[0,1,1],[0,0,-1],[0,0,0]])
xnew = np.array([0,0,0])
#xnuevo=(inv(D)*b+inv(D)*(L+U)*xnuevo)


d = np.array([3,3,3])
l = np.array([1,2,3])

#print(d*l)

C = np.zeros(3)
C = np.dot(1*A.T,b) + 1*C
print(C)


"""for i in range(1000):
    xnew = np.dot(np.linalg.inv(D),b) + np.dot(np.dot(np.linalg.inv(D), (L+U)), xnew ) 

print(xnew)"""

// g++ -std=c++11 csr.cpp -o csr

#include <algorithm> 
#include <iostream> 
#include <vector> 
#include <stdio.h>
#include <stdlib.h>
using namespace std; 
  
typedef std::vector<int> Vector;
typedef vector<vector<int>> Matrix; 

class CSR{
	private:
	Matrix M;
	Vector data;	
	Vector indices;
	Vector ptr = {0};

	public:
	// Constructor
	CSR(char*);
	//Getters
	Matrix getMatrix(){ return M; };
	Vector getData(){ return data; };
	Vector getCOL_INDEX(){ return indices; };
	Vector getROW_INDEX(){ return ptr; };
	// Methods Class CSR
	void ZerosFill_Matrix(Matrix&, int, int );
	void Read_MatrixMTX(Matrix&, char*);
	void PrintVector(const Vector&, char* );
	void PrintMatrix(const Matrix&, char* );
	void Print_SparseMatrix(char* );
	void Sparesify(const Matrix&, Vector&, Vector&, Vector&);
};
// Individual methods
void ZerosFill_Vector(Vector&, int);
void spmv_csr_serial(int, const Vector&, const Vector&, const Vector&, const Vector&, Vector&);

// Code 
int main(){ 
	// Read matrix
	CSR A((char*)"./matrixA.mtx");
	CSR b((char*)"./matrixb.mtx");
	Vector y;

	A.Print_SparseMatrix((char*)"A: ");
	b.Print_SparseMatrix((char*)"b: ");
	if(A.getMatrix()[0].size()==b.getMatrix().size()){
		spmv_csr_serial(A.getMatrix().size(), A.getROW_INDEX(), A.getCOL_INDEX(), A.getData(), b.getData(), y);
		A.PrintVector(y, (char*)"y=A*b\ny = ");
	}else{
		cout << "It is not possible the dot product!" << endl;
	}
	return 0; 
}

CSR::CSR (char* handle){
	Read_MatrixMTX(M, handle);
	Sparesify(M, data, indices, ptr);
}

void CSR::ZerosFill_Matrix(Matrix& M, int m, int n){
	for(int i=0; i < m; i++){
		Vector V;
		for(int j=0; j<n; j++){
			V.push_back(0);
		}
		M.push_back(V);
	}
}

void CSR::Read_MatrixMTX(Matrix& M, char* handle){
	char str1[150];

	FILE * pFile;
	int m=0,n=0,nzv=0,fil,col,data;

	pFile = fopen (handle,"rt");
	fgets(str1,150,pFile);

	fscanf(pFile, "%d %d %d\n", &m,&n,&nzv);

	ZerosFill_Matrix(M, m, n);
	for(int i=0;i<nzv;i++){
		fscanf(pFile, "%d %d %d\n", &fil,&col,&data);

		M[fil-1][col-1] = data;
	}
	fclose(pFile);
}



// Utility Function to print V, COL_INDEX, ROW_INDEX 
// with some decoration. 
void CSR::PrintVector(const Vector& V, char* msg) { 
	cout << msg << "[ "; 
	for_each(V.begin(), V.end(), [](int a) { 
		cout << a << " "; 
	}); 
	cout << "]" << endl; 
} 

// Utility Function to print a Matrix 
void CSR::PrintMatrix(const Matrix& M, char* msg) { 
	int m = M.size();
	int n = M[0].size();
	cout << msg << "\n";
	for (int i = 0; i < m; i++) { 
		for (int j = 0; j < n; j++){
			cout << M[i][j] << " ";
		}        
		cout << endl; 
	}
	cout << endl;
} 

void CSR::Print_SparseMatrix(char* msg){
	PrintMatrix(M, msg);
	PrintVector(data, (char*)"data = "); 
	PrintVector(indices, (char*)"indices = "); 
	PrintVector(ptr, (char*)"ptr = ");
	cout << endl;
}
// Generate the three vectors V, COL_INDEX, ROW_INDEX 
void CSR::Sparesify(const Matrix& M, Vector& data, Vector& indices, Vector& ptr) { 
	int m = M.size(); 
	int n = M[0].size();
	int NNZ = 0; 

	for (int i = 0; i < m; i++) { 
		for (int j = 0; j < n; j++) { 
			if (M[i][j] != 0) { 
				data.push_back(M[i][j]); 
				indices.push_back(j);
				// Count Number of Non Zero  
                		// Elements in row i 
				NNZ++;
			} 
		} 
		ptr.push_back(NNZ); 
	} 
} 

void ZerosFill_Vector(Vector& v, int n){
	for(int i=0; i<n; i++){
		v.push_back(0);	
	}
}
void spmv_csr_serial(int num_rows, const Vector& ptr, const Vector& indices, const Vector& data, const Vector& x, Vector& y){
	ZerosFill_Vector(y, num_rows);
	for (int i = 0; i < num_rows; i++){
		float dot = 0;

		int row_start = ptr[i]; 
		int row_end = ptr[i + 1];
		for (int j = row_start; j < row_end; j++){
			dot += data[j] * x[indices[j]];
		}
		y[i] += dot;
	}
}

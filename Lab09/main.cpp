#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include </usr/include/gsl/gsl_linalg.h>

typedef std::vector<double> Vector;
typedef std::vector< std::vector<double> > Matrix;


int l_index(const int& i, const int& j);
void ij_index(const int& l, int& i, int& j);
void save_vec_to_file(const Vector& vec, std::string name);
void save_rozklad(const Matrix& rozk, std::string name);

void CN(Matrix& A, Matrix& B, Vector& C, Vector& T);	
void ZnajdzRozklad(Vector& T, Matrix& Rozklad); 
void NadajWarunkiBrzegowe(Matrix& A, Matrix& B, Vector& C);
void NadajWarunkiPoczatkowe(Vector& T);

struct Zad {
	static double delta;
	static double d_t;
	static double T_a;
	static double T_b;
	static double T_c;
	static double T_d;
	static double k_B;
	static double k_D;
	static int IT_MAX;
	static int N;
	static int n_x;
	static int n_y;
};
double Zad::delta = 1;
double Zad::d_t = 1;
double Zad::T_a = 40;
double Zad::T_b = 0;
double Zad::T_c = 30;
double Zad::T_d = 0;
double Zad::k_B = 0.1;
double Zad::k_D = 0.6;
int Zad::IT_MAX = 2000;
int Zad::n_x = 40;
int Zad::n_y = 40;
int Zad::N = (Zad::n_x + 1) * (Zad::n_y + 1);



int main() {
	
	Vector row;
	row.resize(Zad::N, 0);
	Matrix Clear_M;
	Clear_M.resize(Zad::N, row);

	Matrix A;
	Matrix B;
	Matrix Rozklad;
	Vector C;
	Vector T;
 
	A = Clear_M;
	B = Clear_M;
	Rozklad = Clear_M;
	C = row;
	T = row;
	Zad::IT_MAX = 100;
	NadajWarunkiBrzegowe(A, B, C);
	NadajWarunkiPoczatkowe(T);
	CN(A, B, C, T);
	ZnajdzRozklad(T, Rozklad);
	save_vec_to_file(T, "zad1a.dat");
	save_rozklad(Rozklad, "zad2a.dat");
  std::cout << "Done1\n"; 


	A = Clear_M;
	B = Clear_M;
	C = row;
	T = row;
	Zad::IT_MAX = 200;
	NadajWarunkiBrzegowe(A, B, C);
	NadajWarunkiPoczatkowe(T);
	CN(A, B, C, T);
	ZnajdzRozklad(T, Rozklad);
	save_vec_to_file(T, "zad1b.dat");
	save_rozklad(Rozklad, "zad2b.dat");
  std::cout << "Done1\n";


	A = Clear_M;
	B = Clear_M;
	C = row;
	T = row;
	Zad::IT_MAX = 500;
	NadajWarunkiBrzegowe(A, B, C);
	NadajWarunkiPoczatkowe(T);
	CN(A, B, C, T);
	ZnajdzRozklad(T, Rozklad);
	save_vec_to_file(T, "zad1c.dat");
	save_rozklad(Rozklad, "zad2c.dat");
  std::cout << "Done1\n";


	A = Clear_M;
	B = Clear_M;
	C = row;
	T = row;
	Zad::IT_MAX = 1000;
	NadajWarunkiBrzegowe(A, B, C);
	NadajWarunkiPoczatkowe(T);
	CN(A, B, C, T);
	ZnajdzRozklad(T, Rozklad);
	save_vec_to_file(T, "zad1d.dat");
	save_rozklad(Rozklad, "zad2d.dat");
  std::cout << "Done1\n";


	A = Clear_M;
	B = Clear_M;
	C = row;
	T = row;
	Zad::IT_MAX = 2000;
	NadajWarunkiBrzegowe(A, B, C);
	NadajWarunkiPoczatkowe(T);
	CN(A, B, C, T);
	ZnajdzRozklad(T, Rozklad);
	save_vec_to_file(T, "zad1e.dat");
	save_rozklad(Rozklad, "zad2e.dat");
  std::cout << "Done1\n";

return 0;
} //end main



void CN(Matrix& A, Matrix& B, Vector& C, Vector& T) {

	gsl_matrix *A_gsl = gsl_matrix_calloc(Zad::N, Zad::N);
	gsl_matrix *B_gsl = gsl_matrix_calloc(Zad::N, Zad::N);
	gsl_vector *C_gsl = gsl_vector_calloc(Zad::N);
	gsl_vector *T_gsl = gsl_vector_calloc(Zad::N);
	gsl_vector *d = gsl_vector_calloc(Zad::N);
	gsl_permutation *p = gsl_permutation_alloc(Zad::N);

	for (int i = 0; i < Zad::N; ++i) {
		for (int j = 0; j < Zad::N; ++j) {
			gsl_matrix_set(A_gsl, i, j, A[i][j]);
			gsl_matrix_set(B_gsl, i, j, B[i][j]);
		}
		gsl_vector_set(C_gsl, i, C[i]);
		gsl_vector_set(T_gsl, i, T[i]);
	}

	double alfa = 1;
	double beta = 0;
	int signum;

	gsl_linalg_LU_decomp(A_gsl, p, &signum);
	for (int k = 0; k < Zad::IT_MAX; ++k) {
	    gsl_blas_dgemv(CblasNoTrans, alfa, B_gsl, T_gsl, beta, d);
	    gsl_blas_daxpy(alfa, C_gsl, d);
	    gsl_linalg_LU_solve(A_gsl, p, d, T_gsl);
	}

	for (int i = 0; i < Zad::N; ++i) {
		for (int j = 0; j < Zad::N; ++j) {
			A[i][j] = gsl_matrix_get(A_gsl, i, j);
			B[i][j] = gsl_matrix_get(B_gsl, i, j);
		}
		C[i] = gsl_vector_get(C_gsl, i);
		T[i] = gsl_vector_get(T_gsl, i);
	}
}



void ZnajdzRozklad(Vector& T, Matrix& Rozklad) {

	int l;
	int l_prawy;
	int l_dolny;
	int l_lewy;
	int l_gorny;

	double poch2_x;
	double poch2_y;
	
	for (int i = 1; i < Zad::n_x; ++i)
		for (int j = 1; j < Zad::n_y; ++j) {

			l = l_index(i, j);
			l_prawy = l_index(i+1, j);
			l_lewy = l_index(i-1, j);
			l_dolny = l_index(i, j-1);
			l_gorny = l_index(i, j+1);

			poch2_x = T[l_prawy] - 2*T[l] + T[l_lewy];

			poch2_y = T[l_gorny] - 2*T[l] + T[l_dolny];

			Rozklad[i][j] = poch2_y + poch2_x;
		}
}



void NadajWarunkiBrzegowe(Matrix& A, Matrix& B, Vector& C) {

	int l;
	double temp;

	temp = Zad::d_t / ( 2 * std::pow(Zad::delta, 2) );

	for(int i = 1; i < Zad::n_x; ++i) //Wnetrze 
		for(int j = 1; j < Zad::n_y; ++j) {
			l = l_index(i, j);
			
			A[l][l - Zad::n_x - 1] = temp;
			A[l][l-1] = temp;
			A[l][l+1] = temp;
			A[l][l + Zad::n_x + 1] = temp;
			A[l][l] = (-4) * temp - 1;

			B[l][l - Zad::n_x - 1] = -temp;
			B[l][l-1] = -temp;
			B[l][l+1] = -temp;
			B[l][l + Zad::n_x + 1] = -temp;
			B[l][l] = 4 * temp - 1;
		}

	for(int j = 0; j <= Zad::n_y; ++j) { //lewy prawy B
		l = l_index(0, j);
		A[l][l] = 1;
		B[l][l] = 1;
		C[l] = 0;
		l = l_index(Zad::n_x, j);
		A[l][l] = 1;
		B[l][l] = 1;
		C[l] = 0;
	}

	temp = 1 / ( Zad::k_B * Zad::delta );
	for(int i = 1; i < Zad::n_x; ++i) {
		l = l_index(i, Zad::n_y);
		A[l][l - Zad::n_x - 1] = -temp;
		A[l][l] = temp  +1;
		C[l] = Zad::T_b;

		for(int k = 0; k < Zad::N; ++k){
			B[l][k] = 0;
		}
	}

	temp = 1 / ( Zad::k_D * Zad::delta );
	for(int i = 1; i < Zad::n_x; ++i) {
		l = l_index(i, 0);
		A[l][l + Zad::n_x + 1] = -temp;
		A[l][l] = temp + 1;
		C[l] = Zad::T_d;

		for(int k = 0; k < Zad::N; ++k){
			B[l][k] = 0;
		}
	}
}



void NadajWarunkiPoczatkowe(Vector& T) {

	int l;

	for(int i = 0; i < Zad::n_x; ++i)
		for (int j = 0; j < Zad::n_y; ++j) {
			l = l_index(i, j);
			T[l] = 0;
		}

	for (int j = 0; j < Zad::n_y; ++j) {
		l = l_index(0, j);
		T[l] = Zad::T_a;
	}

	for(int j = 0; j < Zad::n_y; ++j) {
		l = l_index(Zad::n_x, j);
		T[l] = Zad::T_c;
	}
}



void ij_index(const int& l, int& i, int& j) {
	j = (int) l / (Zad::n_x + 1);
	i = l - j * (Zad::n_x + 1);
}



int l_index(const int& i, const int& j) {
	return j * (Zad::n_x + 1) + i;
}



void save_rozklad(const Matrix& rozk, std::string name) {

	FILE* file;

	file = fopen(name.c_str(), "w");
		for( int i = 1; i < Zad::n_x; ++i) {
			for(int j = 1; j < Zad::n_y; ++j) {
   			fprintf(file, "%10d %10d %15.9lf\n", i, j , rozk[i][j]);
			}
			fprintf(file, "\n");
		}

	fclose(file);
}



void save_vec_to_file(const Vector& vec, std::string name) {

	int i;
	int j;
	double x;
	double y;
	FILE* file;

	file = fopen(name.c_str(), "w");
		for (int l = 0; l < vec.size(); ++l) {
			if ( l % (Zad::n_x + 1) == 0 ) fprintf(file, "\n");
			ij_index(l, i, j);
			x = i * Zad::delta;
			y = j * Zad::delta;
			fprintf(file, "%15.9lf %15.9lf %15.9lf\n", x, y , vec[l]);
		}	
	fclose(file);
}
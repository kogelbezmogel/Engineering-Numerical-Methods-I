#include <math.h>
#include <string>
#include "mgmres.h"
#include <iostream>
#include <fstream>

//void save_vec_to_file(const Vector& vec, std::string name);
struct Zad {
	static int N_x;
	static int N_y;
	static int N;
	static double delta;
	static double epsilon_1;
	static double epsilon_2;
	static double V1;
	static double V2;
	static double V3;
	static double V4;
};

int Zad::N_x = 4;
int Zad::N_y = 4;
int Zad::N = (Zad::N_x+1) * (Zad::N_y+1);
double Zad::delta = 0.1;
double Zad::epsilon_1 = 1;
double Zad::epsilon_2 = 1;
double Zad::V1 = 10;
double Zad::V2 = -10;
double Zad::V3 = 10;
double Zad::V4 = -10;

//deklaracje funkcji

void wypelnij_wektory( double* a, int* ia, int* ja, double* b, double* V );

void znajdz_i_j(int l, int* i, int* j);

void zapisz_V(std::string name, double* V);

double A1(int l);
double A2(int l);
double A3(int l);
double A4(int l);
double A5(int l);

int main() {


	double* a = (double*) new double[5*Zad::N] {0};
	int* ja = (int*) new int[5*Zad::N] {0};
	int* ia = (int*) new int[Zad::N+1] {-1};
	double* b = (double*) new double[Zad::N] {0};
	double* V = (double*) new double[Zad::N] {0};

	int itr_max = 500;
	int mr = 500;
	double tol_abs = 1.e-8;
	double tol_rel = 1.e-8;

/* 									aby to wygenerowac trzeba odkomentowac odpowiednie funkcje Q_1 Q_2
	//1)
	wypelnij_wektory(a, ia, ja, b, V);
	int i;
	int j;

	FILE* file1 = fopen("mat_l_i_j_a.dat", "w");
		for(int l = 0; l < 5 * Zad::N; ++l) {
			znajdz_i_j(l, &i, &j);
			fprintf(file1, "%10d %10d %10d %15.4lf \n", l, i, j, a[l]);
		}
	fclose(file1);

	FILE* file2 = fopen("vec_l_i_j_b.dat", "w");
		int l;
		for(int l = 0; l < Zad::N; ++l) {
			znajdz_i_j(l, &i, &j);
			fprintf(file2, "%10d %10d %10d %15.4lf \n", l, i, j, b[l]);
		}
	fclose(file2);

	delete [] a;
	delete [] ja;
	delete [] ia;
	delete [] b;
	delete [] V;


	//2a) 
	Zad::N_x = 50;
	Zad::N_y = 50;
	Zad::N = (Zad::N_x+1) * (Zad::N_y+1);
	a = (double*) new double[5*Zad::N] {0};
	ja = (int*) new int[5*Zad::N] {0};
	ia = (int*) new int[Zad::N+1] {-1};
	b = (double*) new double[Zad::N] {0};
	V = (double*) new double[Zad::N] {0};
	wypelnij_wektory(a, ia, ja, b, V);
	pmgmres_ilu_cr( Zad::N, ia[Zad::N], ia, ja, a, V, b, itr_max, mr, tol_abs , tol_rel );
	zapisz_V("z1a.dat", V);
	delete [] a;
	delete [] ja;
	delete [] ia;
	delete [] b;
	delete [] V;

	//2b) 
	Zad::N_x = 100;
	Zad::N_y = 100;
	Zad::N = (Zad::N_x+1) * (Zad::N_y+1);
	a = (double*) new double[5*Zad::N] {0};
	ja = (int*) new int[5*Zad::N] {0};
	ia = (int*) new int[Zad::N+1] {-1};
	b = (double*) new double[Zad::N] {0};
	V = (double*) new double[Zad::N] {0};
	wypelnij_wektory(a, ia, ja, b, V);
	pmgmres_ilu_cr( Zad::N, ia[Zad::N], ia, ja, a, V, b, itr_max, mr, tol_abs , tol_rel );
	zapisz_V("z1b.dat", V);
	delete [] a;
	delete [] ja;
	delete [] ia;
	delete [] b;
	delete [] V;

	//2c) 
	Zad::N_x = 200;
	Zad::N_y = 200;
	Zad::N = (Zad::N_x+1) * (Zad::N_y+1);
	a = (double*) new double[5*Zad::N] {0};
	ja = (int*) new int[5*Zad::N] {0};
	ia = (int*) new int[Zad::N+1] {-1};
	b = (double*) new double[Zad::N] {0};
	V = (double*) new double[Zad::N] {0};
	wypelnij_wektory(a, ia, ja, b, V);
	pmgmres_ilu_cr( Zad::N, ia[Zad::N], ia, ja, a, V, b, itr_max, mr, tol_abs , tol_rel );
	zapisz_V("z1c.dat", V);
	delete [] a;
	delete [] ja;
	delete [] ia;
	delete [] b;
	delete [] V;
								aby to wygenerowac trzeba odkomentowac odpowiedznie funkcje Q_1 Q_2*/ 
	//3a) 
	Zad::N_x = 100;
	Zad::N_y = 100;
	Zad::N = (Zad::N_x+1) * (Zad::N_y+1);
	Zad::epsilon_1 = 1;
	Zad::epsilon_2 = 1;
	Zad::V1 = 0;
	Zad::V2 = 0;
	Zad::V3 = 0;
	Zad::V4 = 0;
	a = (double*) new double[5*Zad::N] {0};
	ja = (int*) new int[5*Zad::N] {0};
	ia = (int*) new int[Zad::N+1] {-1};
	b = (double*) new double[Zad::N] {0};
	V = (double*) new double[Zad::N] {0};

	wypelnij_wektory(a, ia, ja, b, V);
	pmgmres_ilu_cr( Zad::N, ia[Zad::N], ia, ja, a, V, b, itr_max, mr, tol_abs , tol_rel );
	zapisz_V("z2a.dat", V);
	delete [] a;
	delete [] ja;
	delete [] ia;
	delete [] b;
	delete [] V;

	//3b) 
	Zad::N_x = 100;
	Zad::N_y = 100;
	Zad::N = (Zad::N_x+1) * (Zad::N_y+1);
	Zad::epsilon_1 = 1;
	Zad::epsilon_2 = 2;
	Zad::V1 = 0;
	Zad::V2 = 0;
	Zad::V3 = 0;
	Zad::V4 = 0;
	a = (double*) new double[5*Zad::N] {0};
	ja = (int*) new int[5*Zad::N] {0};
	ia = (int*) new int[Zad::N+1] {-1};
	b = (double*) new double[Zad::N] {0};
	V = (double*) new double[Zad::N] {0};

	wypelnij_wektory(a, ia, ja, b, V);
	pmgmres_ilu_cr( Zad::N, ia[Zad::N], ia, ja, a, V, b, itr_max, mr, tol_abs , tol_rel );
	zapisz_V("z2b.dat", V);
	delete [] a;
	delete [] ja;
	delete [] ia;
	delete [] b;
	delete [] V;

	//3c) 
	Zad::N_x = 100;
	Zad::N_y = 100;
	Zad::N = (Zad::N_x+1) * (Zad::N_y+1);
	Zad::epsilon_1 = 1;
	Zad::epsilon_2 = 10;
	Zad::V1 = 0;
	Zad::V2 = 0;
	Zad::V3 = 0;
	Zad::V4 = 0;
	a = (double*) new double[5*Zad::N] {0};
	ja = (int*) new int[5*Zad::N] {0};
	ia = (int*) new int[Zad::N+1] {-1};
	b = (double*) new double[Zad::N] {0};
	V = (double*) new double[Zad::N] {0};

	wypelnij_wektory(a, ia, ja, b, V);
	pmgmres_ilu_cr( Zad::N, ia[Zad::N], ia, ja, a, V, b, itr_max, mr, tol_abs , tol_rel );
	zapisz_V("z2c.dat", V);
	delete [] a;
	delete [] ja;
	delete [] ia;
	delete [] b;
	delete [] V;


	return 0;
}

///////////////////////////////////////////definicje

/*double Q_1(int i, int j) {
	return 0;
}

double Q_2(int i, int j) {
	return 0;
}*/

double Q_1(int i, int j) {

	double x_max = Zad::delta * Zad::N_x;
	double y_max = Zad::delta * Zad::N_y;
	double x = i * Zad::delta;  
	double y = j * Zad::delta;
	double sigma = x_max / 10.0;

	double temp = pow( (x - 0.25 * x_max) / sigma, 2);
	temp += pow( (y - 0.5 * y_max) / sigma, 2);
	temp = exp(- temp );

	return temp;
}

double Q_2(int i, int j) {

	double x_max = Zad::delta * Zad::N_x;
	double y_max = Zad::delta * Zad::N_y;
	double x = i * Zad::delta;  
	double y = j * Zad::delta;
	double sigma = x_max / 10.0;

	double temp = pow( (x - 0.75 * x_max) / sigma, 2);
	temp += pow( (y - 0.5 * y_max) / sigma, 2);
	temp = -exp(- temp );

	return temp;
}

void znajdz_i_j(int l, int* i, int* j) {
	*j = (int)  (l / (Zad::N_x+1));
	*i =  (int) (l - (*j)*(Zad::N_x+1));
}

void wypelnij_wektory( double* a, int* ia, int* ja, double* b, double* V ) {

	int k = -1;
	int i = 0;
	int j = 0;

	for(int l = 0; l < Zad::N; ++l) {

		znajdz_i_j(l, &i, &j);

		int brzeg = 0;
		double vb = 0;
	
		if( i == 0 ) {
			brzeg = 1;
			vb = Zad::V1;
		}
		if( j == Zad::N_y ) {
			brzeg = 1;
			vb = Zad::V2;
		}
		if( i == Zad::N_x ) {
			brzeg = 1;
			vb = Zad::V3;
		}
		if( j == 0 ) {
			brzeg = 1;
			vb = Zad::V4;
		}

		b[l] = -(Q_1(i, j) + Q_2(i, j));

		if(brzeg == 1) {
			b[l] = vb;
		}
		ia[l] = - 1;

		if( l - Zad::N_x-1 >= 0 && brzeg == 0 ) {
			++k;
			if( ia[l] < 0 )
				ia[l] = k;
			a[k] = A1(l);
			ja[k] = l - Zad::N_x - 1;
		}

		if(l-1 >= 0 && brzeg == 0 ) {
			++k;
			if( ia[l] < 0 )
				ia[l] = k;
		a[k] = A2(l);
		ja[k] = l -1 ;
		}

		//diagonala
		++k;
		if( ia[l] < 0 )
			ia[l] = k;
		if(brzeg == 0)
			a[k] = A3(l);
		else
			a[k] = 1;
		ja[k] = l;

		if( l < Zad::N && brzeg == 0 ) {
			++k;
			a[k] = A4(l);
			ja[k] = l + 1;
		}

		if( l < Zad::N - Zad::N_x - 1 && brzeg == 0 ) {
			++k;
			a[k] = A5(l);
			ja[k] = l + Zad::N_x + 1;
		}

		ia[Zad::N] = k+1;
	}
}

double A1(int l) {
	int i;
	int j;
	znajdz_i_j(l, &i, &j);
	double e = i < Zad::N_x/2 ? Zad::epsilon_1 : Zad::epsilon_2;
	return e / pow(Zad::delta, 2);
}

double A2(int l) {
	int i;
	int j;
	znajdz_i_j(l, &i, &j);
	double e = i < Zad::N_x/2 ? Zad::epsilon_1 : Zad::epsilon_2;
	return e / pow(Zad::delta, 2);
}

double A3(int l) {
	int i_1;
	int i_2;
	int i_3;
	int j;
	znajdz_i_j(l, &i_1, &j);
	znajdz_i_j(l+1, &i_2, &j);
	znajdz_i_j(l + Zad::N_x + 1, &i_3, &j);
	double e_1 = i_1 < Zad::N_x/2 ? Zad::epsilon_1 : Zad::epsilon_2;
	double e_2 = i_2 < Zad::N_x/2 ? Zad::epsilon_1 : Zad::epsilon_2;
	double e_3 = i_3 < Zad::N_x/2 ? Zad::epsilon_1 : Zad::epsilon_2;
	return - (2*e_1 + e_2 +  e_3) / pow(Zad::delta, 2);
}

double A4(int l) {
	int i;
	int j;
	znajdz_i_j(l+1, &i, &j);
	double e = i < Zad::N_x/2 ? Zad::epsilon_1 : Zad::epsilon_2;
	return e / pow(Zad::delta, 2);
}

double A5(int l) {
	int i;
	int j;
	znajdz_i_j(l + Zad::N_x + 1, &i, &j);
	double e = i < Zad::N_x/2 ? Zad::epsilon_1 : Zad::epsilon_2;
	return e / pow(Zad::delta, 2);
}


void zapisz_V(std::string name, double* V) {

	FILE* file = fopen(name.c_str(), "w");
		int l;
		for(int i = 0; i <= Zad::N_x; ++i) {
			for( int j = 0; j <= Zad::N_y; ++j) {
				l = j * (Zad::N_x + 1) + i;
				fprintf(file, "%5.3lf %5.3lf %15.11lf \n", i*Zad::delta, j*Zad::delta, V[l]);
			}
			fprintf(file, "\n");
		}
	fclose(file);
}
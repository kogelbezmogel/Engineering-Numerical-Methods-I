#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

typedef std::vector< std::vector<double> > Matrix;
typedef std::vector<double> Vector;

struct Zad {
	static double d_t;
	static double delta;
	static double sigma;
	static double D;
	static double xa;
	static double ya;
	static int n_x;
	static int n_y;
	static int i_1;
	static int i_2;
	static int j_1;
};

double Zad::d_t;
double Zad::delta;
double Zad::sigma;
double Zad::D;
double Zad::xa;
double Zad::ya;
int Zad::n_x = 400;
int Zad::n_y = 90;
int Zad::i_1;
int Zad::i_2;
int Zad::j_1;

void MetodaCrankaNicolsona(Matrix& mat_u0, Matrix& mat_u1, Matrix& mat_vx, Matrix& mat_vy, Matrix& mat_psi, Vector& gestosc, Vector& polozenie, Vector& czasu );
double ZnajdzVmax( Matrix& mat_vx, Matrix& mat_vy);
double ZnajdzD_t(double v_max);
void WarunekPoczatkowyU(Matrix& mat_u0);
double WartoscCalkiGestosci(const Matrix& mat_u);
double WartoscPolozeniaSredniego(const Matrix& mat_u);
void WczytajPoleV(std::string name, Matrix& mat_vx, Matrix& mat_vy, Matrix& mat_psi);

void clear_mats( std::initializer_list< Matrix > list ) { }
void save_mat_to_file(const Matrix& mat, std::string name);
void save_vec_to_file(const Vector& vec1, const Vector& vec2, std::string name);

int main( )
{
	Vector row;
	Matrix clear;
	Matrix mat_u1;
	Matrix mat_u0;
	Matrix mat_vx;
	Matrix mat_vy;
	Matrix mat_psi;

	row.resize(Zad::n_y + 1, 0);
	clear.resize(Zad::n_x + 1, row);
	mat_vx = clear;
	mat_vy = clear;
	mat_u1 = clear;
	mat_u0 = clear;
	mat_psi = clear;

	Vector vec_gestosci;
	Vector vec_polozenia;
	Vector vec_czasu;

	Zad::delta = 0.01;
	Zad::sigma = 10 * Zad::delta;
	Zad::D = 0;
	Zad::xa = 0.45;
	Zad::ya = 0.45;
	Zad::n_x = 400;
	Zad::n_y = 90;
	Zad::i_1 = 200;
	Zad::i_2 = 210;
	Zad::j_1 = 50;

	WczytajPoleV("phi.dat", mat_vx, mat_vy, mat_psi);
	//save_mat_to_file( mat_vx, "zad1_vx.dat" );
	//save_mat_to_file( mat_vy, "zad1_vy.dat" );

	Zad::d_t = ZnajdzD_t( ZnajdzVmax(mat_vx, mat_vy) );

	MetodaCrankaNicolsona(mat_u0, mat_u1, mat_vx, mat_vy, mat_psi, vec_gestosci, vec_polozenia, vec_czasu );
	save_vec_to_file( vec_czasu, vec_polozenia, "zad1_x_t.dat");
	save_vec_to_file( vec_czasu, vec_gestosci, "zad1_g_t.dat");
	std::cout << "Done" << std::endl;
	return 0;
}

/////////////////////////////////////////////definicje


void MetodaCrankaNicolsona(Matrix& mat_u0, Matrix& mat_u1, Matrix& mat_vx, Matrix& mat_vy, Matrix& mat_psi, Vector& gestosc, Vector& polozenie, Vector& czasu ) {

	WarunekPoczatkowyU(mat_u0);

	mat_u1 = mat_u0;
	int ip;
	int il;
	int mat_number = 0;

	int ITE_MAX = 13000;

	for( int ite = 0; ite < ITE_MAX; ++ite) {
		
		for( int k = 0; k < 20; ++k) {

			for( int i = 0; i <= Zad::n_x; ++i ) {
				for( int j = 1; j < Zad::n_y; ++j ) {
					if( i >= Zad::i_1 && i <= Zad::i_2 && j <= Zad::j_1 ) {
						continue;
					}
					else {

						i == 0 ? il = Zad::n_x : il = i-1;
						i == Zad::n_x ? ip = 0 : ip = i+1;  

						mat_u1[i][j] = mat_u0[i][j];
						mat_u1[i][j] -= Zad::d_t / 2 * mat_vx[i][j] * ( (mat_u0[ip][j] - mat_u0[il][j]) / (2 * Zad::delta) + (mat_u1[ip][j] - mat_u1[il][j]) / (2 * Zad::delta) );
						mat_u1[i][j] -= Zad::d_t / 2 * mat_vy[i][j] * ( (mat_u0[i][j+1] - mat_u0[i][j-1]) / (2 * Zad::delta) + (mat_u1[i][j+1] - mat_u1[i][j-1]) / (2 * Zad::delta) );
						mat_u1[i][j] += Zad::d_t / 2 * Zad::D * ( mat_u0[ip][j] + mat_u0[il][j] + mat_u0[i][j+1] + mat_u0[i][j-1] - 4 * mat_u0[i][j]) / std::pow(Zad::delta, 2); 
						mat_u1[i][j] += Zad::d_t / 2 * Zad::D * ( mat_u1[ip][j] + mat_u1[il][j] + mat_u1[i][j+1] + mat_u1[i][j-1] ) / std::pow(Zad::delta, 2);
						mat_u1[i][j] *= 1. / (1 + 2 * Zad::D * Zad::d_t / std::pow(Zad::delta, 2));
					}//else
				}//j
			}//i
	
		} //k
		
		gestosc.push_back( WartoscCalkiGestosci(mat_u1) );
		polozenie.push_back( WartoscPolozeniaSredniego(mat_u1) );	
		czasu.push_back( ite * Zad::d_t );

		if(ite % 1000 == 0) {
			std::cout << "ite: " << ite << std::endl;
		}

		if( ite % (ITE_MAX / 5) == 0 ) {
			mat_number++;
			std::string line = "zad1u_";
			line += std::to_string(mat_number);
			line += ".dat";
			save_mat_to_file(mat_u1, line);
			std::cout << "saving " << mat_number << std::endl;
		}

		mat_u0 = mat_u1; 
	} //ite
}


double ZnajdzVmax( Matrix& mat_vx, Matrix& mat_vy) {

	double v_max = -1;
	double temp;

	for(int i = 1; i < Zad::n_x; ++i)
		for(int j = 1; j < Zad::n_y; ++j) {
			temp = std::sqrt( std::pow(mat_vx[i][j], 2) + std::pow(mat_vy[i][j], 2) );
			//std::cout << temp << "  1: " << std::pow(mat_vx[i][j], 2) << "  2:" << std::pow(mat_vy[i][j], 2) << std ::endl;
			if(temp > v_max)
				v_max = temp;
		}
	return v_max;
}

double ZnajdzD_t(double v_max) {
	return Zad::delta / (4 * v_max);
}

void WarunekPoczatkowyU(Matrix& mat_u0) {

	double x;
	double y;
	double maks = -1;

	for(int i = 0; i <= Zad::n_x; ++i) {
		for(int j = 0; j <= Zad::n_y; ++j) {
			x = i * Zad::delta;
			y = j * Zad::delta;
			mat_u0[i][j] = -std::pow(x-Zad::xa, 2) - std::pow(y-Zad::ya, 2);
			mat_u0[i][j] *= 1. / ( 2 * std::pow(Zad::sigma, 2) );
			mat_u0[i][j] = exp(mat_u0[i][j]);
			mat_u0[i][j] *= 1. / ( 2 * M_PI * std::pow(Zad::sigma, 2) );
			if(maks < mat_u0[i][j]) maks = mat_u0[i][j];
		}
	}
}

double WartoscCalkiGestosci(const Matrix& mat_u) {
	
	double c = 0;
	for(int i = 0; i <= Zad::n_x; ++i)
		for(int j = 0; j < Zad::n_y; ++j)
			c += mat_u[i][j] * std::pow(Zad::delta, 2);
	return c;
}

double WartoscPolozeniaSredniego(const Matrix& mat_u) {
	
	double x_sr = 0;
	double x_i;
	for(int i = 0; i <= Zad::n_x; ++i)
		for(int j = 0; j < Zad::n_y; ++j) {
			x_i = i * Zad::delta;
			x_sr += x_i * mat_u[i][j] * std::pow(Zad::delta, 2);
		}
	return x_sr;
}

void WczytajPoleV(std::string name, Matrix& mat_vx, Matrix& mat_vy, Matrix& mat_psi) {
	
	double value;
	int x;
	int y;
	std::ifstream input(name.c_str());

	while( input >> x >> y >> value ) {
		mat_psi[x][y] = value;
	}

	for(int i = 1; i < Zad::n_x; ++i)
		for( int j = 1; j < Zad::n_y; ++j) {
			mat_vx[i][j] = (mat_psi[i][j+1] - mat_psi[i][j-1]) / (2*Zad::delta);
			
		}
	
	for(int i = 1; i < Zad::n_x; ++i)
		for( int j = 1; j < Zad::n_y; ++j)
			mat_vy[i][j] = -(mat_psi[i+1][j] - mat_psi[i-1][j]) / (2*Zad::delta);

	//Warunki brzegowe
	for(int i = Zad::i_1; i <= Zad::i_2; ++i)
		for( int j = 0; j <= Zad::j_1; ++j) {
			mat_vx[i][j] = 0;
			mat_vy[i][j] = 0;
		}

	for(int i = 1; i < Zad::n_x; ++i){
		mat_vx[i][0] = 0;
		mat_vy[i][Zad::n_y] = 0;
	}

	for(int j = 0; j < Zad::n_y; ++j) {
		mat_vx[0][j] = mat_vx[1][j];
		mat_vx[Zad::n_x][j] = mat_vx[Zad::n_x-1][j];
	}
}

void save_vec_to_file(const Vector& vec1, const Vector& vec2, std::string name) {

	FILE* file;
	file = fopen(name.c_str(), "w");
		for (int i = 0; i < vec1.size(); ++i)
			fprintf(file, "%15.9lf %15.9lf\n", vec1[i], vec2[i]);
	fclose(file);
}

void save_mat_to_file(const Matrix& mat, std::string name) {

	FILE* file;
	file = fopen(name.c_str(), "w");
		for(int i = 0; i < mat.size(); ++i){
			for(int j = 0; j < mat[0].size(); ++j){
				fprintf(file, "%5.2lf %5.2lf %15.8lf\n", i*Zad::delta, j*Zad::delta, mat[i][j]);
			}
			fprintf(file, "\n");
		}
	fclose(file);
}

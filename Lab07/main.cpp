#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <cmath>

typedef std::vector< std::vector<double> > Matrix;
typedef std::vector<double> Vector;

struct Zad {
	static double mi;
	static double ro;
	static double delta;
	static int n_y;
	static int n_x;
	static int y_1;
	static int x_1;
	static int IT_MAX;
	static double Q_we;
	static double Q_wy;
};

using std::pow;

int Zad::n_x = 200;
int Zad::n_y = 90;
int Zad::x_1 = 50;
int Zad::y_1 = 55;
int Zad::IT_MAX = 20000;
double Zad::delta = 0.01;
double Zad::ro = 1;
double Zad::mi = 1;
double Zad::Q_we = -1000;
double Zad::Q_wy = 0;

void Q_WY() {
	double yn = Zad::n_y * Zad::delta;
	double y1 = Zad::y_1 * Zad::delta;
	Zad::Q_wy = Zad::Q_we * (pow(yn, 3) - pow(y1, 3) - 3*y1* pow(yn, 2) + 3*pow(y1, 2)*yn) / pow(yn, 3);
}


void RelaksacjaUkladu( Matrix& mat_psi, Matrix& mat_dzeta );
void DzetaBorder( Matrix& mat_psi, Matrix& mat_dzeta );
void PsiBorder( Matrix& mat_psi );
double ErrorControl( Matrix& mat_psi, Matrix& mat_dzeta );

void MapaPredkosicX( Matrix& mat_psi,  Matrix& matV_y);
void MapaPredkosicY( Matrix& mat_psi, Matrix& mapV_x);

double max_mat(const Matrix& mat);
double min_mat(const Matrix& mat);
void save_mat_to_file(const Matrix& mat, std::string name);
void save_vec_to_file(const Vector& vec, std::string name);

int main() {

	Matrix mat_dzeta;
	Vector row;
	row.resize(Zad::n_y+1, 0);
	mat_dzeta.resize(Zad::n_x+1, row);
	Matrix mat_psi = mat_dzeta;
	Matrix clear = mat_dzeta;
	Q_WY();

	Matrix Vx;
	Matrix Vx_clear;
	Vx.resize(Zad::n_x, row);
	row.resize(Zad::n_y, 0);
	Matrix Vy;
	Matrix Vy_clear;
	Vy.resize(Zad::n_x+1, row);

	Vx_clear = Vx;
	Vy_clear = Vy;

	//Zad1

	PsiBorder(mat_psi);
	DzetaBorder(mat_psi, mat_dzeta);
	RelaksacjaUkladu(mat_psi, mat_dzeta);
	MapaPredkosicX(mat_psi, Vx);
	MapaPredkosicY(mat_psi, Vy);
	std::cout << "1 psi)   Max: " << max_mat(mat_psi) << " Min: " << min_mat(mat_psi) << std::endl;
	std::cout << "1 dzeta) Max: " << max_mat(mat_dzeta) << " Min: " << min_mat(mat_dzeta) << std::endl << std::endl;
	save_mat_to_file(Vx, "vx1.dat");
	save_mat_to_file(Vy, "vy1.dat");
	save_mat_to_file(mat_psi, "psi1.dat");
	save_mat_to_file(mat_dzeta, "dzeta1.dat");


	//Zad2

	mat_psi = clear;
	mat_dzeta = clear;
	Zad::Q_we = -4000;
	Q_WY();

	PsiBorder(mat_psi);
	DzetaBorder(mat_psi, mat_dzeta);
	RelaksacjaUkladu(mat_psi, mat_dzeta);
	MapaPredkosicX(mat_psi, Vx);
	MapaPredkosicY(mat_psi, Vy);
	std::cout << "2 psi)   Max: " << max_mat(mat_psi) << " Min: " << min_mat(mat_psi) << std::endl;
	std::cout << "2 dzeta) Max: " << max_mat(mat_dzeta) << " Min: " << min_mat(mat_dzeta) << std::endl << std::endl;
	save_mat_to_file(Vx, "vx2.dat");
	save_mat_to_file(Vy, "vy2.dat");
	save_mat_to_file(mat_psi, "psi2.dat");
	save_mat_to_file(mat_dzeta, "dzeta2.dat");

	//Zad3

	mat_psi = clear;
	mat_dzeta = clear;
	Zad::Q_we = 4000;
	Q_WY();

	PsiBorder(mat_psi);
	DzetaBorder(mat_psi, mat_dzeta);
	RelaksacjaUkladu(mat_psi, mat_dzeta);
	MapaPredkosicX(mat_psi, Vx);
	MapaPredkosicY(mat_psi, Vy);
	std::cout << "3 psi)   Max: " << max_mat(mat_psi) << " Min: " << min_mat(mat_psi) << std::endl;
	std::cout << "3 dzeta) Max: " << max_mat(mat_dzeta) << " Min: " << min_mat(mat_dzeta) << std::endl;
	save_mat_to_file(Vx, "vx3.dat");
	save_mat_to_file(Vy, "vy3.dat");
	save_mat_to_file(mat_psi, "psi3.dat");
	save_mat_to_file(mat_dzeta, "dzeta3.dat");

	std::cout << std::endl << "Done" << std::endl;

	return 0;
}



///////////////////////////////////////Definicje



void RelaksacjaUkladu( Matrix& mat_psi, Matrix& mat_dzeta ) {

	int omega;
	double temp;
	int ite = 0;

	Matrix mat_psi_cp = mat_psi;
	Matrix mat_dzeta_cp = mat_dzeta;

	while( ite < Zad::IT_MAX ) {

		if( ite < 2000 ) omega = 0;
		else omega = 1;

		for( int i = 1; i < Zad::n_x; ++i ) {
			for( int j = 1; j < Zad::n_y; ++j ) {
				if( i > Zad::x_1 || j > Zad::y_1 ) {
					mat_psi_cp[i][j] = mat_psi[i+1][j] + mat_psi[i-1][j] + mat_psi[i][j+1] + mat_psi[i][j-1] - pow(Zad::delta,2) * mat_dzeta[i][j];
					mat_psi_cp[i][j] *= 0.25; 
				}
			}
		}
		mat_psi = mat_psi_cp;

		for( int i = 1; i < Zad::n_x; ++i ) {
			for( int j = 1; j < Zad::n_y; ++j ) {
				if( i > Zad::x_1 || j > Zad::y_1 ) {
					mat_dzeta_cp[i][j] = mat_dzeta[i+1][j] + mat_dzeta[i-1][j] + mat_dzeta[i][j+1] + mat_dzeta[i][j-1]; 
					mat_dzeta_cp[i][j] *= 0.25;

					temp = (mat_psi[i][j+1] - mat_psi[i][j-1]) * (mat_dzeta[i+1][j] - mat_dzeta[i-1][j]);
					temp -= (mat_psi[i+1][j] - mat_psi[i-1][j]) * (mat_dzeta[i][j+1] - mat_dzeta[i][j-1]);
					temp *= omega * Zad::ro / (16 * Zad::mi);

					mat_dzeta_cp[i][j] -= temp;
				}
			}
		}
		mat_dzeta = mat_dzeta_cp;

		DzetaBorder( mat_psi, mat_dzeta );
		//std::cout << "Err: " << ErrorControl(mat_psi, mat_dzeta) << std::endl;
		++ite;
	}
}

double ErrorControl( Matrix& mat_psi, Matrix& mat_dzeta ) {
	double err = 0;
	for( int j = 1; j < Zad::n_y; j+=2) {
		for(int i = 1; i < Zad::n_x; ++i) {
			err += mat_psi[i+1][j] + mat_psi[i-1][j] + mat_psi[i][j+1] + mat_psi[i][j-1] - 4*mat_psi[i][j] - pow(Zad::delta, 2)*mat_dzeta[i][j];
		}
	}
	return err;
}

void MapaPredkosicX( Matrix& mat_psi,  Matrix& matV_x) {
	for(int i = 0; i < Zad::n_x; ++i)
		for(int j = 0; j <= Zad::n_y; ++j)
			if(i > Zad::x_1 || j > Zad::y_1)
				matV_x[i][j] = (mat_psi[i+1][j] - mat_psi[i][j]) / Zad::delta;
}

void MapaPredkosicY( Matrix& mat_psi, Matrix& mapV_y) {
	for(int i = 0; i <= Zad::n_x; ++i)
		for(int j = 0; j < Zad::n_y; ++j)
			if(i > Zad::x_1 || j > Zad::y_1)
				mapV_y[i][j] = (mat_psi[i][j+1] - mat_psi[i][j]) / Zad::delta;
}


void DzetaBorder( Matrix& mat_psi, Matrix& mat_dzeta ) {

	double y;
	double y_1 = Zad::y_1 * Zad::delta;
	double y_n = Zad::n_y * Zad::delta;

	for(int j = Zad::y_1; j <= Zad::n_y; ++j) {//A
		y = Zad::delta*j;
		mat_dzeta[0][j] = Zad::Q_we/(2*Zad::mi)*(2*y - y_1 - y_n);
	}

	for(int j = 0; j <= Zad::n_y; ++j) { //C
		y = Zad::delta*j;
		mat_dzeta[Zad::n_x][j] = Zad::Q_wy / (2*Zad::mi) * (2*y - y_n);
	}

	for(int i = 1; i < Zad::n_x; ++i) { //B
		mat_dzeta[i][Zad::n_y] = 2 / pow(Zad::delta, 2) * ( mat_psi[i][Zad::n_y-1] - mat_psi[i][Zad::n_y] );
	}

	for(int i = Zad::x_1+1; i < Zad::n_x; ++i) { //D
		mat_dzeta[i][0] = 2 / pow(Zad::delta, 2) * ( mat_psi[i][1] - mat_psi[i][0] );
	}

	for(int j = 0; j < Zad::y_1; ++j) { //E
		mat_dzeta[Zad::x_1][j] = 2 / pow(Zad::delta, 2) * ( mat_psi[Zad::x_1+1][j] - mat_psi[Zad::x_1][j] );  
	}

	for(int i = 1; i <= Zad::x_1; ++i) { //F 
		mat_dzeta[i][Zad::y_1] = 2 / pow(Zad::delta, 2) * ( mat_psi[i][Zad::y_1+1] - mat_psi[i][Zad::y_1] );
	}

	mat_dzeta[Zad::x_1][Zad::y_1] = ( mat_dzeta[Zad::x_1-1][Zad::y_1] + mat_dzeta[Zad::x_1][Zad::y_1-1] ) / 2;
}

void PsiBorder( Matrix& mat_psi ) {

	double y;
	double y_1 = Zad::y_1 * Zad::delta;
	double y_n = Zad::n_y * Zad::delta;

	for(int j = Zad::y_1; j <= Zad::n_y; ++j) { //A
		y = Zad::delta*j;
		mat_psi[0][j] = Zad::Q_we / (2 * Zad::mi);
		mat_psi[0][j] *= (pow(y, 3)/3 - pow(y, 2)/2 * (y_1 + y_n) + y * y_1 * y_n);
	}

	for(int j = 0; j <= Zad::n_y; ++j) { //C
		y = Zad::delta*j;
		mat_psi[Zad::n_x][j] = Zad::Q_wy / (2 * Zad::mi);
		mat_psi[Zad::n_x][j] *= pow(y, 3)/3 - pow(y, 2)/2 * y_n;
		mat_psi[Zad::n_x][j] += Zad::Q_we * pow(y_1, 2) * (3*y_n - y_1) / (12 * Zad::mi);
	}

	for( int i = 1; i < Zad::n_x; ++i) { //B
		mat_psi[i][Zad::n_y] = mat_psi[0][Zad::n_y];
	}

	for( int i = Zad::x_1; i < Zad::n_x; ++i) { //D
		mat_psi[i][0] = mat_psi[0][Zad::y_1];
	}

	for( int j = 1; j <= Zad::y_1; ++j) { //E
		mat_psi[Zad::x_1][j] = mat_psi[0][Zad::y_1];
	}

	for( int i = 1; i <= Zad::x_1; ++i) { //F
		mat_psi[i][Zad::y_1] = mat_psi[0][Zad::y_1];
	}
}

void save_vec_to_file(const Vector& vec, std::string name) {

	FILE* file;
	file = fopen(name.c_str(), "w");
		for(int i = 0; i < vec.size(); ++i)
			fprintf(file, "%5d %15.9lf\n", i, vec[i]);
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

double max_mat(const Matrix& mat ) {
	double m = mat[100][40];
	for(int i = 0; i < mat.size(); ++i)
		for(int j = 0; j < mat[0].size(); ++j)
			if(m < mat[i][j] && mat[i][j] != 0)
				m = mat[i][j];
	return m;
}


double min_mat(const Matrix& mat ) {
	double m = mat[100][40];
	for(int i = 0; i < mat.size(); ++i)
		for(int j = 0; j < mat[0].size(); ++j)
			if(m > mat[i][j] && mat[i][j] != 0 )
				m = mat[i][j];
	return m;
}
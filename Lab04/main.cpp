#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>
#include <string>

typedef std::vector< std::vector<double> > Matrix;

struct Zad {
	static double epsilon;
	static double TOL;
	static double delta;
	static double n_x;
	static double n_y;
	static double V1;
	static double V2;
	static double w_g;
	static double w_l;
};

double Zad::epsilon = 1;
double Zad::TOL = 1.0e-8;
double Zad::delta = 0.1;
double Zad::n_x = 150;
double Zad::n_y = 100;
double Zad::V1 = 10;
double Zad::V2 = 0;
double Zad::w_g = 1;
double Zad::w_l = 1;

double rozklad_gestosci_ladunku_1(double x, double y, double x_max, double y_max);
double rozklad_gestosci_ladunku_2(double x, double y, double x_max, double y_max);

void nadaj_warunki_brzegowe_N(Matrix& mat);
void uwzglednij_warunki_brzegowe_N(Matrix& mat);

void relaksacja_globalna(Matrix& siatka_p, Matrix& siatka_v);
void relaksacja_lokalna(Matrix& siatka_p, Matrix& siatka_v);

void blad_relaksacji(Matrix& siatka_p, Matrix& siatka_v, Matrix& siatka_b);

void save_to_file(Matrix mat, std::string name);
void save_vec_to_file(std::vector<double> vec, std::string name);


int main() {


	Matrix siatka_p;
	Matrix siatka_v;
	Matrix siatka_b;
	Matrix siatka_czysta;

	std::vector<double> wiersz;
	wiersz.resize( Zad::n_y, 0 );
	siatka_p.resize( Zad::n_x, wiersz );
	siatka_v = siatka_p;
	siatka_b = siatka_p;
	siatka_czysta = siatka_p;

	double x_max = Zad::delta * Zad::n_x;
	double y_max = Zad::delta * Zad::n_y;
	for(int i = 0; i < Zad::n_x - 1; ++i )
		for(int j = 0; j < Zad::n_y - 1; ++j) {
			siatka_p[i][j] = rozklad_gestosci_ladunku_1(i * Zad::delta, j * Zad::delta, x_max, y_max) + rozklad_gestosci_ladunku_2(i * Zad::delta, j * Zad::delta, x_max, y_max);
		}

	nadaj_warunki_brzegowe_N(siatka_v);
	relaksacja_globalna(siatka_p, siatka_v);
	save_to_file(siatka_v, "V_g.dat");
	//blad_relaksacji(siatka_p, siatka_v, siatka_b);
	//save_to_file(siatka_b, "err_g.dat");

	siatka_v = siatka_czysta;

	nadaj_warunki_brzegowe_N(siatka_v);
	relaksacja_lokalna(siatka_p, siatka_v);
	save_to_file(siatka_v, "V_l_1_0.dat");

	siatka_v = siatka_czysta;
	Zad::w_l = 1.4;

	nadaj_warunki_brzegowe_N(siatka_v);
	relaksacja_lokalna(siatka_p, siatka_v);
	save_to_file(siatka_v, "V_l_1_4.dat");

	siatka_v = siatka_czysta;
	Zad::w_l = 1.8;

	nadaj_warunki_brzegowe_N(siatka_v);
	relaksacja_lokalna(siatka_p, siatka_v);
	save_to_file(siatka_v, "V_l_1_8.dat");

	siatka_v = siatka_czysta;
	Zad::w_l = 1.9;

	nadaj_warunki_brzegowe_N(siatka_v);
	relaksacja_lokalna(siatka_p, siatka_v);
	save_to_file(siatka_v, "V_l_1_9.dat");



	return 0;
}

double rozklad_gestosci_ladunku_1(double x, double y, double x_max, double y_max) {
	return exp( -pow(x - 0.35 * x_max, 2) / pow(0.1 * x_max, 2) - pow(y - 0.5 * y_max, 2) / pow(0.1 * y_max, 2) );
}

double rozklad_gestosci_ladunku_2(double x, double y, double x_max, double y_max) {
	return - exp( -pow(x - 0.65 * x_max, 2) / pow(0.1 * x_max, 2) - pow(y - 0.5 * y_max, 2) / pow(0.1 * y_max, 2) );
}

void relaksacja_globalna(Matrix& siatka_p, Matrix& siatka_v) {

	Matrix siatka_v_s = siatka_v;
	double S_1 = 1;
	double S_2 = 2;

	while( fabs( (S_2 - S_1)/S_1 ) > Zad::TOL ) {

		S_1 = S_2;

		for( int i = 1; i < Zad::n_x - 1; ++i ) 
			for( int j = 1; j < Zad::n_y - 1; ++j ) {
				siatka_v[i][j] = 0.25 * (siatka_v_s[i-1][j] + siatka_v_s[i+1][j]);
				siatka_v[i][j] += 0.25 * (siatka_v_s[i][j-1] + siatka_v_s[i][j+1]);
				siatka_v[i][j] += 0.25 * pow(Zad::delta, 2) * siatka_p[i][j];
			}

		uwzglednij_warunki_brzegowe_N(siatka_v);

		for( int i = 1; i < Zad::n_x - 1; ++i ) 
			for( int j = 1; j < Zad::n_y - 1; ++j ) {
				siatka_v_s[i][j] = siatka_v[i][j] * Zad::w_g + siatka_v_s[i][j] * (1 - Zad::w_g); 	
			}

		S_2 = 0;
		for( int i = 0; i < Zad::n_x - 1; ++i ) 
			for( int j = 0; j < Zad::n_y - 1; ++j ) {
				S_2 += 0.5 * pow( (siatka_v_s[i][j] - siatka_v_s[i+1][j]) / Zad::delta, 2 );
				S_2 += 0.5 * pow( (siatka_v_s[i][j] - siatka_v_s[i][j+1]) / Zad::delta, 2 );
				S_2 -= siatka_p[i][j] * siatka_v_s[i][j];	
			}
		S_2 *= pow(Zad::delta, 2);
		//std::cout << S_2 << " " << S_1 << "    " <<  fabs( (S_2 - S_1)/S_1 ) << std::endl;
	}
	siatka_v = siatka_v_s;
}

void relaksacja_lokalna(Matrix& siatka_p, Matrix& siatka_v) {

	double S_1 = 1;
	double S_2 = 2;
	std::vector<double> S_vec;

	while( fabs( (S_2 - S_1)/S_1 ) > Zad::TOL ) {

		S_1 = S_2;
		
		for( int i = 1; i < Zad::n_x - 1; ++i ) 
			for( int j = 1; j < Zad::n_y - 1; ++j ) {
				siatka_v[i][j] = (1 - Zad::w_l) * siatka_v[i][j] + 0.25 * Zad::w_l * (siatka_v[i-1][j] + siatka_v[i+1][j] + siatka_v[i][j-1] + siatka_v[i][j+1] + pow(Zad::delta, 2) * siatka_p[i][j]);
			}

		uwzglednij_warunki_brzegowe_N(siatka_v);

		S_2 = 0;
		for( int i = 0; i < Zad::n_x - 1; ++i ) 
			for( int j = 0; j < Zad::n_y - 1; ++j ) {
				S_2 += 0.5 * pow( (siatka_v[i][j] - siatka_v[i+1][j]) / Zad::delta, 2 );
				S_2 += 0.5 * pow( (siatka_v[i][j] - siatka_v[i][j+1]) / Zad::delta, 2 );
				S_2 -= siatka_p[i][j] * siatka_v[i][j];	
			}
		S_2 *= pow(Zad::delta, 2);
		S_vec.push_back(S_2);	
	}
	save_vec_to_file(S_vec, "S_it.dat");

}

void uwzglednij_warunki_brzegowe_N(Matrix& mat) {

	for(int i = 1; i < Zad::n_y - 1; i++) {
		mat[0][i] = mat[1][i];
	}

	for(int i = 1; i < Zad::n_y - 1; i++) {
		mat[Zad::n_x - 1][i] = mat[Zad::n_x-2][i];
	}

}

void nadaj_warunki_brzegowe_N(Matrix& mat) {

	for(int i = 0; i < Zad::n_x - 1; ++i) {
		mat[i][0] = Zad::V1;
	}

	for(int i = 0; i < Zad::n_x - 1; ++i) {
		mat[i][ Zad::n_y-1 ] = Zad::V2;
	}
}

void blad_relaksacji(Matrix& siatka_p, Matrix& siatka_v, Matrix& siatka_b) {

	for( int i = 1; i < Zad::n_x - 1; ++i ) 
		for( int j = 1; j < Zad::n_y - 1; ++j ) {
			siatka_b[i][j] = (siatka_v[i-1][j] - 2*siatka_v[i][j] + siatka_v[i+1][j]) / pow(Zad::delta, 2);
			siatka_b[i][j] += (siatka_v[i][j-1] - 2*siatka_v[i][j] + siatka_v[i][j+1]) / pow(Zad::delta, 2);
			siatka_b[i][j] += siatka_p[i][j];
		}


}

void save_to_file(Matrix mat, std::string name) {

	std::fstream file;
	file.open(name, std::fstream::out);
		for(int i = 0; i < Zad::n_x - 1; ++i) {
			for (int j = 0; j < Zad::n_y - 1; ++j) {
				file << i*Zad::delta << " " << j*Zad::delta << " " << mat[i][j] << std::endl;
			}
			file << std::endl;
		}
	file.close();
}

void save_vec_to_file(std::vector<double> vec, std::string name) {
	std::fstream file;
	file.open(name, std::fstream::out);
		for (int j = 0; j < vec.size(); ++j) {
			file << j << " " << vec[j] << std::endl;
		}
	file.close();
}

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include <fstream>

typedef std::vector< std::vector<double> > Matrix;
typedef std::vector<double> Vector;

struct Zad {
	static double x_max;
	static double y_max;
	static double delta;
	static int n_y;
	static int n_x;
	static double TOL;
};

int Zad::n_x = 128;
int Zad::n_y = 128;
double Zad::delta = 0.2;
double Zad::TOL = 1.e-8;
double Zad::x_max = Zad::delta * Zad::n_x;
double Zad::y_max = Zad::delta * Zad::n_y;

//deklaracje funkcji
double V_b1(double y);
double V_b2(double x);
double V_b3(double y);
double V_b4(double y);

void nadaj_warunki_poczatkowe(Matrix& mat_v);
double parametr_warunku_stopu( Matrix& mat_v, int k );
void relaksacja_z_parametrem_zageszczenia( Matrix& mat_v, Vector& vec_calki, int k );

void save_mat_to_file(const Matrix& mat, std::string name);
void save_vec_to_file(const Vector& vec, std::string name);


int main() {

	std::vector<double> wiersz;
	wiersz.resize(Zad::n_y + 1, 0);

	Matrix mat_v;
	Matrix mat_back_up;
	mat_back_up.resize(Zad::n_x + 1, wiersz);
	Vector vec_calki;
	nadaj_warunki_poczatkowe(mat_back_up);

	mat_v = mat_back_up;
	relaksacja_z_parametrem_zageszczenia(mat_v, vec_calki, 16);
	save_vec_to_file(vec_calki, "calka.dat");

	return 0;
}


//definicje

double V_b1(double y) {
	return sin( M_PI * y / Zad::y_max );
}

double V_b2(double x) {
	return -sin( 2 * M_PI * x / Zad::x_max );
}

double V_b3(double y) {
	return sin( M_PI * y / Zad::y_max );
}

double V_b4(double x) {
	return sin( 2 * M_PI * x / Zad::x_max );
}

void nadaj_warunki_poczatkowe(Matrix& mat_v) {

	for(int j = 0; j < Zad::n_y; ++j)
		mat_v[0][j] = V_b1(j * Zad::delta);

	for(int j = 0; j < Zad::n_y; ++j)
		mat_v[Zad::n_x][j] = V_b3(j * Zad::delta);

	for(int i = 0; i < Zad::n_x; ++i)
		mat_v[i][Zad::n_y] = V_b2(i * Zad::delta);

	for(int i = 0; i < Zad::n_x; ++i)
		mat_v[i][0] = V_b4(i * Zad::delta);
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
				fprintf(file, "%5.1lf %5.1lf %15.8lf\n", i*Zad::delta, j*Zad::delta, mat[i][j]);
			}
			fprintf(file, "\n");
		}
	fclose(file);
}

double parametr_warunku_stopu( Matrix& mat_v, int k ) {

	Matrix mat_v2 = mat_v;
	double wynik_cal = 0;
	double temp;

	for(int i = 0; i < Zad::n_x-k; i+=k) {
		for(int j = 0; j < Zad::n_y-k; j+=k){

			temp = (mat_v2[i][j+k] - mat_v2[i][j]) / (2 * k * Zad::delta);
			temp += (mat_v2[i+k][j+k] - mat_v2[i+k][j]) / (2 * k * Zad::delta);
			wynik_cal += pow(k*Zad::delta, 2) / 2 * pow(temp, 2);

			temp = (mat_v2[i+k][j] - mat_v2[i][j]) / (2 * k * Zad::delta);
			temp += (mat_v2[i+k][j+k] - mat_v2[i][j+k]) / (2 * k * Zad::delta);
			wynik_cal += pow(k*Zad::delta, 2) / 2 * pow(temp, 2);
		}
	}

	return wynik_cal;
}

void relaksacja_z_parametrem_zageszczenia( Matrix& mat_v, Vector& vec_calki, int k ) {

	Matrix mat_v2 = mat_v;
	double S_1 = 1;
	double S_2 = 2;
	vec_calki.clear();


	while( k > 0 ) {

		S_1 = 1;
		S_2 = 2;

		while( fabs( (S_2 - S_1) / S_1 ) > Zad::TOL ) {

			for(int i = k; i <= Zad::n_x-k; i+=k) {
				for(int j = k; j <= Zad::n_y-k; j+=k){
					mat_v[i][j] = 0.25 * (mat_v2[i-k][j] + mat_v2[i+k][j] + mat_v2[i][j-k] + mat_v2[i][j+k]);
				}
			}
			mat_v2 = mat_v;


			S_1 = S_2;
			S_2 = parametr_warunku_stopu(mat_v, k);
			vec_calki.push_back(S_2);
		}

		std::string mat_file_name = "mat_" + std::to_string(k) + ".dat";

		FILE* file;
		file = fopen(mat_file_name.c_str(), "w");
		for(int i = 0; i < mat_v.size(); i+=k){
			for(int j = 0; j < mat_v[0].size(); j+=k){
				fprintf(file, "%5.1lf %5.1lf %15.8lf\n", i*Zad::delta, j*Zad::delta, mat_v[i][j]);
			}
			fprintf(file, "\n");
		}
		fclose(file);

		if( k > 1 ) {
			for(int i = 0; i <= Zad::n_x-k; i+=k) {
				for(int j = 0; j <= Zad::n_y-k; j+=k) {


					if( i > 0 && i < Zad::n_x-k && j > 0 && j < Zad::n_y-k ) { //punkt wewnetrzny
						mat_v[i+k/2][j+k/2] = 0.25 * (mat_v[i][j] + mat_v[i+k][j+k] + mat_v[i][j+k] + mat_v[i+k][j]);
						mat_v[i+k][j+k/2] = 0.5 * (mat_v[i+k][j] + mat_v[i+k][j+k]);
						mat_v[i+k/2][j+k] = 0.5 * (mat_v[i][j+k] + mat_v[i+k][j+k]);
						mat_v[i+k/2][j] = 0.5 * (mat_v[i][j] + mat_v[i+k][j]);
						mat_v[i][j+k/2] = 0.5 * (mat_v[i][j] + mat_v[i][j+k]);
					}

					else if( i == 0 && j == 0 ) { //lewy gorny
						mat_v[i+k/2][j+k/2] = 0.25 * (mat_v[i][j] + mat_v[i+k][j+k] + mat_v[i][j+k] + mat_v[i+k][j]);
						mat_v[i+k][j+k/2] = 0.5 * (mat_v[i+k][j] + mat_v[i+k][j+k]);
						mat_v[i+k/2][j+k] = 0.5 * (mat_v[i][j+k] + mat_v[i+k][j+k]);
					}
					else if( i == 0 && j == Zad::n_y-k ) { //lewy dolny
						mat_v[i+k/2][j+k/2] = 0.25 * (mat_v[i][j] + mat_v[i+k][j+k] + mat_v[i][j+k] + mat_v[i+k][j]);
						mat_v[i+k/2][j+k] = 0.5 * (mat_v[i][j+k] + mat_v[i+k][j+k]);
						mat_v[i][j+k/2] = 0.5 * (mat_v[i][j] + mat_v[i][j+k]);
					}

					else if( i == Zad::n_x-k && j == 0 ) { //prawy gorny
						mat_v[i+k/2][j+k/2] = 0.25 * (mat_v[i][j] + mat_v[i+k][j+k] + mat_v[i][j+k] + mat_v[i+k][j]);
						mat_v[i+k][j+k/2] = 0.5 * (mat_v[i+k][j] + mat_v[i+k][j+k]);
						mat_v[i+k/2][j] = 0.5 * (mat_v[i][j] + mat_v[i+k][j]);

					}
					else if( i == Zad::n_x-k && j == Zad::n_y-k ) { //prawy dolny
						mat_v[i+k/2][j+k/2] = 0.25 * (mat_v[i][j] + mat_v[i+k][j+k] + mat_v[i][j+k] + mat_v[i+k][j]);
						mat_v[i+k/2][j] = 0.5 * (mat_v[i][j] + mat_v[i+k][j]);
						mat_v[i][j+k/2] = 0.5 * (mat_v[i][j] + mat_v[i][j+k]);
					}

					else if( i == 0 && j > 0 && j < Zad::n_y-k ) { //gora
						mat_v[i+k/2][j+k/2] = 0.25 * (mat_v[i][j] + mat_v[i+k][j+k] + mat_v[i][j+k] + mat_v[i+k][j]);
						mat_v[i+k/2][j+k] = 0.5 * (mat_v[i][j+k] + mat_v[i+k][j+k]);
						mat_v[i+k/2][j] = 0.5 * (mat_v[i][j] + mat_v[i+k][j]);
						mat_v[i+k][j+k/2] = 0.5 * (mat_v[i+k][j] + mat_v[i+k][j+k]);
					}
					else if( i == Zad::n_x-k && j > 0 && j < Zad::n_y-k ){ //dol
						mat_v[i+k/2][j+k/2] = 0.25 * (mat_v[i][j] + mat_v[i+k][j+k] + mat_v[i][j+k] + mat_v[i+k][j]);
						mat_v[i+k/2][j+k] = 0.5 * (mat_v[i][j+k] + mat_v[i+k][j+k]);
						mat_v[i+k/2][j] = 0.5 * (mat_v[i][j] + mat_v[i+k][j]);
						mat_v[i][j+k/2] = 0.5 * (mat_v[i][j] + mat_v[i][j+k]);
					}
					else if( i > 0 && i < Zad::n_x-k && j == 0 ) { //lewo
						mat_v[i+k/2][j+k/2] = 0.25 * (mat_v[i][j] + mat_v[i+k][j+k] + mat_v[i][j+k] + mat_v[i+k][j]);
						mat_v[i][j+k/2] = 0.5 * (mat_v[i][j] + mat_v[i][j+k]);
						mat_v[i+k][j+k/2] = 0.5 * (mat_v[i+k][j] + mat_v[i+k][j+k]);
						mat_v[i+k/2][j+k] = 0.5 * (mat_v[i][j+k] + mat_v[i+k][j+k]);
					}
					else if( i > 0 && i  < Zad::n_x-k && j == Zad::n_y-k ) {//prawo
						mat_v[i+k/2][j+k/2] = 0.25 * (mat_v[i][j] + mat_v[i+k][j+k] + mat_v[i][j+k] + mat_v[i+k][j]);
						mat_v[i][j+k/2] = 0.5 * (mat_v[i][j] + mat_v[i][j+k]);
						mat_v[i+k][j+k/2] = 0.5 * (mat_v[i+k][j] + mat_v[i+k][j+k]);
						mat_v[i+k/2][j] = 0.5 * (mat_v[i][j] + mat_v[i+k][j]);
					}
				} //for
			} //for

		} //if
		k = k / 2;
	} //while
}


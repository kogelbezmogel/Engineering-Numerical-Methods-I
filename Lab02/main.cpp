#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include <array>
#include <fstream>

typedef std::vector<double> vec;
typedef std::function<double(double, double)> iter_fun;
typedef std::array< std::array< double, 2 >, 2 > Mat2by2;

void z_t(vec& res_u_t, vec& res);
void save_fun_to_file(vec& args, vec& res, std::string name);

//Metoda Trapezow
double fun(double t, double u);
void metoda_trapezow(vec& args, vec& resul, iter_fun f);
double iteracja_picarda(double u_n, double t_n);
double iteracja_newtona(double u_n, double t_n);


//nRK2
void nRK2(vec& args, vec& resul);
void znajdz_U1_U2(double& U1, double& U2, double u_n, Mat2by2& a);
double F1(double U1, double U2, double u_n, Mat2by2& a);
double F2(double U1, double U2, double u_n, Mat2by2& a);

struct WarZad {
	static double beta;
	static double gama;
	static int N;
	static double t_max;
	static double d_t;
	static double TOL;
	static int ite_max;
	static int u_0;
};

double WarZad::beta = 0;
double WarZad::d_t = 0;
double WarZad::gama = 0;
int WarZad::N = 0;
int WarZad::u_0 = 0;
double WarZad::TOL = 0;
double WarZad::t_max = 0;
int WarZad::ite_max = 0;




int main() {

	WarZad::beta = 0.001;
	WarZad::d_t = 0.1;
	WarZad::gama = 0.1;
	WarZad::N = 500;
	WarZad::u_0 = 1;
	WarZad::TOL = 10e-6;
	WarZad::t_max = 100;
	WarZad::ite_max = 20;

	vec args;
	vec res_z;
	vec res_u;

	metoda_trapezow(args, res_u, iteracja_picarda);
	z_t(res_u, res_z);

	save_fun_to_file(args, res_u, "trapezow_pi_u.dat");
	save_fun_to_file(args, res_z, "trapezow_pi_z.dat");

	args.clear();
	res_z.clear();
	res_u.clear();

	metoda_trapezow(args, res_u, iteracja_newtona);
	z_t(res_u, res_z);

	save_fun_to_file(args, res_u, "trapezow_ne_u.dat");
	save_fun_to_file(args, res_z, "trapezow_ne_z.dat");

	args.clear();
	res_z.clear();
	res_u.clear();


	//nRK2

	nRK2(args, res_u);
	z_t(res_u, res_z);
	save_fun_to_file(args, res_u, "nrk2_u.dat");
	save_fun_to_file(args, res_z, "nrk2_z.dat");



	return 0;
}

void z_t(vec& res_u, vec& res_z) {

	for( auto i : res_u ) {
		res_z.push_back(WarZad::N - i);
	}
}

void save_fun_to_file(vec& args, vec& res, std::string name) {

    std::fstream file;
    file.open(name, std::fstream::out);

        for(int i = 0; i < args.size(); i++ ){
            file << std::to_string(args[i]) + " " + std::to_string(res[i]) << std::endl;
        }

    file.close();
}

double fun(double t, double u) {
	return (WarZad::beta * WarZad::N - WarZad::gama) * u - WarZad::beta * u * u;	
}

void metoda_trapezow(vec& args, vec& resul, iter_fun f_ite) {

	double u_n = WarZad::u_0;
	double u_n1 = u_n;
	double t = 0;
	args.push_back(0);
	resul.push_back(u_n);

	while(t < WarZad::t_max) {

	u_n1 = u_n + WarZad::d_t/2 * ( fun(t, u_n) + fun( t + WarZad::d_t, f_ite(u_n, t) ) );
	u_n = u_n1;
	t += WarZad::d_t;
	args.push_back(t);
	resul.push_back(u_n1);

	}
}

double iteracja_picarda(double u_n, double t_n) {
	double u_n1 = u_n;
	double u1_n1 = u_n;
	int ite = 0;

	while( abs(u1_n1 - u_n1) < WarZad::TOL && ite <= WarZad::ite_max) {
		u1_n1 = u_n + WarZad::d_t / 2 *  (fun(t_n, u_n) + fun(t_n + WarZad::d_t, u_n1));
		u_n1 = u1_n1;
		ite++;
	}

	return u1_n1;
}

double iteracja_newtona(double u_n, double t_n) {
	double u_n1 = u_n;
	double u1_n1 = u_n;
	int ite = 0;

	double temp = 0;
	double alpha = WarZad::beta * WarZad::N - WarZad::gama;

	while( abs(u1_n1 - u_n1) < WarZad::TOL && ite <= WarZad::ite_max ) {
		temp = alpha*u_n - WarZad::beta*u_n*u_n;
		temp += alpha*u_n1 - WarZad::beta*u_n1*u_n1;
		temp = u_n1 - u_n - WarZad::d_t/2 * temp;
		u1_n1 = u_n1 - temp / (1 - WarZad::d_t/2 * (alpha - 2*WarZad::beta*u_n1));
		u_n1 = u1_n1;
		ite++;
	}

	return u1_n1;
}

void nRK2(vec& args, vec& resul) {

	Mat2by2 a;
	std::array< double, 2 > c;
	std::array< double, 2 > b;

	c[0] = 0.5 - sqrt(3) / 6;
	c[1] = 0.5 + sqrt(3) / 6;

	b[0] = 0.5;
	b[1] = 0.5;

	a[0][0] = 0.25;
	a[0][1] = 0.25 - sqrt(3) / 6; 
	a[1][0] = 0.25 + sqrt(3) / 6;
	a[1][1] = 0.25;

	double t = 0;
	double u_n = WarZad::u_0;
	double u_n1 = WarZad::u_0;
	double U1 = u_n;
	double U2 = u_n;


	args.push_back(t);
	resul.push_back(u_n);

	while( t < WarZad::t_max ) {

		znajdz_U1_U2(U1, U2, u_n, a);
		u_n1 = u_n + WarZad::d_t * (b[0] * fun(t+c[0]*WarZad::d_t, U1) + b[1] * fun(t+c[1]*WarZad::d_t, U2) );
		u_n = u_n1;
		t += WarZad::d_t;

		args.push_back(t);
		resul.push_back(u_n1);
	}


}

void znajdz_U1_U2(double& U1, double& U2, double u_n, Mat2by2& a) {
	
	double dU1 = 1;
	double dU2 = 1;
	int ite = 0;
	Mat2by2 m;
	double alpha = WarZad::beta * WarZad::N - WarZad::gama;
	U1 = u_n;
	U2 = u_n;

	while( (dU1 > WarZad::TOL || dU2 > WarZad::TOL) && ite <= WarZad::ite_max ) {

	m[0][0] = 1 - a[0][0] * (alpha - 2 * WarZad::beta * U1);
	m[0][1] = -a[0][1] * (alpha - 2 * WarZad::beta * U2);
	m[1][0] = -a[1][0] * (alpha - 2 * WarZad::beta * U1);
	m[1][1] = 1 - a[1][1] * (alpha - 2 * WarZad::beta * U2);

	dU1 = F2(U1, U2, u_n, a) * m[0][1] - F1(U1, U2, u_n, a) * m[1][1];
	dU1 /= m[0][0] * m[1][1] - m[0][1] * m[1][0]; 

	dU2 = F1(U1, U2, u_n, a) * m[1][0] - F2(U1, U2,u_n, a) * m[0][0];
	dU2 /= m[0][0] * m[1][1] - m[0][1] * m[1][0]; 

	U1 += dU1;
	U2 += dU2;
	ite++;
	
	}
}

double F1(double U1, double U2, double u_n, Mat2by2& a) {
	double alpha = WarZad::beta * WarZad::N - WarZad::gama;
	double F = a[0][0] * (alpha * U1 - WarZad::beta * U1 * U1);
	F += a[0][1] * (alpha * U2  - WarZad::beta * U2 * U2);
	F = U1 - u_n - WarZad::d_t * F;
	return  F;
}

double F2(double U1, double U2, double u_n, Mat2by2& a) {
	double alpha = WarZad::beta * WarZad::N - WarZad::gama;
	double F = a[1][0] * (alpha * U1 - WarZad::beta * U1 * U1);
	F += a[1][1] * (alpha * U2  - WarZad::beta * U2 * U2);
	F = U2 - u_n - WarZad::d_t * F;
	return F;
}

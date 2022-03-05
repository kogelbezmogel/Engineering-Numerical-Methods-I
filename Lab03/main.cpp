#include <cmath>
#include <functional>
#include <array>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>

typedef std::function<double(double, double, double, double )> Fun4args;
typedef std::function< void (double x_n, double v_n, double d_t, double alpha, double* x_n1, double* v_n1) > SchematNumeryczny;
typedef std::vector<double> Vec;

//Schematy
void RK2( double x_n, double v_n, double d_t, double alpha, double* x_n1, double* v_n1 );
void MetodaTrapezow( double x_n, double v_n, double d_t, double alpha, double* x_n1, double* v_n1 );

//Kontrola kroku
void KontrolaKroku(Vec& args, Vec& res_x, Vec& res_v, Vec& res_dt, SchematNumeryczny schemat);
double ZmianaKroku(double d_t, double E_x, double E_v);

void save_fun_to_file(Vec& args, Vec& res, std::string name);

struct Zad{
	static double fun_f( double x_n, double v_n );
	static double fun_g( double x_n, double v_n );
	static double alpha;
	static double TOL;
	static double S;
	static double d_t0;
	static double t_0;
	static double x_0;
	static double v_0;
	static double t_end;
	static int p;
};
double Zad::alpha = 5;
double Zad::TOL = 10e-2;
double Zad::S = 0.75;
double Zad::d_t0 = 1;
double Zad::t_0 = 0;
double Zad::x_0 = 0.01;
double Zad::v_0 = 0;
double Zad::t_end = 40;
int Zad::p = 2;


int main() {
	
	Vec args;
	Vec res_x;
	Vec res_v;
	Vec res_dt;

	Zad::TOL = 10e-2;
	KontrolaKroku(args, res_x, res_v, res_dt, RK2);
	save_fun_to_file(args, res_x, "x_t1.dat");
	save_fun_to_file(args, res_v, "v_t1.dat");
	save_fun_to_file(args, res_dt, "dt_t1.dat");

	Zad::TOL = 10e-5;
	KontrolaKroku(args, res_x, res_v, res_dt, RK2);
	save_fun_to_file(args, res_x, "x_t2.dat");
	save_fun_to_file(args, res_v, "v_t2.dat");
	save_fun_to_file(args, res_dt, "dt_t2.dat");

	Zad::TOL = 1.e-2;
	KontrolaKroku(args, res_x, res_v, res_dt, MetodaTrapezow);
	save_fun_to_file(args, res_x, "x_t3.dat");
	save_fun_to_file(args, res_v, "v_t3.dat");
	save_fun_to_file(args, res_dt, "dt_t3.dat");

	Zad::TOL = 1.e-5;
	KontrolaKroku(args, res_x, res_v, res_dt, MetodaTrapezow);
	save_fun_to_file(args, res_x, "x_t4.dat");
	save_fun_to_file(args, res_v, "v_t4.dat");
	save_fun_to_file(args, res_dt, "dt_t4.dat");

}

////////definicje

void RK2( double x_n, double v_n, double d_t, double alpha, double* ref_x_n1, double* ref_v_n1 ) {

	double k_1x = Zad::fun_f(x_n, v_n);
	double k_1v = Zad::fun_g(x_n, v_n);

	double k_2x = Zad::fun_f(x_n + d_t*k_1x, v_n + d_t*k_1v);
	double k_2v = Zad::fun_g(x_n + d_t*k_1x, v_n + d_t*k_1v);
	
	*ref_x_n1 = x_n + d_t/2 * (k_1x + k_2x);
	*ref_v_n1 = v_n + d_t/2 * (k_1v + k_2v);
}

void MetodaTrapezow( double x_n, double v_n, double d_t, double alpha, double* ref_x_n1, double* ref_v_n1 ) {

	double x_n1 = x_n;
	double v_n1 = v_n;

	double d_x = 1;
	double d_v = 1;

	std::array< std::array<double, 2>, 2> a;

	Fun4args Fun_F = [d_t] (double x_n, double v_n, double x_n1, double v_n1) {
		return x_n1 - x_n - d_t/2 * (Zad::fun_f(x_n, v_n) + Zad::fun_f(x_n1, v_n1));
	};

	Fun4args Fun_G = [d_t] (double x_n, double v_n, double x_n1, double v_n1) {
		return v_n1 - v_n - d_t/2 * (Zad::fun_g(x_n, v_n) + Zad::fun_g(x_n1, v_n1));
	};


	while( fabs(d_x) > 1.e-10 || fabs(d_v) > 1.0e-10 ) {
		a[0][0]	= 1;
		a[0][1] = -d_t/2;
		a[1][0] = -d_t/2 *( (-2)* Zad::alpha * x_n1 * v_n1 - 1 ); 
		a[1][1] = 1 - d_t/2 * Zad::alpha * (1 - pow(x_n1,2));

		d_x = (-Fun_F(x_n, v_n, x_n1, v_n1) * a[1][1] + Fun_G(x_n, v_n, x_n1, v_n1) * a[0][1]) / (a[0][0]*a[1][1] - a[0][1]*a[1][0]);
		d_v = (-Fun_G(x_n, v_n, x_n1, v_n1) * a[0][0] + Fun_F(x_n, v_n, x_n1, v_n1) * a[1][0]) / (a[0][0]*a[1][1] - a[0][1]*a[1][0]);

		x_n1 += d_x;
		v_n1 += d_v;

	}

	*ref_x_n1 = x_n + d_t/2 * (Zad::fun_f(x_n, v_n) + Zad::fun_f(x_n1, v_n1));
	*ref_v_n1 = v_n + d_t/2 * (Zad::fun_g(x_n, v_n) + Zad::fun_g(x_n1, v_n1));
}

void KontrolaKroku(Vec& args, Vec& res_x, Vec& res_v, Vec& res_dt, SchematNumeryczny schemat) {

	args.clear();
	res_x.clear();
	res_v.clear();
	res_dt.clear();

	args.push_back(Zad::t_0);
	res_x.push_back(Zad::x_0);
	res_v.push_back(Zad::v_0);
	res_dt.push_back(Zad::d_t0);

	double d_t = Zad::d_t0;
	double t = Zad::t_0;
	double t_end = Zad::t_end;

	double E_x;
	double E_v;

	double x0 = Zad::x_0;
	double v0 = Zad::v_0;

	double x1_sdt = Zad::x_0;
	double v1_sdt = Zad::v_0;

	double x1_ddt = Zad::x_0;
	double v1_ddt = Zad::v_0;

	do {

		x1_ddt = x0;
		x1_sdt = x0;

		v1_ddt = v0;
		v1_sdt = v0;

		schemat(x1_sdt, v1_sdt, d_t, Zad::alpha, &x1_sdt, &v1_sdt);
		schemat(x1_sdt, v1_sdt, d_t, Zad::alpha, &x1_sdt, &v1_sdt);

		schemat(x1_ddt, v1_ddt, 2*d_t, Zad::alpha, &x1_ddt, &v1_ddt);

		E_x = ( x1_sdt - x1_ddt ) / ( pow(2, Zad::p) - 1 );
		E_v = ( v1_sdt - v1_ddt ) / ( pow(2, Zad::p) - 1 );
		E_x = fabs(E_x);
		E_v = fabs(E_v);

		if( std::max(E_v, E_x) < Zad::TOL ) {

			t += 2*d_t;

			res_dt.push_back(d_t);
			res_x.push_back(x1_sdt);
			res_v.push_back(v1_sdt);
			args.push_back(t);

			x0 = x1_sdt;
			v0 = v1_sdt;
		}

		d_t = ZmianaKroku(d_t, E_x, E_v);

	} while( t < t_end );

}

double ZmianaKroku(double d_t, double E_x, double E_v) {
	return pow( Zad::S*Zad::TOL / std::max(E_x, E_v) , 1.0/(Zad::p+1) ) * d_t;
}

double Zad::fun_f( double x_n, double v_n ) {
	return v_n;
}

double Zad::fun_g( double x_n, double v_n ) {
	return Zad::alpha * (1 - pow(x_n, 2)) * v_n - x_n;
}

void save_fun_to_file(Vec& args, Vec& res, std::string name) {

    std::fstream file;
    file.open(name, std::fstream::out);

        for(int i = 0; i < args.size(); i++ ){
            file << args[i] <<  " " << res[i] << std::endl;
        }

    file.close();
}

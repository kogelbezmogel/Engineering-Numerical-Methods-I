#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <functional>

typedef std::function< double(double t, double y, double param) > BaseFun;
typedef std::vector<double> VecD;

//funkcje do pierwszego zadania
double fun(double t, double y, double lambda);
double fun_org( double y, double lambda );
void euler(VecD& args, VecD& res, double d_t, double t_s, double t_e, double y_0, double lambda, BaseFun fun);
void RK2(VecD& args, VecD& res, double d_t, double t_s, double t_e, double y_0, double lambda, BaseFun fun);
void RK4(VecD& args, VecD& res, double d_t, double t_s, double t_e, double y_0, double lambda, BaseFun fun);
void save_fun_to_file(VecD& args, VecD& res, std::string name);
void save_err_to_file(VecD& args, VecD& res, std::string name, double lambda);
//


//funkcje i klasa do drugiego zadania
double fun_f(double t, double q, double i);
double fun_g(double t, double q, double i);
double V_t(double t);
void RK4_v2(VecD& args, VecD& res_q, VecD& res_i, double d_t, double t_s, double t_e, double q_0, double i_0, BaseFun funf, BaseFun fung);

struct Oscylator {
    static double R;
    static double L;
    static double C;
    static double w_0;
    static double w_v;
    static double T;
};

double Oscylator::R = 100.0;
double Oscylator::L = 0.1;
double Oscylator::C = 0.001;
double Oscylator::w_0 = 1/sqrt(Oscylator::L * Oscylator::C);
double Oscylator::T = 2 * M_PI / Oscylator::w_0;
double Oscylator::w_v = 0;
//



int main() {

    double y_0 = 1;
    double y = y_0;
    double lambda = -1;
    double t_s = 0;
    double t_e = 5;

    double d_t1 = 0.01;
    double d_t2 = 0.1;
    double d_t3 = 1;

    VecD results1;
    VecD args1;
    
    VecD results2;
    VecD args2;
    
    VecD results3;
    VecD args3;


    euler(args1, results1, d_t1, t_s, t_e, y_0, lambda, fun);
    euler(args2, results2, d_t2, t_s, t_e, y_0, lambda, fun);
    euler(args3, results3, d_t3, t_s, t_e, y_0, lambda, fun);

    save_fun_to_file(args1, results1, "results1_1.dat");
    save_fun_to_file(args2, results2, "results1_2.dat");
    save_fun_to_file(args3, results3, "results1_3.dat");

    save_err_to_file(args1, results1, "global1_1.dat", lambda);
    save_err_to_file(args2, results2, "global1_2.dat", lambda);
    save_err_to_file(args3, results3, "global1_3.dat", lambda);


    // metoda jawna RK2

    args1.clear();
    args2.clear();
    args3.clear();
    results1.clear();
    results2.clear();
    results3.clear();

    y = y_0;
    double k_1;
    double k_2;

    RK2(args1, results1, d_t1, t_s, t_e, y_0, lambda, fun);
    RK2(args2, results2, d_t2, t_s, t_e, y_0, lambda, fun);
    RK2(args3, results3, d_t3, t_s, t_e, y_0, lambda, fun);

    save_fun_to_file(args1, results1, "results2_1.dat");
    save_fun_to_file(args2, results2, "results2_2.dat");
    save_fun_to_file(args3, results3, "results2_3.dat");

    save_err_to_file(args1, results1, "global2_1.dat", lambda);
    save_err_to_file(args2, results2, "global2_2.dat", lambda);
    save_err_to_file(args3, results3, "global2_3.dat", lambda);


    //metoda RK4

    args1.clear();
    args2.clear();
    args3.clear();
    results1.clear();
    results2.clear();
    results3.clear();

    RK4(args1, results1, d_t1, t_s, t_e, y_0, lambda, fun);
    RK4(args2, results2, d_t2, t_s, t_e, y_0, lambda, fun);
    RK4(args3, results3, d_t3, t_s, t_e, y_0, lambda, fun);

    save_fun_to_file(args1, results1, "results3_1.dat");
    save_fun_to_file(args2, results2, "results3_2.dat");
    save_fun_to_file(args3, results3, "results3_3.dat");

    save_err_to_file(args1, results1, "global3_1.dat", lambda);
    save_err_to_file(args2, results2, "global3_2.dat", lambda);
    save_err_to_file(args3, results3, "global3_3.dat", lambda);


    args1.clear();
    results1.clear();
    double d_t = 10e-4;
    for(double i = t_s; i < t_e; i+=d_t) {
        args1.push_back(i);
        results1.push_back( fun_org(i, lambda) );
    }
    save_fun_to_file(args1, results1, "fun_org.dat");


    //zadanie 2

    double q_0;
    double i_0;

    t_s = 0;
    t_e = 4 * (2 * M_PI / Oscylator::w_0);
    d_t = 10e-4;
    q_0 = 0;
    i_0 = 0;

    args1.clear();
    results1.clear();
    results2.clear();
    Oscylator::w_v = 0.5 * Oscylator::w_0;
    RK4_v2(args1, results1, results2, d_t, t_s, t_e, q_0, i_0, fun_f, fun_g);
    save_fun_to_file(args1, results1, "results4q_1.dat");
    save_fun_to_file(args1, results2, "results4i_1.dat");


    args1.clear();
    results1.clear();
    results2.clear();
    Oscylator::w_v = 0.8 * Oscylator::w_0;
    RK4_v2(args1, results1, results2, d_t, t_s, t_e, q_0, i_0, fun_f, fun_g);
    save_fun_to_file(args1, results1, "results4q_2.dat");
    save_fun_to_file(args1, results2, "results4i_2.dat");


    args1.clear();
    results1.clear();
    results2.clear();
    Oscylator::w_v = Oscylator::w_0;
    RK4_v2(args1, results1, results2, d_t, t_s, t_e, q_0, i_0, fun_f, fun_g);
    save_fun_to_file(args1, results1, "results4q_3.dat");
    save_fun_to_file(args1, results2, "results4i_3.dat");


    args1.clear();
    results1.clear();
    results2.clear();
    Oscylator::w_v = 1.2 * Oscylator::w_0;
    RK4_v2(args1, results1, results2, d_t, t_s, t_e, q_0, i_0, fun_f, fun_g);
    save_fun_to_file(args1, results1, "results4q_4.dat");
    save_fun_to_file(args1, results2, "results4i_4.dat");

}



//-------------------------------------------------------------------definicje


double fun(double t, double y, double lambda){
    return lambda*y;
}

double fun_org( double y, double lambda ) {
    return exp(y*lambda);
}

void euler(VecD& args, VecD& res, double d_t, double t_s, double t_e, double y_0, double lambda, BaseFun fun) {
    
    res.push_back(y_0);
    args.push_back(t_s);

    int steps = (t_e - t_s) / d_t;
    double y = y_0;
    double t = t_s;

    for(int i = 0; i < steps; i++) {
        t += d_t;
        args.push_back(t);
        y = y + d_t * fun(t, y, lambda);
        res.push_back( y );
    }
}

void RK2(VecD& args, VecD& res, double d_t, double t_s, double t_e, double y_0, double lambda, BaseFun fun) {

    res.push_back(y_0);
    args.push_back(t_s);
    double y = y_0;
    double t = t_s;
    double k_1;
    double k_2;

    int steps = (t_e - t_s) / d_t;
    for(int i = 0; i < steps; i++) {

        k_1 = fun( t, y, lambda );
        k_2 = fun( t + d_t, y + d_t*k_1, lambda );

        t += d_t;
        args.push_back( t );
        y = y + d_t/2 * (k_1 + k_2);
        res.push_back( y );
    }
}

void RK4(VecD& args, VecD& res, double d_t, double t_s, double t_e, double y_0, double lambda, BaseFun fun) {

    res.push_back(y_0);
    args.push_back(t_s);
    double y = y_0;
    double t = t_s;

    double k_1;
    double k_2;
    double k_3;
    double k_4;
    
    int steps = (t_e - t_s) / d_t;
    for(int i = 0; i < steps; i++) {

        k_1 = fun( t, y, lambda );
        k_2 = fun( t + 0.5*d_t, y + d_t/2*k_1, lambda );
        k_3 = fun( t + 0.5*d_t, y + d_t/2*k_2, lambda );
        k_4 = fun( t + d_t, y + d_t*k_3, lambda );

        t += d_t;
        args.push_back( t );
        y = y + d_t/6 * (k_1 + 2*k_2 + 2*k_3 + k_4);
        res.push_back( y );
    }
}

void save_fun_to_file(VecD& args, VecD& res, std::string name) {

    std::fstream file;
    file.open(name, std::fstream::out);

        for(int i = 0; i < args.size(); i++ ){
            file << std::to_string(args[i]) + " " + std::to_string(res[i]) << std::endl;
        }

    file.close();
}

void save_err_to_file(VecD& args, VecD& res, std::string name, double lambda) {

    std::fstream file;
    file.open(name, std::fstream::out);

        for(int i = 0; i < args.size(); i++ ){
            file << args[i] << " " << std::abs(res[i] - fun_org(args[i], lambda))  << std::endl;
        }

    file.close();
}             

double fun_f(double t, double q, double i) {
  return i;
}

double fun_g(double t, double q, double i) {
  double to_return = V_t(t) / Oscylator::L;
  to_return -= Oscylator::R * i / Oscylator::L;
  to_return -= q / (Oscylator::C * Oscylator::L);
  return to_return;
}

double V_t(double t) {
  return sin(Oscylator::w_v * t);
}

void RK4_v2(VecD& args, VecD& res_q, VecD& res_i, double d_t, double t_s, double t_e, double q_0, double i_0, BaseFun funf, BaseFun fung) {

  args.push_back(t_s);
  res_q.push_back(q_0);
  res_i.push_back(i_0);

  double i = i_0;
  double q = q_0;
  double t = t_s;
  int steps = (t_e - t_s) / d_t;

  for(int k = 0; k < steps; k++) {

    double kq_1 = funf(t, q, i);
    double ki_1 = fung(t, q, i);
    double kq_2 = funf(t + 0.5*d_t, q + d_t/2 * kq_1, i + d_t/2 * d_t*ki_1);
    double ki_2 = fung(t + 0.5*d_t, q + d_t/2 * kq_1, i + d_t/2 * d_t*ki_1);
    double kq_3 = funf(t + 0.5*d_t, q + d_t/2 * kq_2, i + d_t/2 * d_t*ki_2);
    double ki_3 = fung(t + 0.5*d_t, q + d_t/2 * kq_2, i + d_t/2 * d_t*ki_2);
    double kq_4 = funf(t + d_t, q + d_t * kq_3, i + d_t * ki_3);
    double ki_4 = fung(t + d_t, q + d_t * kq_3, i + d_t * ki_3);

    t += d_t;
    args.push_back( t );
    q = q + d_t/6 * (kq_1 + 2*kq_2 + 2*kq_3 + kq_4);
    i = i + d_t/6 * (ki_1 + 2*ki_2 + 2*ki_3 + ki_4);
    res_q.push_back( q ); 
    res_i.push_back( i );
  }
}
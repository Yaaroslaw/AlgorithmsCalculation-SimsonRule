#pragma once

#include <iomanip>
#include <fstream>




using namespace std;


const double a = 7, b = 20, eps = 10e-2;

using namespace std;

double my_func(double x);
double g(double t);
double Legendre_i(int n, double x);
double g_L(double t, int k);
double Simpson_rule(int k, int n);
double* get_coef(int n);
double result(double x, double* A, int n);

#include "Header.h"

double my_func(double x) {
	return  0.2*log(abs(sin(x / 3)) / 3) + 1;
}


double g(double t) {
	return my_func((b + a) / 2 + t * (b - a) / 2);
}


double Legendre_i(int n, double x) {
	double pi, p1 = x, p0 = 1;
	if (n == 0)
		return 1;
	if (n == 1)
		return x;
	int i = 1;
	while (i < n) {
		pi = (1.0 *(2 * i + 1) / (i + 1)) * x * p1 - (1.0* i / (i + 1)) * p0;
		p0 = p1;
		p1 = pi;
		++i;
	}
	return p1;
}


double g_L(double t, int k) {
	return g(t) * Legendre_i(k, t);
}


double Simpson_rule( int k, int n) {
	double h = 2.0 / n;
	double I, sigma1 = 0, sigma2 = 0;
	double x = -1;

	I = g_L(-1, k) + g_L(1, k);
	x += h;

	for (int i = 0; i < (n - 2) / 2; ++i) {
		sigma1 += g_L(x, k);
		x += h;
		sigma2 += g_L(x, k);
		x += h;
	}
	sigma1 += g_L(x, k);

	return (h / 3) * (I + 4 * sigma1 + 2 * sigma2);

}

double* get_coef(int n) {
	double* A = new  double[n];
	for (int i = 0; i < n; ++i) {
		A[i] = Simpson_rule(i, 5000)* (2 * i + 1) / 2;
	}
	return A;
}


double result(double x, double* A, int n) {
	double  t = (2 * x - (b + a)) / (b - a);
	double sum = 0;
	for (int j = 0; j < n; ++j) {
		sum += A[j] * Legendre_i(j, t);
	}
	return sum;
}



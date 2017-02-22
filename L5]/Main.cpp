
#include "Header.h"


int main() {
	int k = 100, n = 100;
	double step = (b - a) / k, x, y, approx_pres;
	double* coef_a;
	do {
		n++;
		coef_a = get_coef(n);
		approx_pres = 0;
		x = a;
		for (int i = 0; i <= k; ++i) {
			y = my_func(x) - result(x, coef_a, n);
			approx_pres += y * y;
			x += step;
		}
		approx_pres = sqrt(approx_pres) / (k + 1);

	} while (approx_pres > eps);

	ofstream fout1("1.txt");//my function
	x = a;
	for (int i = 0; i < k; i++) {
		fout1 << x << ";" << my_func(x) << "\n";
		x += step;
	}
	fout1.close();

	ofstream fout("2.txt"); //approximations
	x = a;
	for (int i = 0; i <= k; ++i) {
		y = result(x, coef_a, n);
		fout << x << " ; ";
		fout << setprecision(3) << y << "\n";
		x += step;
	}
	fout.close();

	return 0;
}
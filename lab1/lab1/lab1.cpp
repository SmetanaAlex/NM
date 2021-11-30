#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

using namespace std;

double f1(double x) {
	return x * x * x - 7 * x * x + 7 * x + 15;
}

double f2(double x) {
	return x * x * x - 5 * x * x - 4 * x + 20;
}

double f3(double x) {
	return x * x * x - 4 * x * x - 4 * x + 16;
}

double df1(double x) {
	return 3 * x*x - 14 * x + 7;
}

double df2(double x) {
	return 3 * x * x - 10 * x - 4;
}

double df3(double x) {
	return 3 * x * x - 8 * x - 4;
}

double ddf1(double x) {
	return 6 * x - 14;
}

double ddf2(double x) {
	return 6 * x - 10;
}

double ddf3(double x) {
	return 6 * x - 8;
}

bool relaxMethod(double(f)(double), double(df)(double), double& a, double b, double eps = 0.001)
{
	double m, M, tau;
	double prev_x;
	if (abs(df(a)) > abs(df(b)))
	{
		m = abs(df(b));
		M = abs(df(a));
	}
	else
	{
		m = abs(df(a));
		M = abs(df(b));
	}

	tau = 2 / (M + m);
	double q = (M - m) / (M + m);
	
	if (tau * m <= 0 || tau * M >= 2 || q >= 1)
	{
		cout << "The condition of convergence are not met" << endl;
		return false;
	}

	double n = int(log(abs(a - b) / eps) / log(1 / q)) + 1;
	int i = 0;
	cout << "Step: " << setw(2) << i << "; Number: " << setprecision(12) << a << endl;
	do {
		i++;
		prev_x = a;
		a = prev_x - tau * f(prev_x);
		cout << "Step: " << setw(2) << i << "; Number: " << setprecision(12) << a << endl;
	} while (i < n);
	return true;
}

bool newtonMethod(double(f)(double), double(df)(double), double(ddf)(double),
	double& a, double b, double eps = 0.001)
{
	if (f(a) * f(b) > 0)
	{
		cout << "Function values have the same sign" << endl;
		return false;
	}
	if (df(a) * df(b) < 0)
	{
		cout << "The first derivative is not linear" << endl;
		return false;
	}
	if (ddf(a) * ddf(b) < 0)
	{
		cout << "Theorem 2.2 failed" << endl;
		return false;
	}
	double M = max(ddf(a), ddf(b));
	double m = min(df(a), df(b));
	double q = M * abs(a - b) / (2 * m);
	if (q >= 1)
	{
		cout << "q >= 1" << endl;
		return false;
	}
	int n = 1 + int(log2(1 + log(abs(a - b) / eps) / log(1 / q)));
	int i = 0;
	cout << "Step: " << setw(2) << i << "; Number: " << setprecision(12) << a << endl;
	do {
		i++;
		a = a - f(a) / df(a);
		cout << "Step: " << setw(2) << i << "; Number: " << setprecision(12) << a << endl;
	} while (i < n);

	return true;
}

bool crossLineMethod(double(f)(double), double(df)(double), double& a_1, double a_0, double b)
{
	if (f(a_1) * f(b) > 0)
	{
		cout << "Function values have the same sign" << endl;
		return false;
	}
	if (f(a_1) * f(a_0) < 0)
	{
		cout << "The solution can`t be between x_0 and x_1" << endl;
		return false;
	}
	if (df(a_1) * df(a_0) < 0)
	{
		cout << "The first derivative is not linear" << endl;
		return false;
	}
	if (df(a_1) * df(b) < 0)
	{
		cout << "The first derivative is not linear" << endl;
		return false;
	}
	int n = 10;
	int i = 0;
	do {
		i++;
		a_1 = a_1 - f(a_1)*(a_1-a_0) / (f(a_1) - f(a_0));
		cout << "Step: " << setw(2) << i << "; Number: " << setprecision(12) << a_1 << endl;
	} while (i < n);
	return true;
}

int main() {
	double a1 = -3.0;
	double a2 = -1.25;
	double a3_0 = -3.0;
	double a3_1 = -2.5;
	double b1 = -1.0;
	double b2 = -0.75;
	double b3 = -1;
	
	if (relaxMethod(f3, df3, a1, b1))
		cout << endl << "Relax method result:  x = " << a1 << "\n\n";

	if (newtonMethod(f1, df1, ddf1, a2, b2))
		cout << endl << "Newton method result:  x = " << a2 << "\n\n";

	if (crossLineMethod(f2, df2, a3_1, a3_0, b3))
		cout << endl << "Cross line method result:  x = " << a3_1 << "\n\n";
}
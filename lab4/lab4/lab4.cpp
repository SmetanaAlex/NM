#include <iostream>
#include <math.h>
#include <string>
#include <iomanip>
#include <vector>
#include <algorithm>

using namespace std;

double f1(double x)
{
    return 6 * pow(x, 14) + 8 * pow(x, 13) + 3 * pow(x, 12);
}

double f2(double x)
{
    return 6 * pow(x, 14) + 5 * pow(x, 11) + 7 * pow(x, 9);
}

double getCofficientX_0(double a, vector<double> x)
{
    double res = 1.0;
    for (auto i : x)
        if (a != i)
            res *= i;
    return res;
}

double getCofficientX_1(double a, vector<double> x)
{
    double res = 0.0;
    for (auto i : x)
        if (a != i)
        {
            double cur_res = 1.0;
            for (auto j : x)
                if (a != j && i != j)
                    cur_res *= j;
            res += cur_res;
        }
    return (-1) * res;
}

double getCofficientX_2(double a, vector<double> x)
{
    double res = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = i + 1; j < x.size(); j++)
            if (a != x[i] && a != x[j])
            {
                double cur_res = 1.0;
                for (auto k : x)
                    if (a != k && x[i] != k && x[j] != k)
                        cur_res *= k;
                res += cur_res;
            }
    return res;
}

double getCofficientX_3(double a, vector<double> x)
{
    double res = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = i + 1; j < x.size(); j++)
            for (size_t k = j + 1; k < x.size(); k++)
                if (a != x[i] && a != x[j] && a != x[k])
                {
                    double cur_res = 1.0;
                    for (auto m : x)
                        if (a != m && x[i] != m && x[j] != m && x[k] != m)
                            cur_res *= m;
                    res += cur_res;
                }
    return (-1) * res;
}

double getCofficientX_4(double a, vector<double> x)
{
    double res = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = i + 1; j < x.size(); j++)
            for (size_t k = j + 1; k < x.size(); k++)
                for (size_t m = k + 1; m < x.size(); m++)
                    if (a != x[i] && a != x[j] && a != x[k] && a != x[m])
                    {
                        double cur_res = 1.0;
                        for (auto n : x)
                            if (a != n && x[i] != n && x[j] != n && x[k] != n && x[m] != n)
                                cur_res *= n;
                        res += cur_res;
                    }
    return res;
}

double getCofficientX_5(double a, vector<double> x)
{
    double res = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = i + 1; j < x.size(); j++)
            for (size_t k = j + 1; k < x.size(); k++)
                for (size_t m = k + 1; m < x.size(); m++)
                    for (size_t n = m + 1; n < x.size(); n++)
                        if (a != x[i] && a != x[j] && a != x[k] && a != x[m] && a != x[n])
                        {
                            double cur_res = 1.0;
                            for (auto l : x)
                                if (a != l && x[i] != l && x[j] != l && x[k] != l && x[m] != l && x[n] != l)
                                    cur_res *= l;
                            res += cur_res;
                        }
    return (-1) * res;
}

double getCofficientX_6(double a, vector<double> x)
{
    double res = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = i + 1; j < x.size(); j++)
            for (size_t k = j + 1; k < x.size(); k++)
                for (size_t m = k + 1; m < x.size(); m++)
                    for (size_t n = m + 1; n < x.size(); n++)
                        for (size_t l = n + 1; l < x.size(); l++)
                            if (a != x[i] && a != x[j] && a != x[k] && a != x[m] && a != x[n] && a != x[l])
                                res += i * j * k * m * n * l;
    return res;
}

double getCofficientX_7(double a, vector<double> x)
{
    double res = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = i + 1; j < x.size(); j++)
            for (size_t k = j + 1; k < x.size(); k++)
                for (size_t m = k + 1; m < x.size(); m++)
                    for (size_t n = m + 1; n < x.size(); n++)
                        if (a != x[i] && a != x[j] && a != x[k] && a != x[m] && a != x[n])
                            res += i * j * k * m * n;
    return (-1) * res;
}

double getCofficientX_8(double a, vector<double> x)
{
    double res = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = i + 1; j < x.size(); j++)
            for (size_t k = j + 1; k < x.size(); k++)
                for (size_t m = k + 1; m < x.size(); m++)
                    if (a != x[i] && a != x[j] && a != x[k] && a != x[m])
                        res += i * j * k * m;
    return res;
}

double getCofficientX_9(double a, vector<double> x)
{
    double res = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = i + 1; j < x.size(); j++)
            for (size_t k = j + 1; k < x.size(); k++)
                if (a != x[i] && a != x[j] && a != x[k])
                    res += i * j * k;
    return (-1) * res;
}

double getCofficientX_10(double a, vector<double> x)
{
    double res = 0.0;
    for (size_t i = 0; i < x.size(); i++)
        for (size_t j = i + 1; j < x.size(); j++)
            if (a != x[i] && a != x[j])
                res += x[i] * x[j];
    return res;
}

double getCofficientX_11(double a, vector<double> x)
{
    double res = 0.0;
    for (auto i : x)
        if (a != i)
            res += i;
    return (-1) * res;
}

double getDenominator(double a, vector<double> x) 
{
    double res = 1.0;
    for (auto i : x)
        if (a != i)
            res *= (a - i);
    return res;
}

void printPolinom(const vector<double>& x)
{
    bool first = true;
    for (size_t i = 0; i < x.size(); i++)
    {
        if (x[i] != 0.0)
        {
            if (first)
            {
                first = false;
                if (x[i] == 1.0)
                    cout << "x^(" << to_string(x.size() - 1 - i) << ')';
                else if (x[i] == -1.0)
                    cout << "-x^(" << to_string(x.size() - 1 - i) << ')';
                else
                    cout << x[i] << "x^(" << to_string(x.size() - 1 - i) << ')';
            }
            else
            {
                if (x[i] > 0.0)
                    cout << " + ";
                else cout << " - ";
                if (abs(x[i]) == 1.0)
                    cout << "x^(" << to_string(x.size() - 1 - i) << ')';
                else
                    cout << abs(x[i]) << "x^(" << to_string(x.size() - 1 - i) << ')';
            }
        }
    }
    cout << endl;
}

void polinomLagrange(double(f)(double),
    double start, double interval,
    int numOfPoint)
{
    vector<double> x;
    for (int i = 0; i < numOfPoint; i++)
        x.push_back(start + i * interval);
    vector<double> coefficients(numOfPoint, 0.0);

    for (auto i : x)
    {
        double coefficient = f(i) / getDenominator(i, x);
        coefficients[0] += coefficient;
        coefficients[1] += coefficient * getCofficientX_11(i, x);
        coefficients[2] += coefficient * getCofficientX_10(i, x);
        coefficients[3] += coefficient * getCofficientX_9(i, x);
        coefficients[4] += coefficient * getCofficientX_8(i, x);
        coefficients[5] += coefficient * getCofficientX_7(i, x);
        coefficients[6] += coefficient * getCofficientX_6(i, x);
        coefficients[7] += coefficient * getCofficientX_5(i, x);
        coefficients[8] += coefficient * getCofficientX_4(i, x);
        coefficients[9] += coefficient * getCofficientX_3(i, x);
        coefficients[10] += coefficient * getCofficientX_2(i, x);
        coefficients[11] += coefficient * getCofficientX_1(i, x);
        coefficients[12] += coefficient * getCofficientX_0(i, x);
    }

    printPolinom(coefficients);
}

void printM(const vector<vector<double>>& A_matrix)
{
    for (auto i : A_matrix)
    {
        for (auto j : i)
            cout << setw(15) << setprecision(13) << j << ' ';
        cout << endl;
    }
}

vector<double> GaussMethod(vector<vector<double>>& A_matrix)
{
    for (size_t i = 0; i < A_matrix.size(); i++)
    {
        for (size_t j = i + 1; j < A_matrix.size() + 1; j++)
        {
            A_matrix[i][j] = A_matrix[i][j] / A_matrix[i][i];
            for (size_t k = i + 1; k < A_matrix.size(); k++)
            {
                A_matrix[k][j] = A_matrix[k][j] - A_matrix[k][i] * A_matrix[i][j];
            }
        }
        A_matrix[i][i] = 1;
        for (size_t k = i + 1; k < A_matrix.size(); k++)
        {
            A_matrix[k][i] = 0;
        }
        cout << "Step = " << i + 1 << endl;
        printM(A_matrix);
        cout << endl;
    }
    vector<double> result(A_matrix.size());

    for (int i = static_cast<int>(A_matrix.size()) - 1; i >= 0; i--)
    {
        result[i] = A_matrix[i][A_matrix.size()];
        for (size_t j = i + 1; j < A_matrix.size(); j++)
        {
            result[i] -= A_matrix[i][j] * result[j];
        }
    }
    return result;
}

void polinomPower(double(f)(double),
    double start, double interval,
    int numOfPoint)
{
    vector<vector<double>> A(numOfPoint, vector<double>(numOfPoint + 1, 0.0));
    vector<double> x;
    for (int i = 0; i < numOfPoint; i++)
    {
        for (int j = 0; j < numOfPoint; j++)
        {
            A[i][j] = pow(start + i * interval, j);
        }
        A[i][numOfPoint] = f(start + i * interval);
    }
    vector<double> cofficients = GaussMethod(A);
    reverse(cofficients.begin(), cofficients.end());
    printPolinom(cofficients);
}

int main()
{
    polinomLagrange(f1, -3.0, 0.5, 13);
    cout << endl << endl;
    polinomPower(f2, -3.0, 0.5, 12);
}


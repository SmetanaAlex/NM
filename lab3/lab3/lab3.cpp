#include <iostream>
#include <vector>
#include <math.h>
#include <functional>
#include <iomanip>

using namespace std;

double f1(double x, double y)
{
    return (1 + y * y) / y;
}

double f2(double x, double y)
{
    return 5 / (x * x + 1);
}


void simpleIterationMethod(double(f_x)(double, double), double(f_y)(double, double), double& x, double& y)
{
    double x_new = 0.0;
    double y_new = 0.0;
    double eps = 0.001;
    bool first = true;
    int n = 0;
    do
    {
        if (first)
        {
            first = false;
        }
        else
        {
            x = x_new;
            y = y_new;
        }
        x_new = f_x(x, y);
        y_new = f_y(x, y);
        cout << "Iteration n = " << n << endl;
        cout << "X = " << x << endl;
        cout << "Y = " << y << endl;

        n++;
    } while (n < 20);
    x = x_new;
    y = y_new;

}

bool isSymmetrical(const vector<vector<double>>& A)
{
    if (A.size() == 0)
        return false;
    for (int i = 0; i < static_cast<int>(A.size()); i++)
        if (A.size() != A[i].size())
            return false;
    for (int i = 0; i < static_cast<int>(A.size()); i++)
        for (int j = i; j < static_cast<int>(A.size()); j++)
            if (A[i][j] != A[j][i])
                return false;
    return true;
}

double normMatrixInfinity(const vector<vector<double>>& A)
{
    double m = 0.0;
    for (int i = 0; i < static_cast<int>(A.size()); i++)
    {
        double m_i = 0.0;
        for (int j = 0; j < static_cast<int>(A.size()); j++)
            m_i += abs(A[i][j]);
        if (m_i > m)
            m = m_i;
    }
    return m;
}

void normalizeVector(vector<double>& a)
{
    double dist = 0.0;
    for (int i = 0; i < static_cast<int>(a.size()); i++)
        dist += a[i] * a[i];
    dist = sqrt(dist);
    for (int i = 0; i < static_cast<int>(a.size()); i++)
        a[i] /= dist;
}

double powerMethod(const vector<vector<double>>& A)
{
    if (isSymmetrical(A))
    {
        double n = normMatrixInfinity(A);
        vector<vector<double>> B(A.size(), vector<double>(A.size(), 0.0));
        for (int i = 0; i < static_cast<int>(A.size()); i++)
        {
            B[i][i] = n;
            for (int j = 0; j < static_cast<int>(A.size()); j++)
                B[i][j] -= A[i][j];
        }
        vector<double> x(A.size(), 1.0);
        double m = 1000.0;
        double m_new = 999.0;
        double eps = 0.001;
        while(abs(m - m_new) > eps) 
        {
            normalizeVector(x);
            vector<double> x_new;
            for (int i = 0; i < static_cast<int>(A.size()); i++)
            {
                double x_i_new = 0.0;
                for (int j = 0; j < static_cast<int>(A.size()); j++)
                    x_i_new += B[i][j] * x[j];
                x_new.push_back(x_i_new);
            }
            m = m_new;
            m_new = 0;
            for (int j = 0; j < static_cast<int>(A.size()); j++)
                m_new += x_new[j] * x[j];
            x = x_new;
        }
        return n - m_new;
    }
}

int main()
{
    {
        double delta = 1;
        double i = 0.0;
        for (; i < 10.0; i += delta)
        {
            cout << "A = " << i << endl;
            vector<vector<double>> A = { {2.0, 1.0, 0.0},
                                         {1.0, 2.0, i},
                                         {0.0, i, 2.0} };
            double res = powerMethod(A);
            cout << "Result of algorithm: " << res << endl;
            if (res < 0.3)
            {
                i -= delta;
                delta /= 10;
                if (delta < 0.001)
                    break;
            }
        }
        cout << "A = " << i << endl;
        vector<vector<double>> A = { {2.0, 1.0, 0.0},
                                         {1.0, 2.0, i},
                                         {0.0, i, 2.0} };
        double res = powerMethod(A);
        cout << "Result of algorithm: " << res << endl;
    }/*/
    {
        double x = 0.12;
        double y = -0.705;
        simpleIterationMethod(f1_1, f2_2, x, y);
        cout << "X = " << x << endl;
        cout << "Y = " << y << endl;
    }
    {
        double x = -0.12;
        double y = 0.705;
        simpleIterationMethod(f1_1, f2_1, x, y);
        cout << "X = " << x << endl;
        cout << "Y = " << y << endl;
    }
    {
        double x = -1.0;
        double y = -0.5;
        simpleIterationMethod(f1_2, f2_2, x, y);
        cout << "X = " << x << endl;
        cout << "Y = " << y << endl;
    }
    {
        double x = 1.0;
        double y = 0.5;
        simpleIterationMethod(f1_1, f2_1, x, y);
        cout << "X = " << x << endl;
        cout << "Y = " << y << endl;
    }*/
    {
        double x = 4;
        double y = 2;
        simpleIterationMethod(f1, f2, x, y);
    }
}
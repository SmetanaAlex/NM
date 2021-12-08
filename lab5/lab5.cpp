#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace std;

double f(double x) {
    return pow(x, 8) - pow(x, 5) + 3 * pow(x, 3) - 2;
}

vector<double> calculatePart(double(f)(double), double x_start, double x_end)
{
    double denominator = x_end - x_start;
    double y_start = f(x_start);
    double y_end = f(x_end);
    vector<double> result;
    double coefficient_x1 = (y_end - y_start) / denominator;
    double coefficient_x0 = (y_start * x_end - y_end * x_start) / denominator;
    result.push_back(coefficient_x1);
    result.push_back(coefficient_x0);
    return result;
}

map<pair<double, double>, vector<double>> calculateSplaynLinear(double(f)(double),
    double x_start,
    double x_end,
    int number_of_nodes)
{
    map<pair<double, double>, vector<double>> result;
    double interval = (x_end - x_start + 1) / number_of_nodes;
    for (int i = 1; i < number_of_nodes; i++)
    {
        double x_s = x_start + (i - 1) * interval;
        double x_f = x_start + i * interval;
        result[make_pair(x_s, x_f)] = calculatePart(f, x_s, x_f);
    }
    return result;
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
}

void printResult(const map<pair<double, double>, vector<double>>& result)
{
    for (auto i : result)
    {
        printPolinom(i.second);
        cout << "\t\tx in [" << i.first.first << "; " << i.first.second << "];\n";
    }
}

int main()
{
    printResult(calculateSplaynLinear(f, 1.0, 6.0, 6));
}
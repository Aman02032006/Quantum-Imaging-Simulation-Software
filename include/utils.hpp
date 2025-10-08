#ifndef UTILS_H
#define UTILS_H

#pragma once
#include <cmath>
#include <cstdlib>
#include <limits>

using namespace std;

// Constants

const double INF = numeric_limits<double>::infinity();
const double PI = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees)
{
    return degrees * PI / 180.0;
}

// Generates random doubles in [0.0, 1.0)
inline double random_double()
{
    return std ::rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max)
{
    return min + (max - min) * random_double();
}

inline double sq(double a)
{
    return a * a;
}

enum class FieldType
{
    PLANE,
    GAUSSIAN,
    LG,
    HG
};

inline double factorial(int n)
{
    double res = 1.0;
    for (int i = 2; i <= n; ++i)
        res *= i;
    return res;
}

double genLaguerre(int p, int l, double x)
{
    double sum = 0.0;
    for (int m = 0; m <= p; ++m)
    {
        double term = pow(-1, m) * factorial(p + l) /
                      (factorial(p - m) * factorial(l + m) * factorial(m)) *
                      pow(x, m);
        sum += term;
    }
    return sum;
}

double hermitePol(int n, double x)
{
    if (n == 0)
        return 1.0;
    if (n == 1)
        return 2.0 * x;
    double Hnm2 = 1.0, Hnm1 = 2.0 * x, Hn;
    for (int i = 2; i <= n; ++i)
    {
        Hn = 2.0 * x * Hnm1 - 2.0 * (i - 1) * Hnm2;
        Hnm2 = Hnm1;
        Hnm1 = Hn;
    }
    return Hn;
}

#endif
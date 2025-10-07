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

#endif
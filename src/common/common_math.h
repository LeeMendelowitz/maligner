#ifndef MALIGNER_MATH_H
#define MALIGNER_MATH_H
#include <vector>

double mad(const std::vector<double>& in);
double mad(const std::vector<double>& in, double m); // For when median m already known.
double median(const std::vector<double>&in);
double median_sorted(const std::vector<double>&in);

#endif
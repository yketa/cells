#include "tools.hpp"

#include <cmath>
#include <numbers>
#include <vector>

double cross2(std::vector<double> const& a, std::vector<double> const& b)
    { return a[0]*b[1] - a[1]*b[0]; }

std::vector<double> cross2z(std::vector<double> const& a)
    { return {a[1], -a[0]}; }

double hexagonEdgeLength(double const& hexagonArea)
    { return std::sqrt((2*hexagonArea)/(3*tan(std::numbers::pi/3.))); }


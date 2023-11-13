#include "tools.hpp"

#include <vector>                                                               

double const cross2(
    std::vector<double> const& a, std::vector<double> const& b)
    { return a[0]*b[1] - a[1]*b[0]; }                                           
                                                                                
std::vector<double> const cross2z(std::vector<double> const& a)
    { return {a[1], -a[0]}; }                                                   


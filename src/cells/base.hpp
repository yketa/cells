/*
Base types for base types.
*/

#ifndef BASE_HPP
#define BASE_HPP

#include <map>
#include <vector>

typedef std::map<std::string, double> ParametersType;
typedef std::map<long int, std::vector<double>> ForcesType;
typedef ForcesType VelocitiesType;

#endif


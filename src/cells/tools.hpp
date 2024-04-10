#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

/*
 *  TEMPLATED FUNCTIONS
 *
 */

template<class T> int const sign(T const& val)                                  // sign function (https://stackoverflow.com/questions/1903954)
    { return (T(0) < val) - (val < T(0)); }

template<class T, class TT> T const pmod(T const& i, TT const& n)               // positive modulo (remainder) i%n (https://stackoverflow.com/questions/14997165)
    { return std::fmod(std::fmod(i, n) + n, n); }
template<class T> long int const qpmod(T const& i, T const& n)                  // quotient `i/n' associated to positive remainder `i%n = pmod(i, n)'
    { return (i - pmod<T>(i, n))/n; }

template<class T> void eraseInVec(std::vector<T>& vec, T const& val) {          // remove ONE occurence of value in vector (https://stackoverflow.com/questions/3385229)
    typename std::vector<T>::iterator position =
        std::find(vec.begin(), vec.end(), val);
    if (position != vec.end()) { vec.erase(position); } // if element was found
}

template<class T> T norm2(std::vector<T> vec)                                   // norm 2 of vector
    { T nsq = 0; for (T el : vec) nsq += el*el; return sqrt(nsq); }

template<class T>
bool const inVec(std::vector<T> const& vec, T const& val)                       // is value in vector
    { return (std::find(vec.begin(), vec.end(), val) != vec.end()); }
template<class T, class TT>
bool const inMap(std::map<T, TT> const& map, T const& key)                      // is key in map
    { return (map.find(key) != map.end()); }

template<class T, class TT> T const maxKey(std::map<T, TT> const& map)          // return max key of map (https://stackoverflow.com/questions/1660195)
    { return map.rbegin()->first; }
template<class T, class TT> void printMap(std::map<T, TT> const& map) {         // print content of map
    for (auto it=map.begin(); it != map.end(); ++it) {
        std::cout << it->first << ": " << it->second << std::endl;
    }
}

template<class T> double const angle2(T const& x, T const& y)                   // angle of 2D vector (x, y) with respect to horizontal axis
    { double angle = acos(x/sqrt(x*x + y*y)); return (y > 0 ? 1 : - 1)*angle; }
template<class T> double const angle2(std::vector<T> const& vec)
    { return angle2(vec[0], vec[1]); }
template<class T> double const angle2(T* const& array)
    { return angle2(array[0], array[1]); }

/*
 *  FUNCTION PROTOTYPES
 *
 */

double cross2(std::vector<double> const& a, std::vector<double> const& b);  // cross product of 2D vector

std::vector<double> cross2z(std::vector<double> const& a);                  // cross product of 2D (x, y[, 0]) vector with (0, 0, 1) projected in the xy-plane

/*
 *  CLASSES
 *
 */

class Counter {
/*
Return integers in order at each call.
*/

    private:

        long int counter;

    public:

        Counter(long int const& initial=0) : counter(initial) {}    // initialise counter

        Counter& operator=(long int const& current) // lvalue assignment: set counter
            { counter = current; return *this; }
        operator long int() const                   // rvalue assignment: return counter
            { return counter; }
        long int operator()()                       // increment and return counter
            { counter++; return counter - 1; }

};

#endif


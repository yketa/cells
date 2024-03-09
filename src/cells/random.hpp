/*
Random numbers generator.
*/

#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>
#include <algorithm>
#include <ostream>

#define RNDG std::mt19937 // Mersenne Twister

class Random {
/*
This simple class is a wrapper for the built-in C++ random number generator.
Adapted from code initially written by Robert L. Jack.
*/

    private:

        RNDG generator;
        // std::default_random_engine generator; // (old) period = 1,686,629,820

        int const _max = 0x7fffffff;
        std::uniform_int_distribution<int>* const intmax
            = new std::uniform_int_distribution<int>(0, _max);
        std::uniform_real_distribution<double>* const real01
            = new std::uniform_real_distribution<double>(0.0, 1.0);
        std::normal_distribution<double>* const normal
            = new std::normal_distribution<double>(0.0, 1.0);

    public:

        // CONSTRUCTORS AND DESTRUCTORS

        Random(long int const& seed=0) : generator() { generator.seed(seed); }
        Random(RNDG const& rndeng) : generator(rndeng) {}
        Random(Random const& random_) : Random(random_.getGenerator()) {}   // copy constructor

        ~Random() { delete intmax; delete real01; delete normal; }

        // GETTERS AND SETTERS

        RNDG const& getGenerator() const { return generator; }
        void setGenerator(RNDG const& rndeng) { generator = rndeng; }
        // std::normal_distribution<double>* const getNormal() { return normal; }

        // METHODS

        double const random01()                                     // get random double in [0, 1] w/ uniform distribution
            { return (*real01)(generator); }
        int const randomInt(int const& max)                         // get random integer in [|0, max - 1|] w/ uniform distribution
            { return (*intmax)(generator) % max; }
        double const gauss()                                        // get random real number w/ normal distribution N(0, 1)
            { return (*normal)(generator); }
        double const gauss(double const& mean, double const& std)   // get random real number w/ normal distribution N(mean, std)
            { return std::normal_distribution<double>(mean, std)(generator); }

        template<class RandomAccessIterator>                        // randomly rearrange elements (https://cplusplus.com/reference/algorithm/shuffle/)
            void shuffle(RandomAccessIterator first, RandomAccessIterator last) 
            { std::shuffle(first, last, generator); }
        template<class T> void shuffle(std::vector<T>& vec)         // randomly rearrange elements of std::vector
            { shuffle(vec.begin(), vec.end()); }
        template<class T> T const pick(std::vector<T> const& vec)   // randomly pick element of std::vector
            { return vec[randomInt(vec.size())]; }

        // OPERATORS

        friend std::ostream& operator<<(std::ostream& os, Random const& rnd)    // overload << operator
            { os << rnd.generator; return os; }
        friend std::istream& operator>>(std::istream& is, Random& rnd)          // overload >> operator
            { is >> rnd.generator; return is; }

        void operator=(Random const& rnd) { setGenerator(rnd.getGenerator()); } // copy operator

};

#endif


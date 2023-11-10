#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <map>
#include <vector>
#include <iostream>
#include <assert.h>
#include <string>
#include <algorithm>
#include <math.h>

/*
 * FUNCTIONS
 *
 */

template<class T> int const sign(T const& val)                                  // sign function (https://stackoverflow.com/questions/1903954)
    { return (T(0) < val) - (val < T(0)); }

template<class T> T const pmod(T const& i, T const& n) { return (i%n + n)%n; }  // positive modulo i%n (https://stackoverflow.com/questions/14997165)

double const cross2(std::vector<double> const& a, std::vector<double> const& b) // cross product of 2D vector
    { return a[0]*b[1] - a[1]*b[0]; }

std::vector<double> const cross2z(std::vector<double> const& a)                 // cross product of 2D (x, y[, 0]) vector with (0, 0, 1) projected in the xy-plane
    { return {a[1], -a[0]}; }

template<class T> void eraseInVec(std::vector<T> vec, T const& val) {           // remove ONE occurence of value in vector (https://stackoverflow.com/questions/3385229)
    typename std::vector<T>::iterator position =
        std::find(vec.begin(), vec.end(), val);
    if (position != vec.end()) { vec.erase(position); } // if element was found
}

template<class T> bool const inVec(std::vector<T> const& vec, T const& val)     // is value in vector
    { return std::find(vec.begin(), vec.end(), val) == vec.end(); }

template<class T> double const angle2(T x, T y)                                 // angle of 2D vector (x, y) with respect to horizontal axis
    { double angle = acos(x/sqrt(x*x + y*y)); return (y > 0 ? 1 : - 1)*angle; }
template<class T> double const angle2(std::vector<T> vec)
    { return angle2(vec[0], vec[1]); }
template<class T> double const angle2(T* array)
    { return angle2(array[0], array[1]); }

/*
 * CLASSES
 *
 */

class Counter;
template<class T> class MultiIntKeyDict;

class Counter {
/*
Return integers in order at each call.
*/

    private:

        long int counter;

    public:

        Counter(long int const& initial=0) : counter(initial - 1) {}    // initialise counter

        Counter& operator=(long int const& current) // set counter
            { counter = current; return *this; }
        long int operator()()                       // increment and return counter
            { counter++; return counter; }

};

template<class T> class MultiIntKeyDict {
/*
Dictionary-like (map-like) object in which multiple integer keys are associated
to the same value.
*/

    private:

        std::map<long int, long int> keys;  // secondary dictionary which associates multiple keys to an unique key of data
        std::map<long int, T> data;         // dictionary values with unique identification
        Counter maxIndex;                   // largest value of unique indices

        /*
        uses proxy classes to set values
        https://stackoverflow.com/questions/20168543
        https://stackoverflow.com/questions/994488
        https://stackoverflow.com/questions/39926941
        */

        class Proxy {
            private:
                MultiIntKeyDict* mikd;
                long int const key;
            public:
                Proxy(MultiIntKeyDict& mikd_, long int const& key_)
                    : mikd(&mikd_), key(key_) {}
                void operator=(T const& value) {    // lvalue assignment: set value
                    mikd->erase(key);                           // delete entry corresponding to key
                    (mikd->keys)[key] = mikd->maxIndex();       // create new index for key
                    (mikd->data)[(mikd->keys)[key]] = value;    // set value
                }
                operator T() const {                // rvalue assignment: return value
                    assert(mikd->in(key));                      // check that key is in dictionary
                    return (mikd->data)[(mikd->keys)[key]];
                }
        };

        class ProxyBis {
            private:
                MultiIntKeyDict* mikd;
                std::vector<long int> keys;
            public:
                ProxyBis(MultiIntKeyDict& mikd_,
                    std::vector<long int> const& keys_)
                    : mikd(&mikd_), keys(keys_) {}
                void operator=(T const& value) {    // lvalue assignment: set values
                    if (keys.size() >= 1) {
                        Proxy(*mikd, keys.at(0)) = value;   // set first one
                        for (long int i=1; i < keys.size(); i++) {
                            mikd->erase(keys.at(i));        // delete others (if existing)
                            (mikd->keys)[keys.at(i)] =      // then set identical to first one
                                (mikd->keys)[keys.at(0)];
                        }
                    }
                }
        };

    public:

        MultiIntKeyDict() {}

        // SET

        Proxy operator[] (long int const& key)                  // set and return value from index
            { return Proxy(*this, key); }
        ProxyBis operator[] (std::vector<long int> const& keys) // set from vector or initialiser list of indices
            { return ProxyBis(*this, keys); }

        bool in(long int const& key) { return (keys.find(key) != keys.end()); } // is key in the dictionary?

        // REMOVE

        void erase(long int const& key) {               // remove entry
            if (in(key)) { // key is in dictionary
                long int index = keys[key];
                keys.erase(key);
                // https://stackoverflow.com/questions/4263640
                for (auto it=keys.begin(); it != keys.end(); ++it) {
                    if (it->second == index) { return; }    // there is at least one other key with same value then quit
                }
                data.erase(index);                          // there is no other key for the same value then delete value
            }
        }
        void erase(std::vector<long int> const& keys) { // remove entries
            for (long int key : keys) { erase(key); }
        }
        void clear() {                                  // remove all entries
            keys.clear();
            data.clear();
            maxIndex = -1;
        }

        // CONSULT

        std::vector<T*> getValues() {
            std::vector<T*> values(0);
            for (auto it=data.begin(); it != data.end(); ++it) {
                values.push_back(&(it->second));
            }
            return values;
        }

        // OUTPUT

        void print(std::string const& comment="") {
            std::cout << "[KEYS] " << comment << std::endl;
            for (auto it=keys.begin(); it != keys.end(); it++) {
                std::cout << it->first << ": " << it->second << std::endl;
            }
            std::cout << "[DATA] " << comment << std::endl;
            for (auto it=data.begin(); it != data.end(); it++) {
                std::cout << it->first << ": " << it->second << std::endl;
            }
            std::cout << std::endl;
        }
};

#endif


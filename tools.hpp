#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <map>
#include <vector>
#include <iostream>
#include <assert.h>
#include <string>
#include <algorithm>

/*
 * FUNCTIONS
 *
 */

template<class T> int sign(T const& val)                                        // sign function (https://stackoverflow.com/questions/1903954)
    { return (T(0) < val) - (val < T(0)); }

double const cross2(std::vector<double> const& a, std::vector<double> const& b) // cross product of 2D vector
    { return a[0]*b[1] - a[1]*b[0]; }

template<class T> void eraseInVec(std::vector<T> vec, T const& val) {           // remove ONE occurence of value in vector (https://stackoverflow.com/questions/3385229)
    typename std::vector<T>::iterator position =
        std::find(vec.begin(), vec.end(), val);
    if (position != vec.end()) { vec.erase(position); } // if element was found
}

template<class T> bool const inVec(std::vector<T> const& vec, T const& val)     // is value in vector
    { return std::find(vec.begin(), vec.end(), val) == vec.end(); }

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
                long int key;
            public:
                Proxy(MultiIntKeyDict& mikd_, long int& key_)
                    : mikd(&mikd_), key(key_) {}
                void operator=(T value) {   // lvalue assignment: set value
                    mikd->erase(key);                           // delete entry corresponding to key
                    (mikd->keys)[key] = mikd->maxIndex();       // create new index for key
                    (mikd->data)[(mikd->keys)[key]] = value;    // set value
                }
                operator T() const {        // rvalue assignment: return value
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
                    std::vector<long int>& keys_)
                    : mikd(&mikd_), keys(keys_) {}
                void operator=(T value) {   // lvalue assignment: set values
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

        Proxy operator[] (long int key)                     // set and return value from index
            { return Proxy(*this, key); }
        ProxyBis operator[] (std::vector<long int> keys)    // set from vector or initialiser list of indices
            { return ProxyBis(*this, keys); }

        bool in(long int const& key) { return (keys.find(key) != keys.end()); } // is key in the dictionary?

        // REMOVE

        void erase(long int key) {                  // remove entry
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
        void erase(std::vector<long int> keys) {    // remove entries
            for (long int key : keys) { erase(key); }
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

        void print(std::string comment="") {
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


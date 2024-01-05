/*
A simple implementation of a class factory.
*/

#ifndef CLASS_FACTORY_HPP
#define CLASS_FACTORY_HPP

#include <map>
#include <string>
#include <memory>

template<class BaseType> class ClassFactory {
/*
The assumption is that one has a class hierarchy with a base class (BaseType)
and its children (DerivedType). Class factory then creates a mapping (map)
between a string (class name) and a smart pointer to an object of a given base
or derived class.
Adapted from code initially written by Rastko Sknepnek.
*/

    protected:

        std::map<std::string, std::shared_ptr<BaseType>> factory_map;   // map of strings to objects in the factory

        /*                                                                      
        ClassFactory objects contain std::shared_ptr which can be copied. Would 
        they have been std::unique_ptr these could not be copied and thus e.g.  
        not be passed as arguments.                                             
        https://stackoverflow.com/questions/41060167                            
        https://stackoverflow.com/questions/6876751                             
        */

    public:

        // CONSTRUCTOR AND DESTRUCTOR

        ClassFactory() = default;
        ~ClassFactory() = default;

        // CONSULT

        bool in(std::string const& key) const
        /*
        Indicate if key is in the factory.
        */
            { return (factory_map.find(key) != factory_map.end()); }

        std::shared_ptr<BaseType>& get(std::string const& key)
        /*
        Return pointer to the object associated with the key.
        */
            { return factory_map[key]; }

        std::map<std::string, std::shared_ptr<BaseType>> const& map() const
        /*
        Return the entire factory map.
        */
            { return factory_map; }

        std::map<std::string, std::shared_ptr<BaseType>>::iterator
            begin() { return factory_map.begin(); }
        std::map<std::string, std::shared_ptr<BaseType>>::iterator
            end() { return factory_map.end(); }
        std::map<std::string, std::shared_ptr<BaseType>>::const_iterator
            const_begin() const { return factory_map.begin(); }
        std::map<std::string, std::shared_ptr<BaseType>>::const_iterator
            const_end() const { return factory_map.end(); }

        // SET

        template<class DerivedType, typename... Args> void add(
            std::string const& key, Args... args)
        /*
        Add new object to the factory.
        */
            { factory_map[key] = std::make_shared<DerivedType>(args...); }

        // REMOVE

        void remove(std::string const& key)
        /*
        Remove object from the factory.
        */
            { factory_map.erase(key); }

        void clear()
        /*
        Remove all entries.
        */
            { factory_map.clear(); }

};

#endif


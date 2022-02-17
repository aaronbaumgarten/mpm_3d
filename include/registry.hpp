//
// Created by aaron on 4/24/18.
// registry.hpp
//

#ifndef MPM_V3_REGISTRY_HPP
#define MPM_V3_REGISTRY_HPP

#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <map>
#include <memory>

/*----------------------------------------------------------------------------*/
//Registry object templated on a base class which stores derived classes by string
//the purpose of this object is to allow runtime determination of which
//derived class to instantiate.
/*----------------------------------------------------------------------------*/

template<typename T>
class Registry {
public:
    Registry();

    std::map<std::string,std::unique_ptr<T>(*)()> object;

    std::unique_ptr<T> get_object(std::string strIn) {
        auto search = object.find(strIn);
        if(search == object.end()) {
            std::cerr << "Undefined registry string: \"" << strIn << "\". Are you sure this class has been added to the registry in registry.cpp?" << std::endl;
            throw std::invalid_argument("String \"" + strIn + "\" not found in registry.");
        } else {
            return search->second();
        }
    }
};

#endif //MPM_V3_REGISTRY_HPP

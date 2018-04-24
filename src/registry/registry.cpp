//
// Created by aaron on 4/24/18.
// registry.cpp
//

#include "registry.hpp"

template<typename Base, typename Derived> std::unique_ptr<B> createInstance() { return std::unique_ptr<Base>(new Derived); }

/*template<>
Registry<Base>::Registry() {
    object["DerivedA"] = &createInstance<Base,DerivedA>;
    object["DerivedB"] = &createInstance<Base,DerivedB>;
}*/
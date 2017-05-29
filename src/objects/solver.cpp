//
// Created by aaron on 5/13/17.
// solver.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <dlfcn.h>

#include "stringparser.hpp"
#include "job.hpp"

#include "serializer.hpp"
#include "solver.hpp"

Solver::Solver() {
    filename = "";
    filepath = "";
    fp64_props = std::vector<double>();
    int_props = std::vector<int>();
    str_props = std::vector<std::string>();

    handle = NULL;

    solverInit = NULL;
    solverStep = NULL;

    solverSaveState = NULL;
    solverLoadState = NULL;
}

Solver::~Solver() {
    if (handle){
        dlclose(handle);
    }
}

void Solver::solverSetPlugin(Job* job, std::string nameIN, std::string pathIN, std::vector<double> fp64IN, std::vector<int> intIN, std::vector<std::string> strIN){
    filename = nameIN;
    filepath = StringParser::stringMakeDirectory(pathIN);
    fp64_props = fp64IN;
    int_props = intIN;
    str_props = strIN;

    handle = dlopen((filepath+filename).c_str(), RTLD_LAZY);

    solverSetFnPointers(handle);

    return;
}

void Solver::solverSetFnPointers(void* handle){
    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';

        solverInit = NULL;
        solverStep = NULL;

        solverSaveState = NULL;
        solverLoadState = NULL;

    } else {
        dlerror();
        solverInit = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "solverInit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'solverInit': " << dlsym_error <<
                      '\n';
        }

        solverStep = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "solverStep"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'solverStep': " << dlsym_error <<
                      '\n';
        }

        solverSaveState = reinterpret_cast<std::string (*)(Job*,Serializer*,std::string)>(dlsym(handle, "solverSaveState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'solverSaveState': " << dlsym_error <<
                      '\n';
        }
        solverLoadState = reinterpret_cast<int (*)(Job*,Serializer*,std::string)>(dlsym(handle, "solverLoadState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'solverLoadState': " << dlsym_error <<
                      '\n';
        }
    }
    return;
}
//
// Created by aaron on 5/13/17.
// driver.cpp
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
#include "driver.hpp"

Driver::Driver() {
    fullpath = "";
    filename = "";
    filepath = "";
    fp64_props = std::vector<double>();
    int_props = std::vector<int>();
    str_props = std::vector<std::string>();

    handle = NULL;

    driverInit = NULL;
    driverRun = NULL;
    driverGenerateGravity = NULL;
    driverApplyGravity = NULL;

    driverSaveState = NULL;
    driverLoadState = NULL;
}

Driver::~Driver() {
    if (handle){
        dlclose(handle);
    }
}

void Driver::driverSetPlugin(Job* job, std::string pathIN, std::string nameIN, std::vector<double> fp64IN, std::vector<int> intIN, std::vector<std::string> strIN){
    filename = nameIN;
    fullpath = StringParser::stringMakeDirectory(pathIN);
    fp64_props = fp64IN;
    int_props = intIN;
    str_props = strIN;

    handle = dlopen((fullpath+filename).c_str(), RTLD_LAZY);

    driverSetFnPointers(handle);

    return;
}

void Driver::driverSetFnPointers(void* handle){
    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';

        driverInit = NULL;
        driverRun = NULL;
        driverGenerateGravity = NULL;
        driverApplyGravity = NULL;

        driverSaveState = NULL;
        driverLoadState = NULL;
    } else {
        dlerror();
        driverInit = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "driverInit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'driverInit': " << dlsym_error <<
                      '\n';
        }

        driverRun = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "driverRun"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'driverRun': " << dlsym_error <<
                      '\n';
        }

        driverGenerateGravity = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "driverGenerateGravity"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'driverGenerateGravity': " << dlsym_error <<
                      '\n';
        }

        driverApplyGravity = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "driverApplyGravity"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'driverApplyGravity': " << dlsym_error <<
                      '\n';
        }

        driverSaveState = reinterpret_cast<std::string (*)(Job*,Serializer*,std::string)>(dlsym(handle, "driverSaveState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'driverSaveState': " << dlsym_error <<
                      '\n';
        }
        driverLoadState = reinterpret_cast<int (*)(Job*,Serializer*,std::string)>(dlsym(handle, "driverLoadState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'driverLoadState': " << dlsym_error <<
                      '\n';
        }
    }
    return;
}
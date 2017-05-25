//
// Created by aaron on 5/19/17.
// boundary.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <dlfcn.h>

#include "job.hpp"

#include "serializer.hpp"

#include "body.hpp"
#include "boundary.hpp"

Boundary::Boundary() {
    filename = "";
    filepath = "";
    fp64_props = std::vector<double>();
    int_props = std::vector<int>();
    str_props = std::vector<std::string>();

    handle = NULL;

    boundaryInit = NULL;
    boundaryGenerateRules = NULL;
    boundaryApplyRules = NULL;

    boundaryWriteFrame = NULL;
    boundarySaveState = NULL;
    boundaryLoadState = NULL;
}

Boundary::~Boundary() {
    if (handle){
        dlclose(handle);
    }
}

void Boundary::boundarySetPlugin(Job* job, Body* body, std::string nameIN, std::string pathIN, std::vector<double> fp64IN, std::vector<int> intIN, std::vector<std::string> strIN){
    filename = nameIN;
    filepath = pathIN;
    fp64_props = fp64IN;
    int_props = intIN;
    str_props = strIN;

    handle = dlopen((filepath+filename).c_str(), RTLD_LAZY);

    boundarySetFnPointers(handle);

    return;
}

void Boundary::boundarySetFnPointers(void* handle){
    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';

        boundaryInit = NULL;
        boundaryGenerateRules = NULL;
        boundaryApplyRules = NULL;

        boundaryWriteFrame = NULL;
        boundarySaveState = NULL;
        boundaryLoadState = NULL;

    } else {
        dlerror();
        boundaryInit = reinterpret_cast<void (*)(Job *, Body*)>(dlsym(handle, "boundaryInit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'boundaryInit': " << dlsym_error <<
                      '\n';
        }
        boundaryGenerateRules = reinterpret_cast<void (*)(Job *, Body *)>(dlsym(handle, "boundaryGenerateRules"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'boundaryGenerateRules': " << dlsym_error <<
                      '\n';
        }
        boundaryApplyRules = reinterpret_cast<void (*)(Job *, Body *)>(dlsym(handle, "boundaryApplyRules"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'boundaryApplyRules': " << dlsym_error <<
                      '\n';
        }

        boundaryWriteFrame = reinterpret_cast<void (*)(Job*,Body*,Serializer*)>(dlsym(handle, "boundaryWriteFrame"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'boundaryWriteFrame': " << dlsym_error <<
                      '\n';
        }
        boundarySaveState = reinterpret_cast<std::string (*)(Job*,Body*,Serializer*,std::string)>(dlsym(handle, "boundarySaveState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'boundarySaveState': " << dlsym_error <<
                      '\n';
        }
        boundaryLoadState = reinterpret_cast<int (*)(Job*,Body*,Serializer*,std::string)>(dlsym(handle, "boundaryLoadState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'boundaryLoadState': " << dlsym_error <<
                      '\n';
        }
    }
    return;
}
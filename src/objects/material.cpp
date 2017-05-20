//
// Created by aaron on 5/19/17.
// material.cpp
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
#include "material.hpp"

Material::Material() {
    filename = "";
    filepath = "";
    fp64_props = std::vector<double>();
    int_props = std::vector<int>();

    handle = NULL;

    materialInit = NULL;
    materialCalculateStress = NULL;
    materialAssignStress = NULL;
    materialAssignPressure = NULL;

    materialWriteFrame = NULL;
    materialSaveState = NULL;
    materialLoadState = NULL;
}

Material::~Material() {
    if (handle){
        dlclose(handle);
    }
}

void Material::materialSetPlugin(Job* job, Body* body, std::string nameIN, std::string pathIN, std::vector<double> fp64IN, std::vector<int> intIN){
    filename = nameIN;
    filepath = pathIN;
    fp64_props = fp64IN;
    int_props = intIN;

    handle = dlopen((filepath+filename).c_str(), RTLD_LAZY);

    materialSetFnPointers(handle);

    return;
}

void Material::materialSetFnPointers(void* handle){
    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';

        materialInit = NULL;
        materialCalculateStress = NULL;
        materialAssignStress = NULL;
        materialAssignPressure = NULL;

        materialWriteFrame = NULL;
        materialSaveState = NULL;
        materialLoadState = NULL;

    } else {
        dlerror();
        materialInit = reinterpret_cast<void (*)(Job *, Body*)>(dlsym(handle, "materialInit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'materialInit': " << dlsym_error <<
                      '\n';
        }
        materialCalculateStress = reinterpret_cast<void (*)(Job *, Body *, int=1)>(dlsym(handle, "materialCalculateStress"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'materialCalculateStress': " << dlsym_error <<
                      '\n';
        }
        materialAssignStress = reinterpret_cast<void (*)(Job *, Body *, Eigen::MatrixXd, int)>(dlsym(handle, "materialAssignStress"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'materialAssignStress': " << dlsym_error <<
                      '\n';
        }
        materialAssignPressure = reinterpret_cast<void (*)(Job *, Body *, double, int)>(dlsym(handle, "materialAssignPressure"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'materialAssignPressure': " << dlsym_error <<
                      '\n';
        }

        materialWriteFrame = reinterpret_cast<void (*)(Job*,Body*,Serializer*)>(dlsym(handle, "materialWriteFrame"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'materialWriteFrame': " << dlsym_error <<
                      '\n';
        }
        materialSaveState = reinterpret_cast<std::string (*)(Job*,Body*,Serializer*,std::string)>(dlsym(handle, "materialSaveState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'materialSaveState': " << dlsym_error <<
                      '\n';
        }
        materialLoadState = reinterpret_cast<int (*)(Job*,Body*,Serializer*,std::string)>(dlsym(handle, "materialLoadState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'materialLoadState': " << dlsym_error <<
                      '\n';
        }
    }
    return;
}
//
// Created by aaron on 5/12/17.
// serializer.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <dlfcn.h>

#include "job.hpp"

#include "serializer.hpp"

Serializer::Serializer() {
    filename = "";
    filepath = "";
    fp64_props = std::vector<double>();
    int_props = std::vector<int>();

    handle = NULL;

    serializerInit = NULL;

    serializerWriteFrame = NULL;
    serializerWriteScalar = NULL;
    serializerWriteVector = NULL;
    serializerWriteTensor = NULL;

    serializerSaveState = NULL;
    serializerLoadState = NULL;
}

Serializer::~Serializer() {
    if (handle){
        dlclose(handle);
    }
}

void Serializer::serializerSetPlugin(Job* job, std::string nameIN, std::string pathIN, std::vector<double> fp64IN, std::vector<int> intIN){
    filename = nameIN;
    filepath = pathIN;
    fp64_props = fp64IN;
    int_props = intIN;

    handle = dlopen((filepath+filename).c_str(), RTLD_LAZY);

    serializerSetFnPointers(handle);

    return;
}

void Serializer::serializerSetFnPointers(void* handle){
    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';

        serializerInit = NULL;

        serializerWriteFrame = NULL;
        serializerWriteScalar = NULL;
        serializerWriteVector = NULL;
        serializerWriteTensor = NULL;

        serializerSaveState = NULL;
        serializerLoadState = NULL;
    } else {
        dlerror();
        serializerInit = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "serializerInit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerInit': " << dlsym_error <<
                      '\n';
        }

        serializerWriteFrame = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "serializerWriteFrame"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerWriteFrame': " << dlsym_error <<
                      '\n';
        }
        serializerWriteScalar = reinterpret_cast<void (*)(Eigen::Matrix*,std::string)>(dlsym(handle, "serializerWriteScalar"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerWriteScalar': " << dlsym_error <<
                      '\n';
        }
        serializerWriteVector = reinterpret_cast<void (*)(Eigen::Matrix*,std::string)>(dlsym(handle, "serializerWriteVector"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerWriteVector': " << dlsym_error <<
                      '\n';
        }
        serializerWriteTensor = reinterpret_cast<void (*)(Eigen::Matrix*, std::string)>(dlsym(handle, "serializerWriteTensor"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerWriteTensor': " << dlsym_error <<
                      '\n';
        }

        serializerSaveState = reinterpret_cast<std::string (*)(Job*)>(dlsym(handle, "serializerSaveState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerSaveState': " << dlsym_error <<
                      '\n';
        }
        serializerLoadState = reinterpret_cast<int (*)(Job*,std::string)>(dlsym(handle, "serializerLoadState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerLoadState': " << dlsym_error <<
                      '\n';
        }
    }
    return;
}
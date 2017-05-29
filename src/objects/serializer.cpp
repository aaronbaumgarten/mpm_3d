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

#include "stringparser.hpp"
#include "job.hpp"

#include "serializer.hpp"

Serializer::Serializer() {
    filename = "";
    filepath = "";
    fp64_props = std::vector<double>();
    int_props = std::vector<int>();
    str_props = std::vector<std::string>();

    handle = NULL;

    serializerInit = NULL;

    serializerWriteFrame = NULL;
    serializerWriteScalarArray = NULL;
    serializerWriteVectorArray = NULL;
    serializerWriteTensorArray = NULL;

    serializerSaveState = NULL;
    serializerLoadState = NULL;
}

Serializer::~Serializer() {
    if (handle){
        dlclose(handle);
    }
}

void Serializer::serializerSetPlugin(Job* job, std::string nameIN, std::string pathIN, std::vector<double> fp64IN, std::vector<int> intIN, std::vector<std::string> strIN){
    filename = nameIN;
    filepath = StringParser::stringMakeDirectory(pathIN);
    fp64_props = fp64IN;
    int_props = intIN;
    str_props = strIN;

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
        serializerWriteScalarArray = NULL;
        serializerWriteVectorArray = NULL;
        serializerWriteTensorArray = NULL;

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

        serializerWriteFrame = reinterpret_cast<int (*)(Job *)>(dlsym(handle, "serializerWriteFrame"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerWriteFrame': " << dlsym_error <<
                      '\n';
        }
        serializerWriteScalarArray = reinterpret_cast<void (*)(Eigen::VectorXd&,std::string)>(dlsym(handle, "serializerWriteScalarArray"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerWriteScalarArray': " << dlsym_error <<
                      '\n';
        }
        serializerWriteVectorArray = reinterpret_cast<void (*)(Eigen::MatrixXd&,std::string)>(dlsym(handle, "serializerWriteVectorArray"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerWriteVectorArray': " << dlsym_error <<
                      '\n';
        }
        serializerWriteTensorArray = reinterpret_cast<void (*)(Eigen::MatrixXd&, std::string)>(dlsym(handle, "serializerWriteTensorArray"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'serializerWriteTensorArray': " << dlsym_error <<
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

/*----------------------------------------------------------------------------*/

void Serializer::serializerSetMainPath(Job* job, std::string program){
    std::vector<std::string> svec;
    svec = StringParser::stringSplitString(program,'/');
    std::string filepath = "";
    for (size_t i=0; i<(svec.size()-1);i++){
        filepath += svec[i];
        filepath += "/";
    }
    mainpath = filepath;
    return;
}
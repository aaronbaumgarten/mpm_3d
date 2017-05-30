//
// Created by aaron on 5/14/17.
// contact.cpp
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
#include "contact.hpp"

Contact::Contact() {
    id = 0;
    name = "";
    filename = "";
    filepath = "";
    fp64_props = std::vector<double>();
    int_props = std::vector<int>();
    str_props = std::vector<std::string>();

    handle = NULL;

    contactInit = NULL;
    contactGenerateRules = NULL;
    contactApplyRules = NULL;

    contactWriteFrame = NULL;
    contactSaveState = NULL;
    contactLoadState = NULL;
}

Contact::~Contact() {
    if (handle){
        dlclose(handle);
    }
}

void Contact::contactSetPlugin(Job* job, std::string pathIN, std::string nameIN, std::vector<double> fp64IN, std::vector<int> intIN, std::vector<std::string> strIN){
    filename = nameIN;
    filepath = StringParser::stringMakeDirectory(pathIN);
    fp64_props = fp64IN;
    int_props = intIN;
    str_props = strIN;

    handle = dlopen((filepath+filename).c_str(), RTLD_LAZY);

    contactSetFnPointers(handle);

    return;
}

void Contact::contactSetFnPointers(void* handle){
    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';

        contactInit = NULL;
        contactGenerateRules = NULL;
        contactApplyRules = NULL;

        contactWriteFrame = NULL;
        contactSaveState = NULL;
        contactLoadState = NULL;

    } else {
        dlerror();
        contactInit = reinterpret_cast<void (*)(Job *, Contact*)>(dlsym(handle, "contactInit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contactInit': " << dlsym_error <<
                      '\n';
        }
        contactGenerateRules = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "contactGenerateRules"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contactGenerateRules': " << dlsym_error <<
                      '\n';
        }
        contactApplyRules = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "contactApplyRules"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contactApplyRules': " << dlsym_error <<
                      '\n';
        }

        contactWriteFrame = reinterpret_cast<void (*)(Job*,Serializer*)>(dlsym(handle, "contactWriteFrame"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contactWriteFrame': " << dlsym_error <<
                      '\n';
        }
        contactSaveState = reinterpret_cast<std::string (*)(Job*,Serializer*,std::string)>(dlsym(handle, "contactSaveState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contactSaveState': " << dlsym_error <<
                      '\n';
        }
        contactLoadState = reinterpret_cast<int (*)(Job*,Serializer*,std::string)>(dlsym(handle, "contactLoadState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contactLoadState': " << dlsym_error <<
                      '\n';
        }
    }
    return;
}
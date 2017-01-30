//
// Created by aaron on 12/23/16.
// contact.cpp
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include "contact.hpp"
#include "process.hpp"
#include <dlfcn.h>

Contact::Contact(){
    id = 0;
    bodyIDs = std::vector<int>();
    num_fp64_props = 0;
    num_int_props = 0;
    fp64_props = std::vector<double>();
    int_props = std::vector<int>();
    contact_filename = "";
    handle = NULL;

    contact_init = NULL;
    resolve_contact = NULL;
}

Contact::Contact(std::string filename, std::string filepath, size_t idIn, std::vector<int> bodyIDsIn, std::vector<double> fp64props, std::vector<int> intprops){
    id = idIn;
    bodyIDs = bodyIDsIn;
    //std::string filepath = "src/contacts/";
    filepath += filename;

    use_builtin = 0;
    handle = dlopen(filepath.c_str(), RTLD_LAZY);
    contact_filename = filename;

    num_fp64_props = fp64props.size();
    num_int_props = intprops.size();
    fp64_props = fp64props;
    int_props = intprops;

    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';
        contact_init = NULL;
        resolve_contact = NULL;
    } else {
        dlerror();
        contact_init = reinterpret_cast<void (*)(job_t *, size_t)>(dlsym(handle, "contact_init"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contact_init': " << dlsym_error <<
            '\n';
        }
        resolve_contact = reinterpret_cast<void (*)(job_t *, size_t)>(dlsym(handle, "resolve_contact"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'resolve_contact': " << dlsym_error <<
            '\n';
        }
    }
}

Contact::~Contact() {
    //std::cout << "closing material" << std::endl;
    if(handle){
        dlclose(handle);
    }
}

void Contact::setContact(std::string filename, std::string filepath, size_t idIn, std::vector<int> bodyIDsIn, std::vector<double> fp64props, std::vector<int> intprops){
    this->id = idIn;
    this->bodyIDs = bodyIDsIn;
    //std::string filepath = "src/contacts/";
    filepath += filename;

    this->use_builtin = 0;
    this->handle = dlopen(filepath.c_str(), RTLD_LAZY);
    this->contact_filename = filename;

    this->num_fp64_props = fp64props.size();
    this->num_int_props = intprops.size();
    this->fp64_props = fp64props;
    this->int_props = intprops;

    char* dlsym_error;
    if (!this->handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';
        this->contact_init = NULL;
        this->resolve_contact = NULL;
    } else {
        dlerror();
        this->contact_init = reinterpret_cast<void (*)(job_t *, size_t)>(dlsym(this->handle, "contact_init"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contact_init': " << dlsym_error <<
            '\n';
        }
        this->resolve_contact = reinterpret_cast<void (*)(job_t *, size_t)>(dlsym(this->handle, "resolve_contact"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'resolve_contact': " << dlsym_error <<
            '\n';
        }
    }
}
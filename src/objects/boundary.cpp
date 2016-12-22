//
// Created by aaron on 10/26/16.
// boundary.cpp
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "boundary.hpp"
#include "process.hpp"
#include <dlfcn.h>

//namespace boundary {
//#include "floorBC.cpp"
//#include "floorBC2D.cpp"
//#include "boxBC2d.cpp"
//#include "boxBC.cpp"
//}

Boundary::Boundary(){
    use_builtin = 0;
    handle = NULL;
    boundary_filename = "";

    num_fp64_props = 0;
    num_int_props = 0;
    fp64_props = NULL;
    int_props = NULL;

    bc_init = NULL;
    bc_validate = NULL;
    bc_time_varying = NULL;

    generate_dirichlet_bcs = NULL;
    generate_node_number_override = NULL;
    bc_momentum = NULL;
    bc_force = NULL;
}

Boundary::Boundary(std::string filename, size_t nfp64, size_t nint, double *fp64props, int *intprops){
    std::string filepath = "src/boundaries/";
    filepath += filename;

    use_builtin = 0;
    handle = dlopen(filepath.c_str(), RTLD_LAZY);
    boundary_filename = filename;

    num_fp64_props = nfp64;
    num_int_props = nint;
    fp64_props = fp64props;
    int_props = intprops;

    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';
        bc_init = NULL;
        bc_validate = NULL;
        bc_time_varying = NULL;

        generate_dirichlet_bcs = NULL;
        generate_node_number_override = NULL;
        bc_momentum = NULL;
        bc_force = NULL;
    } else {
        dlerror();
        bc_init = reinterpret_cast<void (*)(job_t*)>(dlsym(handle,"bc_init"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'bc_init': " << dlsym_error <<
            '\n';
        }
        bc_validate = reinterpret_cast<void (*)(job_t*)>(dlsym(handle,"bc_validate"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'bc_validate': " << dlsym_error <<
            '\n';
        }
        bc_time_varying = reinterpret_cast<void (*)(job_t*)>(dlsym(handle,"bc_time_varying"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'bc_time_varying': " << dlsym_error <<
            '\n';
        }

        generate_dirichlet_bcs = reinterpret_cast<void (*)(job_t*)>(dlsym(handle,"generate_dirichlet_bcs"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'generate_dirichlet_bcs': " << dlsym_error <<
            '\n';
        }
        generate_node_number_override = reinterpret_cast<void (*)(job_t*)>(dlsym(handle,"generate_node_number_override"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'generate_node_number_override': " << dlsym_error <<
            '\n';
        }
        bc_momentum = reinterpret_cast<void (*)(job_t*)>(dlsym(handle,"bc_momentum"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'bc_momentum': " << dlsym_error <<
            '\n';
        }
        bc_force = reinterpret_cast<void (*)(job_t*)>(dlsym(handle,"bc_force"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'bc_force': " << dlsym_error <<
            '\n';
        }
    }
}

Boundary::~Boundary(){
    if(handle){
        dlclose(handle);
    }
}


void Boundary::setBoundary(std::string filename, size_t nfp64, size_t nint, double *fp64props, int *intprops){
    std::string filepath = "src/boundaries/";
    filepath += filename;

    this->use_builtin = 0;
    this->handle = dlopen(filepath.c_str(), RTLD_LAZY);
    this->boundary_filename = filename;

    this->num_fp64_props = nfp64;
    this->num_int_props = nint;
    this->fp64_props = fp64props;
    this->int_props = intprops;

    char* dlsym_error;
    if (!this->handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';
        this->bc_init = NULL;
        this->bc_validate = NULL;
        this->bc_time_varying = NULL;

        this->generate_dirichlet_bcs = NULL;
        this->generate_node_number_override = NULL;
        this->bc_momentum = NULL;
        this->bc_force = NULL;
    } else {
        dlerror();
        this->bc_init = reinterpret_cast<void (*)(job_t*)>(dlsym(this->handle,"bc_init"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'bc_init': " << dlsym_error <<
            '\n';
        }
        this->bc_validate = reinterpret_cast<void (*)(job_t*)>(dlsym(this->handle,"bc_validate"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'bc_validate': " << dlsym_error <<
            '\n';
        }
        this->bc_time_varying = reinterpret_cast<void (*)(job_t*)>(dlsym(this->handle,"bc_time_varying"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'bc_time_varying': " << dlsym_error <<
            '\n';
        }

        this->generate_dirichlet_bcs = reinterpret_cast<void (*)(job_t*)>(dlsym(this->handle,"generate_dirichlet_bcs"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'generate_dirichlet_bcs': " << dlsym_error <<
            '\n';
        }
        this->generate_node_number_override = reinterpret_cast<void (*)(job_t*)>(dlsym(this->handle,"generate_node_number_override"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'generate_node_number_override': " << dlsym_error <<
            '\n';
        }
        this->bc_momentum = reinterpret_cast<void (*)(job_t*)>(dlsym(this->handle,"bc_momentum"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'bc_momentum': " << dlsym_error <<
            '\n';
        }
        this->bc_force = reinterpret_cast<void (*)(job_t*)>(dlsym(this->handle,"bc_force"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'bc_force': " << dlsym_error <<
            '\n';
        }
    }
}
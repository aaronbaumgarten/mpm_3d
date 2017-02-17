//
// Created by aaron on 12/22/16.
// material.cpp
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include "material.hpp"
#include "process.hpp"
#include <dlfcn.h>

Material::Material(){
        num_fp64_props = 0;
        num_int_props = 0;
        fp64_props = std::vector<double>();
        int_props = std::vector<int>();
        material_filename = "";
        handle = NULL;

        material_init = NULL;
        calculate_stress = NULL;
        calculate_stress_threaded = NULL;
        calculate_stress_implicit = NULL;
        volumetric_smoothing = NULL;
        volumetric_smoothing_implicit = NULL;
}

Material::Material(std::string filename, std::string filepath, std::vector<double> fp64props, std::vector<int> intprops){
    //std::string filepath = "src/materials/";
    filepath += filename;

    use_builtin = 0;
    handle = dlopen(filepath.c_str(), RTLD_LAZY);
    material_filename = filename;

    num_fp64_props = fp64props.size();
    num_int_props = intprops.size();
    fp64_props = fp64props;
    int_props = intprops;

    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';
        material_init = NULL;
        calculate_stress = NULL;
        calculate_stress_threaded = NULL;
        calculate_stress_implicit = NULL;
        volumetric_smoothing = NULL;
        volumetric_smoothing_implicit = NULL;
    } else {
        dlerror();
        material_init = reinterpret_cast<void (*)(Body *)>(dlsym(handle, "material_init"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'material_init': " << dlsym_error <<
            '\n';
        }
        calculate_stress = reinterpret_cast<void (*)(Body *, double)>(dlsym(handle, "calculate_stress"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'calculate_stress': " << dlsym_error <<
            '\n';
        }
        calculate_stress_threaded = reinterpret_cast<void (*)(threadtask_t*, Body*, double)>(dlsym(handle, "calculate_stress_threaded"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'calculate_stress_threaded': " << dlsym_error <<
            '\n';
        }
        calculate_stress_implicit = reinterpret_cast<void (*)(Body *, double)>(dlsym(handle, "calculate_stress_implicit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'calculate_stress_implicit': " << dlsym_error <<
            '\n';
        }
        volumetric_smoothing = reinterpret_cast<void (*)(Body *, Eigen::VectorXd, Eigen::VectorXd)>(dlsym(handle, "volumetric_smoothing"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'volumetric_smoothing': " << dlsym_error <<
            '\n';
        }
        volumetric_smoothing_implicit = reinterpret_cast<void (*)(Body *, Eigen::VectorXd, Eigen::VectorXd)>(dlsym(handle, "volumetric_smoothing_implicit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'volumetric_smoothing_implicit': " << dlsym_error <<
            '\n';
        }
    }
}

Material::~Material() {
    //std::cout << "closing material" << std::endl;
    if(handle){
        dlclose(handle);
    }
}

void Material::setMaterial(std::string filename, std::string filepath, std::vector<double> fp64props, std::vector<int> intprops){
    //std::string filepath = "src/materials/";
    filepath += filename;

    this->use_builtin = 0;
    this->handle = dlopen(filepath.c_str(), RTLD_LAZY);
    this->material_filename = filename;

    this->num_fp64_props = fp64props.size();
    this->num_int_props = intprops.size();
    this->fp64_props = fp64props;
    this->int_props = intprops;

    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';
        this->material_init = NULL;
        this->calculate_stress = NULL;
        this->calculate_stress_threaded = NULL;
        this->calculate_stress_implicit = NULL;
        this->volumetric_smoothing = NULL;
        this->volumetric_smoothing_implicit = NULL;
    } else {
        dlerror();
        this->material_init = reinterpret_cast<void (*)(Body *)>(dlsym(this->handle, "material_init"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'material_init': " << dlsym_error <<
            '\n';
        }
        this->calculate_stress = reinterpret_cast<void (*)(Body *, double)>(dlsym(this->handle, "calculate_stress"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'calculate_stress': " << dlsym_error <<
            '\n';
        }
        this->calculate_stress_threaded = reinterpret_cast<void (*)(threadtask_t *, Body *, double)>(dlsym(this->handle, "calculate_stress_threaded"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'calculate_stress_threaded': " << dlsym_error <<
            '\n';
        }
        this->calculate_stress_implicit = reinterpret_cast<void (*)(Body *, double)>(dlsym(this->handle, "calculate_stress_implicit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'calculate_stress_implicit': " << dlsym_error <<
            '\n';
        }
        this->volumetric_smoothing = reinterpret_cast<void (*)(Body *,Eigen::VectorXd, Eigen::VectorXd)>(dlsym(handle, "volumetric_smoothing"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'volumetric_smoothing': " << dlsym_error <<
            '\n';
        }
        this->volumetric_smoothing_implicit = reinterpret_cast<void (*)(Body *, Eigen::VectorXd, Eigen::VectorXd)>(dlsym(handle, "volumetric_smoothing_implicit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'volumetric_smoothing_implicit': " << dlsym_error <<
            '\n';
        }
    }
}
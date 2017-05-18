//
// Created by aaron on 5/14/17.
// grid.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <dlfcn.h>

#include "job.hpp"

#include "serializer.hpp"
#include "grid.hpp"

Grid::Grid() {
    filename = "";
    filepath = "";
    fp64_props = std::vector<double>();
    int_props = std::vector<int>();

    handle = NULL;

    gridInit = NULL;
    gridWhichElement = NULL;
    gridNodeIDToPosition = NULL;
    gridEvaluateShapeFnValue = NULL;
    gridEvaluateShapeFnGradient = NULL;

    gridWriteFrame = NULL;
    gridSaveState = NULL;
    gridLoadState = NULL;
}

Grid::~Grid() {
    if (handle){
        dlclose(handle);
    }
}

void Grid::gridSetPlugin(Job* job, std::string nameIN, std::string pathIN, std::vector<double> fp64IN, std::vector<int> intIN){
    filename = nameIN;
    filepath = pathIN;
    fp64_props = fp64IN;
    int_props = intIN;

    handle = dlopen((filepath+filename).c_str(), RTLD_LAZY);

    gridSetFnPointers(handle);

    return;
}

void Grid::gridSetFnPointers(void* handle){
    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';

        gridInit = NULL;
        gridWhichElement = NULL;
        gridNodeIDToPosition = NULL;
        gridEvaluateShapeFnValue = NULL;
        gridEvaluateShapeFnGradient = NULL;

        gridWriteFrame = NULL;
        gridSaveState = NULL;
        gridLoadState = NULL;

    } else {
        dlerror();
        gridInit = reinterpret_cast<void (*)(Job *)>(dlsym(handle, "gridInit"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridInit': " << dlsym_error <<
                      '\n';
        }
        gridWhichElement = reinterpret_cast<int (*)(Job *, Eigen::VectorXd)>(dlsym(handle, "gridWhichElement"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridWhichElement': " << dlsym_error <<
                      '\n';
        }
        gridNodeIDToPosition = reinterpret_cast<Eigen::VectorXd (*)(Job *, int)>(dlsym(handle, "gridNodeIDToPosition"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridNodeIDToPosition': " << dlsym_error <<
                      '\n';
        }
        gridEvaluateShapeFnValue = reinterpret_cast<void (*)(Job *, Eigen::VectorXd, Eigen::VectorXi*, Eigen::VectorXd*)>(dlsym(handle, "gridEvaluateShapeFnValue"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridEvaluateShapeFnValue': " << dlsym_error <<
                      '\n';
        }
        gridEvaluateShapeFnGradient = reinterpret_cast<void (*)(Job *, Eigen::VectorXd, Eigen::VectorXi*, Eigen::MatrixXd*)>(dlsym(handle, "gridEvaluateShapeFnGradient"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridEvaluateShapeFnGradient': " << dlsym_error <<
                      '\n';
        }

        gridWriteFrame = reinterpret_cast<void (*)(Job*,Serializer*)>(dlsym(handle, "contactWriteFrame"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contactWriteFrame': " << dlsym_error <<
                      '\n';
        }
        gridSaveState = reinterpret_cast<std::string (*)(Job*,Serializer*,std::string)>(dlsym(handle, "contactSaveState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contactSaveState': " << dlsym_error <<
                      '\n';
        }
        gridLoadState = reinterpret_cast<int (*)(Job*,Serializer*,std::string)>(dlsym(handle, "contactLoadState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'contactLoadState': " << dlsym_error <<
                      '\n';
        }
    }
    return;
}
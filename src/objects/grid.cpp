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

#include "stringparser.hpp"
#include "job.hpp"

#include "serializer.hpp"
#include "grid.hpp"

Grid::Grid() {
    fullpath = "";
    filename = "";
    filepath = "";
    fp64_props = std::vector<double>();
    int_props = std::vector<int>();
    str_props = std::vector<std::string>();

    node_count = 0;
    element_count = 0;

    handle = NULL;

    gridInit = NULL;
    gridWhichElement = NULL;
    gridInDomain = NULL;
    gridNodeIDToPosition = NULL;
    gridEvaluateShapeFnValue = NULL;
    gridEvaluateShapeFnGradient = NULL;
    gridNodalVolume = NULL;
    gridElementVolume = NULL;

    gridWriteFrame = NULL;
    gridSaveState = NULL;
    gridLoadState = NULL;
}

Grid::Grid(const Grid& obj){
    fullpath = obj.fullpath;
    filename = obj.filename;
    filepath = obj.filepath;
    fp64_props = obj.fp64_props;
    int_props = obj.int_props;
    str_props = obj.str_props;

    node_count = obj.node_count;
    element_count = obj.element_count;


    handle = NULL;

    gridInit = NULL;
    gridWhichElement = NULL;
    gridInDomain = NULL;
    gridNodeIDToPosition = NULL;
    gridEvaluateShapeFnValue = NULL;
    gridEvaluateShapeFnGradient = NULL;
    gridNodalVolume = NULL;
    gridElementVolume = NULL;

    gridWriteFrame = NULL;
    gridSaveState = NULL;
    gridLoadState = NULL;
}

Grid::~Grid() {
    if (handle){
        dlclose(handle);
    }
}

void Grid::gridSetPlugin(Job* job, std::string pathIN, std::string nameIN, std::vector<double> fp64IN, std::vector<int> intIN, std::vector<std::string> strIN){
    filename = nameIN;
    fullpath = StringParser::stringMakeDirectory(pathIN);
    fp64_props = fp64IN;
    int_props = intIN;
    str_props = strIN;

    gridSetFnPointers();

    return;
}

void Grid::gridSetFnPointers(){
    if (!handle) {
        handle = dlopen((fullpath + filename).c_str(), RTLD_LAZY);
    }

    char* dlsym_error;
    if (!handle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';

        gridInit = NULL;
        gridWhichElement = NULL;
        gridInDomain = NULL;
        gridNodeIDToPosition = NULL;
        gridEvaluateShapeFnValue = NULL;
        gridEvaluateShapeFnGradient = NULL;
        gridNodalVolume = NULL;
        gridElementVolume = NULL;

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
        gridInDomain = reinterpret_cast<bool (*)(Job *, Eigen::VectorXd)>(dlsym(handle, "gridInDomain"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridInDomain': " << dlsym_error <<
                      '\n';
        }
        gridNodeIDToPosition = reinterpret_cast<Eigen::VectorXd (*)(Job *, int)>(dlsym(handle, "gridNodeIDToPosition"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridNodeIDToPosition': " << dlsym_error <<
                      '\n';
        }
        gridEvaluateShapeFnValue = reinterpret_cast<void (*)(Job *, Eigen::VectorXd, std::vector<int>&, std::vector<double>&)>(dlsym(handle, "gridEvaluateShapeFnValue"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridEvaluateShapeFnValue': " << dlsym_error <<
                      '\n';
        }
        gridEvaluateShapeFnGradient = reinterpret_cast<void (*)(Job *, Eigen::VectorXd, std::vector<int>&, std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>&)>(dlsym(handle, "gridEvaluateShapeFnGradient"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridEvaluateShapeFnGradient': " << dlsym_error <<
                      '\n';
        }
        gridNodalVolume = reinterpret_cast<double (*)(Job*,int)>(dlsym(handle, "gridNodalVolume"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridNodalVolume': " << dlsym_error <<
                      '\n';
        }
        gridElementVolume = reinterpret_cast<double (*)(Job*,int)>(dlsym(handle, "gridElementVolume"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridElementVolume': " << dlsym_error <<
                      '\n';
        }

        gridWriteFrame = reinterpret_cast<void (*)(Job*,Serializer*)>(dlsym(handle, "gridWriteFrame"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridWriteFrame': " << dlsym_error <<
                      '\n';
        }
        gridSaveState = reinterpret_cast<std::string (*)(Job*,Serializer*,std::string)>(dlsym(handle, "gridSaveState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridSaveState': " << dlsym_error <<
                      '\n';
        }
        gridLoadState = reinterpret_cast<int (*)(Job*,Serializer*,std::string)>(dlsym(handle, "gridLoadState"));
        dlsym_error = dlerror();
        if (dlsym_error) {
            std::cerr << "Cannot load symbol 'gridLoadState': " << dlsym_error <<
                      '\n';
        }
    }
    return;
}
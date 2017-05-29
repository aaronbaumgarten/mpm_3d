//
// Created by aaron on 5/6/17.
// boundary.hpp
//

#ifndef MPM_3D_BOUNDARY_HPP
#define MPM_3D_BOUNDARY_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

#include "runtimedef.hpp"

class Job;
class Serializer;
class Body;

class Boundary: public RunTimeDef{
public:
    //boundary properties here
    std::string filename; //name of file
    std::string filepath; //directory of file for access
    std::vector<double> fp64_props; //double properties
    std::vector<int> int_props; //integer properties
    std::vector<std::string> str_props; //string properties
    void *handle; //.so file handle

    //boundary specific functions
    Boundary();
    ~Boundary();
    void boundarySetPlugin(Job*, Body*, std::string, std::string, std::vector<double>, std::vector<int>, std::vector<std::string>); // set .so file
    void boundarySetFnPointers(void*); //set function pointers to .so file handle

    void (*boundaryWriteFrame)(Job*, Body*, Serializer*); //write frame to serializer
    std::string (*boundarySaveState)(Job*, Body*, Serializer*, std::string); //save state to serializer folder with returned filename
    int (*boundaryLoadState)(Job*, Body*, Serializer*, std::string); //read state from given full path

    void (*boundaryInit)(Job*,Body*); //initialize boundary object
    void (*boundaryGenerateRules)(Job*,Body*); //generate the rules given job and body state
    void (*boundaryApplyRules)(Job*,Body*); //apply the rules given job and body state

    //other functions
};

#endif //MPM_3D_BOUNDARY_HPP

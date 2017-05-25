//
// Created by aaron on 5/6/17.
// material.hpp
//

#ifndef MPM_3D_MATERIAL_HPP
#define MPM_3D_MATERIAL_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

class Job;
class Serializer;
class Body;

class Material{
public:
    //material properties here
    std::string filename; //name of file
    std::string filepath; //directory of file for access
    std::vector<double> fp64_props; //double properties
    std::vector<int> int_props; //integer properties
    std::vector<std::string> str_props; //string properties
    void *handle; //.so file handle

    //material specific functions
    Material();
    ~Material();
    void materialSetPlugin(Job*, Body*, std::string, std::string, std::vector<double>, std::vector<int>, std::vector<std::string>); //set .so material file
    void materialSetFnPointers(void*); //set function pointers to .so file handle

    void (*materialWriteFrame)(Job*, Body*, Serializer*); //write frame to serializer
    std::string (*materialSaveState)(Job*, Body*, Serializer*, std::string); //save state to serializer folder with returned filename
    int (*materialLoadState)(Job*, Body*, Serializer*, std::string); //read state from given full path

    void (*materialInit)(Job*,Body*); // initialize material
    void (*materialCalculateStress)(Job*,Body*,int=1); // calculate stress given body and job state, int=1 updates internal variables
    void (*materialAssignStress)(Job*,Body*,Eigen::MatrixXd,int); // assign stress state to specific particle id
    void (*materialAssignPressure)(Job*,Body*,double,int); // assign pressure state to specific particle

    //other functions
};


#endif //MPM_3D_MATERIAL_HPP

//
// Created by aaron on 5/11/17.
// contact.hpp
//

#ifndef MPM_3D_CONTACT_HPP
#define MPM_3D_CONTACT_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

#include "runtimedef.hpp"

class Job;
class Serializer;

class Contact: public RunTimeDef{
public:
    //boundary properties here
    int id;
    std::string name;
    /*std::string fullpath; //long path to file
    std::string filename; //name of file
    std::string filepath; //directory of file for access
    std::vector<double> fp64_props; //double properties
    std::vector<int> int_props; //integer properties
    std::vector<std::string> str_props; //string properties
    void *handle; //.so file handle*/

    //boundary specific functions
    Contact();
    ~Contact();
    void contactSetPlugin(Job*, std::string, std::string, std::vector<double>, std::vector<int>, std::vector<std::string>); // set .so file
    void contactSetFnPointers(void*); //set function pointers to .so file handle

    void (*contactWriteFrame)(Job*, Serializer*); //write frame to serializer
    std::string (*contactSaveState)(Job*, Serializer*, std::string); //save state to serializer folder with returned filename
    int (*contactLoadState)(Job*, Serializer*, std::string); //read state from given full path

    void (*contactInit)(Job*, Contact*); //initialize contact (need access to it initially)
    void (*contactGenerateRules)(Job*); //generate contact rules
    void (*contactApplyRules)(Job*); //apply rules

    //other functions
};

#endif //MPM_3D_CONTACT_HPP

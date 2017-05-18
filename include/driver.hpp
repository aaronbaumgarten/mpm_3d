//
// Created by aaron on 5/11/17.
// driver.hpp
//

#ifndef MPM_3D_DRIVER_HPP
#define MPM_3D_DRIVER_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

class Job;
class Serializer;

class Driver{
public:
    //driver properties here
    std::string filename; //name of file
    std::string filepath; //directory of file for access
    std::vector<double> fp64_props; //double properties
    std::vector<int> int_props; //integer properties
    void *handle; //.so file handle

    //boundary specific functions
    Driver();
    ~Driver();
    void driverSetPlugin(Job*, std::string, std::string, std::vector<double>, std::vector<int>); //assign .so plugin file
    void driverSetFnPointers(void*); //set function pointers to .so file handle

    void (*driverInit)(Job*); //initialize driver
    void (*driverRun)(Job*); //run simulation

    std::string (*driverSaveState)(Job*,Serializer*, std::string); //save driver state to returned filename in serializer folder
    int (*driverLoadState)(Job*,Serializer*,std::string); //load state from given full path

    //other functions
};

#endif //MPM_3D_DRIVER_HPP

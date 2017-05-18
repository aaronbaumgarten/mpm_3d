//
// Created by aaron on 5/5/17.
// serializer.hpp
//

#ifndef MPM_3D_SERIALIZER_HPP
#define MPM_3D_SERIALIZER_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

class Job;

class Serializer{
public:
    //io properties here
    std::string filename; //name of file
    std::string filepath; //directory of file for access
    std::vector<double> fp64_props; //double properties
    std::vector<int> int_props; //integer properties
    void *handle; //.so file handle

    //objects here

    //io object specific functions
    Serializer();
    ~Serializer();
    void serializerSetPlugin(Job*, std::string, std::string, std::vector<double>, std::vector<int>); //assign .so plugin for functions
    void serializerSetFnPointers(void*); //set function pointers to .so file handle

    //input output (not configuration)
    void (*serializerWriteFrame)(Job*); //initialize writing of frame from job
    void (*serializerWriteScalar)(Eigen::Matrix*, std::string); //call to functions for dumping state into frame file
    void (*serializerWriteVector)(Eigen::Matrix*, std::string); //pass name of vector
    void (*serializerWriteTensor)(Eigen::Matrix*, std::string);

    void (*serializerInit)(Job*); //initialize serializer
    std::string (*serializerSaveState)(Job*); //save job state (output string for output directory)
    int (*serializerLoadState)(Job*,std::string); //load state from given full path

};

#endif //MPM_3D_SERIALIZER_HPP

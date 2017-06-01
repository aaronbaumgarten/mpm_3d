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

#include "runtimedef.hpp"

class Job;

class Serializer: public RunTimeDef{
public:
    //io properties here
    /*std::string fullpath; //long path to file
    std::string filename; //name of file
    std::string filepath; //directory of file for access
    std::vector<double> fp64_props; //double properties
    std::vector<int> int_props; //integer properties
    std::vector<std::string> str_props; //string properties
    void *handle; //.so file handle*/

    //other info
    std::string mainpath;

    //io object specific functions
    Serializer();
    ~Serializer();
    void serializerSetPlugin(Job*, std::string, std::string, std::vector<double>, std::vector<int>, std::vector<std::string>); //assign .so plugin for functions
    void serializerSetFnPointers(void*); //set function pointers to .so file handle
    void serializerSetMainPath(Job*,std::string); //set path to main program

    //input output (not configuration)
    int (*serializerWriteFrame)(Job*); //initialize writing of frame from job (return 1 if frame written)
    void (*serializerWriteScalarArray)(Eigen::VectorXd&, std::string); //call to functions for dumping state into frame file
    void (*serializerWriteVectorArray)(Eigen::MatrixXd&, std::string); //pass name of vector
    void (*serializerWriteTensorArray)(Eigen::MatrixXd&, std::string);

    void (*serializerInit)(Job*); //initialize serializer
    std::string (*serializerSaveState)(Job*); //save job state (output string for output directory)
    int (*serializerLoadState)(Job*,std::string); //load state from given full path

};

#endif //MPM_3D_SERIALIZER_HPP

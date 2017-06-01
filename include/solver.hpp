//
// Created by aaron on 5/11/17.
// solver.hpp
//

#ifndef MPM_3D_SOLVER_HPP
#define MPM_3D_SOLVER_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

#include "runtimedef.hpp"

class Job;
class Serializer;

class Solver: public RunTimeDef{
public:
    //solver properties here
    /*std::string fullpath; //long path to file
    std::string filename; //name of file
    std::string filepath; //directory of file for access
    std::vector<double> fp64_props; //double properties
    std::vector<int> int_props; //integer properties
    std::vector<std::string> str_props; //string properties
    void *handle; //.so file handle*/

    //solver specific functions
    Solver();
    ~Solver();
    void solverSetPlugin(Job*, std::string, std::string, std::vector<double>, std::vector<int>, std::vector<std::string>); //assign .so plugin file
    void solverSetFnPointers(void*); //set function pointers to .so file handle

    void (*solverInit)(Job*); //initialize solver
    void (*solverStep)(Job*); //one forward mpm step

    std::string (*solverSaveState)(Job*,Serializer*,std::string); //save solver state to returned filename in serializer folder
    int (*solverLoadState)(Job*,Serializer*,std::string); //load state from given full path

    //other functions
};

#endif //MPM_3D_SOLVER_HPP

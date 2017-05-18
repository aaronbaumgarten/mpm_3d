//
// Created by aaron on 5/5/17.
// job.hpp
//

#ifndef MPM_3D_JOB_HPP
#define MPM_3D_JOB_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

class Serializer;
class Driver;
class Solver;
class Grid;
class Body;
class Contact;

class Job{
public:
    //job properties here
    int DIM;
    int XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ;

    std::vector<int> activeBodies; //reference map for active bodies (1 = active, 0 = inactive)
    std::vector<int> activeContacts; //reference map for active contacts

    //objects here
    Serializer serializer; // input output
    Driver driver; // simulation driver
    Solver solver; // solver type (explicit, implicit newton-cg, etc.)
    Grid grid; // element types and mappings
    std::vector<Body> bodies; // vector of simulated bodies
    std::vector<Contact> contacts; // vector of contacts to monitor

    //job object specific functions
    Job();
    int jobInit();

    //templates of vectors and tensors for consistency
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> jobVector();
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> jobVectorArray(int);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> jobTensor();
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> jobTensorArray(int);
};

#endif //MPM_3D_JOB_HPP

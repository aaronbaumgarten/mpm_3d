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
#include <iostream>
#include <fstream>

class Serializer;
class Driver;
class Solver;
class Grid;
class Body;
class Contact;

class Job{
public:
    //static integers
    static int ZERO = 0;
    static int ONES = 1;
    static int IDENTITY = 2;
    //job properties here
    int DIM;
    int XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ;
    int X, Y, Z;
    double t, dt;

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
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::RowMajor> jobVector();
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::RowMajor> jobVector(T*);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::RowMajor> jobVector(int);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jobVectorArray(int);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jobTensor();
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jobTensor(T*);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jobTensor(int);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jobTensorArray(int);

    void jobScalarArrayToFile(Eigen::Matrix &, std::ostream &);//write matrix ref to file ref
    void jobVectorArrayToFile(Eigen::Matrix &, std::ostream &);
    void jobTensorArrayToFile(Eigen::Matrix &, std::ostream &);

    void jobScalarArrayFromFile(Eigen::Matrix &, std::istream &);//read matrix from file ref
    void jobVectorArrayFromFile(Eigen::Matrix &, std::istream &);
    void jobTensorArrayFromFile(Eigen::Matrix &, std::istream &);
};

#endif //MPM_3D_JOB_HPP

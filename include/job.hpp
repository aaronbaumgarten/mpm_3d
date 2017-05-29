//
// Created by aaron on 5/5/17.
// job.hpp
//

#ifndef MPM_3D_JOB_HPP
#define MPM_3D_JOB_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <iostream>
#include <fstream>

#include "serializer.hpp"
#include "driver.hpp"
#include "solver.hpp"
#include "grid.hpp"

class Body;
class Contact;

class Job{
public:
    //static integers
    static const int ZERO = 0;
    static const int ONES = 1;
    static const int IDENTITY = 2;
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
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> jobVector();
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> jobVector(T*);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> jobVector(int);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jobVectorArray(int);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jobTensor();
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jobTensor(T*);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jobTensor(int);
    template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> jobTensorArray(int);

    template <typename Derived> void jobScalarArrayToFile(Eigen::MatrixBase<Derived>& , std::ostream &);//write matrix ref to file ref
    template <typename Derived> void jobVectorArrayToFile(Eigen::MatrixBase<Derived>& , std::ostream &);
    template <typename Derived> void jobTensorArrayToFile(Eigen::MatrixBase<Derived>& , std::ostream &);

    template <typename Derived> void jobScalarArrayFromFile(Eigen::MatrixBase<Derived>& , std::istream &);//read matrix from file ref
    template <typename Derived> void jobVectorArrayFromFile(Eigen::MatrixBase<Derived>& , std::istream &);
    template <typename Derived> void jobTensorArrayFromFile(Eigen::MatrixBase<Derived>& , std::istream &);
};


template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> Job::jobVector(){
    return Eigen::Matrix<T, Eigen::Dynamic, 1>(DIM);
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> Job::jobVector(T* data){
    Eigen::Matrix<T, Eigen::Dynamic, 1> tmpVec(DIM);
    for (size_t i=0;i<tmpVec.size();i++){
        tmpVec(i) = data[i];
    }
    return tmpVec;
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> Job::jobVector(int SPEC){
    if (SPEC == ZERO) {
        return Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(DIM);
    } else if (SPEC == ONES) {
        return Eigen::Matrix<T, Eigen::Dynamic, 1>::Ones(DIM);
    } else {
        std::cerr << "jobVector(SPEC): Unknown SPEC [" << SPEC << "]" << std::endl;
        return Eigen::Matrix<T, Eigen::Dynamic, 1>(DIM);
    }
};


template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Job::jobVectorArray(int len){
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(DIM,len);
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Job::jobTensor(){
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(DIM,DIM);
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Job::jobTensor(T* data){
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmpMat(DIM,DIM);
    for (size_t i=0;i<tmpMat.size();i++){
        tmpMat(i) = data[i];
    }
    return tmpMat;
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Job::jobTensor(int SPEC){
    if (SPEC == ZERO) {
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(DIM, DIM);
    } else if (SPEC == ONES) {
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Ones(DIM, DIM);
    } else if (SPEC == IDENTITY){
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Identity(DIM, DIM);
    } else {
        std::cerr << "jobTensor(SPEC): Unknown SPEC [" << SPEC << "]" << std::endl;
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(DIM, DIM);
    }
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Job::jobTensorArray(int len){
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>(DIM*DIM,len);
};

template <typename Derived> void Job::jobScalarArrayToFile(Eigen::MatrixBase<Derived>& x, std::ostream &ffile){
    assert(x.cols() == 1);

    for (size_t i = 0; i < x.rows(); i++) {
        ffile << x(i) << "\n";
    }
    return;
}

template <typename Derived> void Job::jobVectorArrayToFile(Eigen::MatrixBase<Derived>& x, std::ostream &ffile){
    assert(x.cols() >= DIM);

    for (size_t i = 0; i < x.rows(); i++) {
        for (size_t pos = 0; pos < DIM; pos++) {
            ffile << x(i, pos) << " ";
        }
        ffile << "\n";
    }
    return;
}

template <typename Derived> void Job::jobTensorArrayToFile(Eigen::MatrixBase<Derived>& x, std::ostream &ffile){
    assert(x.cols() >= DIM*DIM);

    for (size_t i = 0; i < x.rows(); i++) {
        for (size_t pos = 0; pos < DIM*DIM; pos++) {
            ffile << x(i, pos) << " ";
        }
        ffile << "\n";
    }
    return;
}

template <typename Derived> void Job::jobScalarArrayFromFile(Eigen::MatrixBase<Derived>& x, std::istream &ffile){
    //x must be preformed to correct dimensions
    assert(x.cols()==1);

    std::string line;
    std::stringstream ss;

    for (size_t i=0;i<x.rows();i++){
        if (std::getline(ffile, line)){
            ss = std::stringstream(line);
            ss >> x(i);
        } else {
            std::cerr << "\nJob::job<field>FromFile Error:\nFile ended before array was filled!\n" << std::endl;
            return;
        }
    }
    return;
}

template <typename Derived> void Job::jobVectorArrayFromFile(Eigen::MatrixBase<Derived>& x, std::istream &ffile){
    //x must be preformed to correct dimensions
    assert(x.cols()>=DIM);

    std::string line;
    std::stringstream ss;

    for (size_t i=0;i<x.rows();i++){
        if (std::getline(ffile, line)){
            ss = std::stringstream(line);
            for (size_t pos=0;pos<DIM;pos++) {
                ss >> x(i,pos);
            }
        } else {
            std::cerr << "\nJob::job<field>FromFile Error:\nFile ended before array was filled!\n" << std::endl;
            return;
        }
    }
    return;
}

template <typename Derived> void Job::jobTensorArrayFromFile(Eigen::MatrixBase<Derived>& x, std::istream &ffile){
    //x must be preformed to correct dimensions
    assert(x.cols()>=DIM*DIM);

    std::string line;
    std::stringstream ss;

    for (size_t i=0;i<x.rows();i++){
        if (std::getline(ffile, line)){
            ss = std::stringstream(line);
            for (size_t pos=0;pos<DIM*DIM;pos++) {
                ss >> x(i,pos);
            }
        } else {
            std::cerr << "\nJob::job<field>FromFile Error:\nFile ended before array was filled!\n" << std::endl;
            return;
        }
    }
    return;
}

#endif //MPM_3D_JOB_HPP

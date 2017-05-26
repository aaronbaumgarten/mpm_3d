//
// Created by aaron on 5/11/17.
// job.cpp
//

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "stringparser.hpp"

#include "job.hpp"

#include "serializer.hpp"
#include "driver.hpp"
#include "solver.hpp"

#include "grid.hpp"
#include "contact.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "points.hpp"

#include "material.hpp"
#include "boundary.hpp"

Job::Job():
        activeBodies(0),
        activeContacts(0),
        bodies(0),
        contacts(0)
{
    DIM = 3;
    XX = 0; XY = 1; XZ = 2;
    YX = 3; YY = 4; YZ = 5;
    ZX = 6; ZY = 7; ZZ = 8;

    X = 0; Y = 1; Z = 2;

    t = 0; dt = 1e-3;

    serializer = Serializer();
    driver = Driver();
    solver = Solver();

    grid = Grid();
}

int Job::jobInit(){

    serializer.serializerInit(this);
    driver.driverInit(this);
    solver.solverInit(this);

    grid.gridInit(this);

    for (size_t i=0; i<contacts.size(); i++){
        contacts[i].contactInit(this);
    }
    for (size_t i=0; i<bodies.size(); i++){
        bodies[i].bodyInit(this);
    }

    return 1;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::RowMajor> Job::jobVector(){
    return Eigen::Matrix<T, Eigen::Dynamic, 1>(DIM);
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::RowMajor> Job::jobVector(T* data){
    return Eigen::Matrix<T, Eigen::Dynamic, 1>(data,DIM);
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::RowMajor> Job::jobVector(int SPEC){
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
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(DIM,len);
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Job::jobTensor(){
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(DIM,DIM);
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Job::jobTensor(T* data){
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(data,DIM,DIM);
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Job::jobTensor(int SPEC){
    if (SPEC == ZERO) {
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(DIM, DIM);
    } else if (SPEC == ONES) {
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Ones(DIM, DIM);
    } else if (SPEC == IDENTITY){
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(DIM, DIM);
    } else {
        std::cerr << "jobTensor(SPEC): Unknown SPEC [" << SPEC << "]" << std::endl;
        return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(DIM, DIM);
    }
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Job::jobTensorArray(int len){
    return Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(DIM*DIM,len);
};

void Job::jobScalarToFile(Eigen::Matrix& x,std::ostream& ffile){
    assert(x.cols() == 1);

    for (size_t i = 0; i < x.rows(); i++) {
        ffile << x(i) << "\n";
    }
    return;
}

void Job::jobVectorToFile(Eigen::Matrix& x,std::ostream& ffile){
    assert(x.cols() >= DIM);

    for (size_t i = 0; i < x.rows(); i++) {
        for (size_t pos = 0; pos < DIM; pos++) {
            ffile << x(i, pos) << " ";
        }
        ffile << "\n";
    }
    return;
}

void Job::jobTensorToFile(Eigen::Matrix& x,std::ostream& ffile){
    assert(x.cols() >= DIM*DIM);

    for (size_t i = 0; i < x.rows(); i++) {
        for (size_t pos = 0; pos < DIM*DIM; pos++) {
            ffile << x(i, pos) << " ";
        }
        ffile << "\n";
    }
    return;
}

void Job::jobScalarFromFile(Eigen::Matrix& x,std::istream& ffile){
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

void Job::jobVectorFromFile(Eigen::Matrix& x,std::istream& ffile){
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

void Job::jobTensorFromFile(Eigen::Matrix& x,std::istream& ffile){
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
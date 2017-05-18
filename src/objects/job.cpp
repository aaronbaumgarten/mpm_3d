//
// Created by aaron on 5/11/17.
// job.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

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

    serializer = Serializer();
    driver = Driver();
    solver = Solver();

    grid = Grid();
}

int Job::jobInit(){
    //job initialization stuff
    //kind of useless right now
    //use for stuff that can't be in constructor
    return 1;
}

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> Job::jobVector(){
    Eigen::Matrix<T, Eigen::Dynamic, 1> jvec(DIM);
    jvec.setZero();
    return jvec;
};


template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Job::jobVectorArray(int len){
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> jvec(DIM,len);
    jvec.setZero();
    return jvec;
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Job::jobTensor(){
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> jmat(DIM,DIM);
    jmat.setZero();
    return jmat;
};

template <typename T> Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Job::jobTensorArray(int len){
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> jmat(DIM*DIM,len);
    jmat.setZero();
    return jmat;
};
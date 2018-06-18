//
// Created by aaron on 5/15/18.
// default_nodes.cpp
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>
#include <sstream>
#include <Eigen/Core>
#include <ctime>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"
#include "nodes.hpp"

/*----------------------------------------------------------------------------*/
//initialize node object assuming that other higher objects are initialized first
//job, serializer, driver, solver, body, grid all need to be initialized
void DefaultNodes::init(Job* job, Body* body) {
    //assume that grid object has been created first
    int len = job->grid->node_count;

    //initialize vectors to length of grid object
    if (x.size() != len) {
        x = KinematicVectorArray(len, job->JOB_TYPE);
        u = KinematicVectorArray(len, job->JOB_TYPE);
        x_t = KinematicVectorArray(len, job->JOB_TYPE);
        diff_x_t = KinematicVectorArray(len, job->JOB_TYPE);
        m.resize(len);
        mx_t = KinematicVectorArray(len, job->JOB_TYPE);
        f = KinematicVectorArray(len, job->JOB_TYPE);
        active.resize(len);
    }

    x.setZero();
    u.setZero();
    x_t.setZero();
    diff_x_t.setZero();

    m.setZero();
    mx_t.setZero();
    f.setZero();

    for (int i = 0; i < len; i++) {
        x(i) = job->grid->nodeIDToPosition(job, i);
    }

    active.setOnes(); //all nodes active

    std::cout << "Nodes Initialized: [" << body->name << "]." << std::endl;

    return;
}


/*----------------------------------------------------------------------------*/
//write frame to serializer
void DefaultNodes::writeFrame(Job* job, Body* body, Serializer* serializer) {
    //serializer will use x-position to create format for file
    //serializer.serializerWriteVectorArray(&x, "position")
    serializer->writeVectorArray(u, "displacement");
    serializer->writeVectorArray(x_t, "velocity");
    serializer->writeVectorArray(diff_x_t, "diff_velocity");
    serializer->writeScalarArray(m, "mass");
    serializer->writeVectorArray(mx_t, "momentum");
    serializer->writeVectorArray(f, "force");
    //need to make double
    Eigen::VectorXd tmpVec = active.cast<double>();
    serializer->writeScalarArray(tmpVec, "active");

    return;
}


/*----------------------------------------------------------------------------*/
//not implemented yet
std::string DefaultNodes::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int DefaultNodes::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//
void DefaultNodes::generateLoads(Job* job, Body* body){
    //do nothing
    return;
}

void DefaultNodes::applyLoads(Job* job, Body* body){
    //do nothing
    return;
}
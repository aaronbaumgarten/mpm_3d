//
// Created by aaron on 5/16/18.
// default_body.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "bodies.hpp"

/*----------------------------------------------------------------------------*/
//
void DefaultBody::init(Job* job){
    points->init(job,this);
    nodes->init(job,this);

    if (activeMaterial != 0) {
        material->init(job, this);
    }
    if (activeBoundary != 0) {
        boundary->init(job, this);
    }

    //assign vector type
    S = MPMScalarSparseMatrix(nodes->x.size(), points->x.size());
    gradS = KinematicVectorSparseMatrix(nodes->x.size(), points->x.size(), job->JOB_TYPE);
    std::cout << "Body Initialized: [" << name << "]." << std::endl;
    return;
}


/*----------------------------------------------------------------------------*/
//
std::string DefaultBody::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}


/*----------------------------------------------------------------------------*/
//
int DefaultBody::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//
void DefaultBody::generateMap(Job *job, int SPEC) {
    S.clear();
    gradS.clear();
    points->generateMap(job, this, SPEC);
}


/*----------------------------------------------------------------------------*/
//
void DefaultBody::generateLoads(Job* job){
    points->generateLoads(job, this);
    nodes->generateLoads(job, this);
    return;
}

void DefaultBody::applyLoads(Job* job){
    points->applyLoads(job, this);
    nodes->applyLoads(job, this);
    return;
}
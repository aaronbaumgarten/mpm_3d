//
// Created by aaron on 11/29/18.
// assigned_velocity.cpp
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
#include "boundaries.hpp"

/*----------------------------------------------------------------------------*/
//initialize boundary assuming other objects have been constructed correctly
void AssignedVelocity::init(Job* job, Body* body){
    if (fp64_props.size() < job->grid->GRID_DIM) {
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 1 properties defined (<velocity>).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //set velocity
        v_set = KinematicVector(job->JOB_TYPE);
        for (int i=0; i<job->grid->GRID_DIM; i++){
            v_set[i] = fp64_props[i];
        }

        //print grid properties
        std::cout << "Boundary properties (velocity = " << EIGEN_MAP_OF_KINEMATIC_VECTOR(v_set).transpose() << ")." << std::endl;
    }

    //set initial velocity of points
    for (int i=0; i<body->points->x.size(); i++){
        body->points->x_t(i) = v_set;
        body->points->mx_t(i) = body->points->m(i) * v_set;
    }

    std::cout << "Boundary Initialized: [" << body->name << "]." << std::endl;

    return;
}


/*----------------------------------------------------------------------------*/
//generate boundary conditions
void AssignedVelocity::generateRules(Job* job, Body* body){
    //nothing to do here
    return;
}


/*----------------------------------------------------------------------------*/
//apply boundary rules
void AssignedVelocity::applyRules(Job* job, Body* body){
    for (int i=0;i<body->nodes->x_t.size();i++){
        if (body->nodes->m(i) > 0) {
            //set velocity
            body->nodes->x_t(i) = v_set;
            body->nodes->mx_t(i) = body->nodes->m(i) * v_set;
            body->nodes->f(i).setZero();
        }
    }
    return;
}


/*----------------------------------------------------------------------------*/
//write frame and save file
void AssignedVelocity::writeFrame(Job* job, Body* body, Serializer* serializer){
    //nothing to report
    return;
}

std::string AssignedVelocity::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int AssignedVelocity::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}

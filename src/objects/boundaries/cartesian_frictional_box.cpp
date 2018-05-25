//
// Created by aaron on 5/23/18.
// cartesian_frictional_box.cpp
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

void CartesianFrictionalBox::init(Job* job, Body* body){
    if (job->grid->object_name.compare("CartesianLinear") != 0){
        std::cout << "\nBOUNDARY CONDITION WARNING!" << std::endl;
        std::cout << "\"CartesianFrictionalBox\" boundary expects \"CartesianLinear\" grid NOT \"" << job->grid->object_name << "\"!\n" << std::endl;
    }

    //check that contact properties are set
    if ((body->boundary->fp64_props.size() < 1)){
        //need to coefficient of friction and bodies
        std::cout << body->boundary->fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 1 property defined (mu_f).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        mu_f = body->boundary->fp64_props[0];
        std::cout << "Boundary properties ( mu_f = " << mu_f << ")." << std::endl;
    }

    //find bounds of box
    KinematicVector Lx = KinematicVector(job->JOB_TYPE);
    Lx.setZero();
    for (int i=0; i < body->nodes->x.size(); i++){
        for (int pos=0; pos < body->nodes->x.DIM; pos++){
            if (body->nodes->x(i,pos) > Lx(pos)){
                Lx(pos) = body->nodes->x(i,pos);
            }
        }
    }

    //set bounding mask
    double len = body->nodes->x.size();
    bcNodalMask = KinematicVectorArray(len,job->JOB_TYPE); //should be integer array but its not...
    bcNodalForce = KinematicVectorArray(len,job->JOB_TYPE);
    bcNodalMask.setZero();
    bcNodalForce.setZero();

    bool is_edge = false;
    for (int i=0;i<len;i++){
        is_edge = false;
        for (int pos=0;pos<body->nodes->x.DIM;pos++){
            if (body->nodes->x(i,pos) == 0){
                bcNodalMask(i,pos) = -1; //lower bound
                is_edge = true;
            } else if (body->nodes->x(i,pos) == Lx(pos)) {
                bcNodalMask(i,pos) = 1; //upper bound
                is_edge = true;
            }
        }
        for (int pos=0;pos<body->nodes->x.DIM;pos++){
            if (is_edge && bcNodalMask(i,pos) != 1 && bcNodalMask(i,pos) != -1){
                bcNodalMask(i,pos) = 2;
            }
        }
    }

    std::cout << "Boundary Initialized: [" << body->name << "]." << std::endl;

    return;
}


void CartesianFrictionalBox::generateRules(Job* job, Body* body){
    tmp = body->points->T;
    for (int i=0;i<tmp.size();i++){
        tmp(i) *= body->points->v(i);
    }
    return;
}

void CartesianFrictionalBox::applyRules(Job* job, Body* body){
    double f_normal;
    KinematicVector delta_momentum = KinematicVector(job->JOB_TYPE);

    //map body force to nodes (avoids possible contact forces)
    pvec = body->points->b;
    for (int i=0;i<pvec.size();i++){
        pvec(i) *= body->points->m(i);
    }
    bcNodalForce = body->S * pvec;

    //map divergence of stress
    nvec = body->gradS.left_multiply_by_tensor(tmp);
    for (int i=0;i<nvec.size();i++){
        bcNodalForce(i) -= KinematicVector(nvec[i],bcNodalForce.VECTOR_TYPE);
    }

    for (int i=0;i<body->nodes->x_t.size();i++){
        for (int pos=0;pos<body->nodes->x_t.DIM;pos++){
            if (bcNodalMask(i,pos) == -1){
                //zero out velocity on lower boundary
                //only add friction for closing velocity
                delta_momentum(pos) = -std::min(0.0, body->nodes->mx_t(i,pos) + job->dt * bcNodalForce(i,pos));
                body->nodes->x_t(i,pos) = 0;
                body->nodes->mx_t(i,pos) = 0;
                body->nodes->f(i,pos) = 0;
            } else if (bcNodalMask(i,pos) == 1){
                //zero out velocity on upper boundary
                //only add friction for closing velocity
                delta_momentum(pos) = -std::max(0.0, body->nodes->mx_t(i,pos) + job->dt * bcNodalForce(i,pos));
                body->nodes->x_t(i,pos) = 0;
                body->nodes->mx_t(i,pos) = 0;
                body->nodes->f(i,pos) = 0;
            } else {
                delta_momentum(pos) = 0;
            }
        }
        f_normal = delta_momentum.norm() / job->dt;

        //apply friction to resulting motion
        delta_momentum = body->nodes->mx_t(i) + job->dt * body->nodes->f(i);
        if (delta_momentum.norm() > 0) {
            body->nodes->f(i) -=
                    std::min(delta_momentum.norm() / job->dt, f_normal * mu_f) * delta_momentum / delta_momentum.norm();
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void CartesianFrictionalBox::writeFrame(Job* job, Body* body, Serializer* serializer){
    return;
}

std::string CartesianFrictionalBox::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int CartesianFrictionalBox::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}
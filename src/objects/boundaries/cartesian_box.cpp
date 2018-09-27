//
// Created by aaron on 5/15/18.
// cartesian_box.cpp
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
void CartesianBox::init(Job* job, Body* body){
    if (job->grid->object_name.compare("CartesianLinear") != 0){
        std::cout << "\nBOUNDARY CONDITION WARNING!" << std::endl;
        std::cout << "\"CartesianBox\" boundary expects \"CartesianLinear\" grid NOT \"" << job->grid->object_name << "\"!\n" << std::endl;
    }

    //find bounds of box
    KinematicVector Lx = KinematicVector(job->JOB_TYPE);
    Lx.setZero();
    for (int i=0; i < body->nodes->x.size(); i++){
        for (int pos=0; pos < job->grid->GRID_DIM; pos++){
            if (body->nodes->x(i,pos) > Lx(pos)){
                Lx(pos) = body->nodes->x(i,pos);
            }
        }
    }
    for (int i=job->grid->GRID_DIM; i<Lx.DIM; i++){
        Lx(i) = 0;
    }

    //set bounding mask
    double len = body->nodes->x.size();
    bcNodalMask = KinematicVectorArray(len,job->JOB_TYPE); //should be integer array but its not...
    bcNodalMask.setZero();

    for (int i=0;i<len;i++){
        for (int pos=0;pos<job->grid->GRID_DIM;pos++){
            if (body->nodes->x(i,pos) == 0 || body->nodes->x(i,pos) == Lx(pos)) {
                bcNodalMask(i).setOnes();
                break;
            }
        }
    }

    std::cout << "Boundary Initialized: [" << body->name << "]." << std::endl;

    return;
}


/*----------------------------------------------------------------------------*/
//generate boundary conditions
void CartesianBox::generateRules(Job* job, Body* body){
    //nothing to do here
    return;
}


/*----------------------------------------------------------------------------*/
//apply boundary rules
void CartesianBox::applyRules(Job* job, Body* body){
    for (int i=0;i<body->nodes->x_t.size();i++){
        for (int pos=0; pos < body->nodes->x_t.DIM; pos++){
            if (bcNodalMask(i,pos) == 1){
                //zero out velocity on boundary
                body->nodes->x_t(i,pos) = 0;
                body->nodes->mx_t(i,pos) = 0;
                body->nodes->f(i,pos) = 0;
            }
        }
    }
    return;
}


/*----------------------------------------------------------------------------*/
//write frame and save file
void CartesianBox::writeFrame(Job* job, Body* body, Serializer* serializer){
    serializer->writeVectorArray(bcNodalMask,"bc_nodal_mask");
    return;
}

std::string CartesianBox::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int CartesianBox::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}

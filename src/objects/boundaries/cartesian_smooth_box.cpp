//
// Created by aaron on 5/23/18.
// cartesian_smooth_box.cpp
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
void CartesianSmoothBox::init(Job* job, Body* body){
    if (job->grid->object_name.compare("CartesianLinear") != 0){
        std::cout << "\nBOUNDARY CONDITION WARNING!" << std::endl;
        std::cout << "\"CartesianSmoothBox\" boundary expects \"CartesianLinear\" grid NOT \"" << job->grid->object_name << "\"!\n" << std::endl;
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
    bcNodalMask.setZero();

    for (int i=0;i<len;i++){
        for (int pos=0;pos<body->nodes->x.DIM;pos++){
            if (body->nodes->x(i,pos) == 0 || body->nodes->x(i,pos) == Lx(pos)) {
                bcNodalMask(i,pos) = 1; //only lock normal direction
            }
        }
    }

    std::cout << "Boundary Initialized: [" << body->name << "]." << std::endl;

    return;
}


/*----------------------------------------------------------------------------*/
//generate boundary conditions
void CartesianSmoothBox::generateRules(Job* job, Body* body){
    //nothing to do here
    return;
}


/*----------------------------------------------------------------------------*/
//apply boundary rules
void CartesianSmoothBox::applyRules(Job* job, Body* body){
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
void CartesianSmoothBox::writeFrame(Job* job, Body* body, Serializer* serializer){
    serializer->writeVectorArray(bcNodalMask,"bc_nodal_mask");
    return;
}

std::string CartesianSmoothBox::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int CartesianSmoothBox::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}

//
// Created by aaron on 12/27/19.
// general_custom_boundary.cpp
//

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <cmath>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"
#include "boundaries.hpp"

/* input to this file denotes the type of boundary for each domain face (-x, +x, -y, +y, -z, +z)
 * DO NOT USE WITH AXISYMMETRIC JOBS!
 * -1 -- FREE_BOUNDARY
 * 0 -- NO-SLIP
 */


void GeneralCustomBoundary::init(Job* job, Body* body){
    //check that boundary properties are set
    if ((int_props.size() < 1)){
        //assume all boundaries are no-slip
        limit_props = Eigen::VectorXi(0);               //initialize zero size limit_properties vector
        limit_vals = std::vector<KinematicVector>(0);   //initialize zero size value vector
    } else {
        //assign size of limit properties vector
        limit_props = Eigen::VectorXi(int_props[0]);
        limit_vals = std::vector<KinematicVector>(int_props[0]);

        //define limit props
        for (int i=0;i<limit_props.size();i++){
            limit_props[i] = int_props[i+1];
            if (i < fp64_props.size()/job->grid->GRID_DIM){
                limit_vals[i] = KinematicVector(job->JOB_TYPE);
                for (int ii=0; ii < job->grid->GRID_DIM; ii++){
                    limit_vals[i][ii] = fp64_props[job->grid->GRID_DIM*i + ii];
                }
            } else {
                limit_vals[i] = KinematicVector(job->JOB_TYPE);
                limit_vals[i].setZero();
            }
        }

        std::cout << "Boundary properties (";
        for (int i=0;i<limit_props.size();i++) {
            std::cout << " " << limit_props[i];
        }
        std::cout << ")." << std::endl;
    }

    //initialize bc tags
    bcTag = Eigen::VectorXi(job->grid->node_count);
    bcTag.setConstant(-1);
    bcValues = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    bcValues.setZero();
    for (int i=0; i<bcTag.rows(); i++){
        if (job->grid->nodeTag(job,i) > -1){
            //node is on boundary
            if (job->grid->nodeTag(job,i) < limit_props.rows()) {
                //boundary condition defined by input
                //check that boundary property is valid
                if (limit_props(job->grid->nodeTag(job,i)) <= MAX_PROP_VALUE){
                    bcTag(i) = limit_props(job->grid->nodeTag(job,i));
                    bcValues[i] = limit_vals[job->grid->nodeTag(job,i)];
                } else {
                    bcTag(i) = 0;
                }
            }
        }
    }

    //if bcTag == SLIP_BOUNDARY, then bcValues are outward facing unit normal
    //loop over str-props and assign relevant flags
    std::vector<std::string> options = {"USE_Y_AXIS_ROTATION"};
    for (int i=0; i<str_props.size(); i++){
        switch (Parser::findStringID(options, str_props[i])){
            case 0:
                //USE_Y_AXIS_ROTATION
                USE_Y_AXIS_SYMMETRY = true;
                std::cout << "GeneralCustomBoundary using symmetric y-axis (i.e. y-direction motion only)." << std::endl;
                break;
            default:
                //do nothing
                break;
        }
    }
    if (USE_Y_AXIS_SYMMETRY){
        //loop over points for points on y-axis
        for (int i=0; i<job->grid->node_count; i++){
            if (body->nodes->x(i,0) == 0. && body->nodes->x(i,2) == 0.){
                //node is on y-axis
                bcTag(i) = Y_AXIS_SLIP;
            }
        }
    }

    std::cout << "Boundary Initialized: [" << body->name << "]." << std::endl;

    return;
}

void GeneralCustomBoundary::generateRules(Job* job, Body* body){
    //nothing to do here
    return;
}

void GeneralCustomBoundary::applyRules(Job* job, Body* body){
    for (int i=0;i<body->nodes->x_t.size();i++){
        //assign force that will cancel out force
        if (bcTag(i) == DIRICHLET && body->nodes->m(i) > 0){
            body->nodes->f[i] = (bcValues[i]*body->nodes->m(i) - body->nodes->mx_t[i])/job->dt;
        } else if (bcTag(i) == SLIP_BOUNDARY && body->nodes->m(i) > 0){
            body->nodes->f[i] = -(bcValues[i]*body->nodes->mx_t[i].dot(bcValues[i])/job->dt);
        } else if (bcTag(i) == Y_AXIS_SLIP && body->nodes->m(i) > 0){
            body->nodes->f(i,0) = -body->nodes->mx_t(i,0)/job->dt;
            body->nodes->f(i,2) = -body->nodes->mx_t(i,2)/job->dt;
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void GeneralCustomBoundary::writeFrame(Job* job, Body* body, Serializer* serializer){
    Eigen::VectorXd tmp = bcTag.cast<double>();
    serializer->writeScalarArray(tmp,"bc_tags");
    return;
}

std::string GeneralCustomBoundary::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int GeneralCustomBoundary::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}

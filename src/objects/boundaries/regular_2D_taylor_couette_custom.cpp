//
// Created by aaron on 5/31/18.
// regular_2D_taylor_couette_custom.cpp
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
#include "objects/grids/grids.hpp"

/* input to this file denotes the type of boundary for each domain face (-x, +x, -y, +y, -z, +z)
 * DO NOT USE WITH AXISYMMETRIC JOBS!
 * 0 -- NO-SLIP
 * 1 -- FRICTION-LESS
 * 2 -- FRICTIONAL (need mu_f defined)
 * 3 -- DRIVEN WALL (VELOCITY, need vel defined)
 */

/*----------------------------------------------------------------------------*/
//
void Regular2DTaylorCouetteCustom::init(Job* job, Body* body){
    if (job->grid->object_name.compare("Regular2DTaylorCouetteCell") != 0){
        std::cout << "\nBOUNDARY CONDITION WARNING!" << std::endl;
        std::cout << "\"Regular2DTaylorCouetteCustom\" boundary expects \"Regular2DTaylorCouetteCell\" grid NOT \"" << job->grid->object_name << "\"!\n" << std::endl;
    }

    //check that boundary properties are set
    if ((int_props.size() < 2)){
        //need 2 coefficients for direction
        std::cout << int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined (inner_boundary, outer_boundary).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //define limit props
        inner_int_prop = int_props[0];
        outer_int_prop = int_props[1];

        bool enuf_props = true;

        if (inner_int_prop == FRICTIONAL_WALL || outer_int_prop == FRICTIONAL_WALL){
            if (fp64_props.size() < 2){
                enuf_props = false;
            } else {
                inner_fp64_prop = fp64_props[0];
                outer_fp64_prop = fp64_props[1];
            }
        } else if (inner_int_prop == DRIVEN_VELOCITY || outer_int_prop == DRIVEN_VELOCITY){
            if (fp64_props.size() < 2){
                enuf_props = false;
            } else {
                inner_fp64_prop = fp64_props[0];
                outer_fp64_prop = fp64_props[1];
            }
        }

        if (!enuf_props){
            std::cout << int_props.size() << "\n";
            fprintf(stderr,
                    "%s:%s: Need at least 4 properties defined (inner_boundary, outer_boundary).\n",
                    __FILE__, __func__);
            exit(0);
        }
    }

    //set bounding mask
    int len = body->nodes->x.size();
    bcNodalMask = Eigen::VectorXi(len);
    bcNodalMask.setZero();

    for (int i=0; i<len; i++){
        //-1 inner boundary
        //1 outer boundary
        bcNodalMask(i) = job->grid->nodeTag(job, i);
    }

    //initialize boundary traction measure
    bcReactionForce = KinematicVectorArray(len, job->JOB_TYPE);

    std::cout << "Boundary Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//
void Regular2DTaylorCouetteCustom::generateRules(Job* job, Body* body){
    //setup friction
    tmp = body->points->T;
    for (int i=0;i<tmp.size();i++){
        tmp(i) *= body->points->v(i);
    }
    return;
}

/*----------------------------------------------------------------------------*/
//
void Regular2DTaylorCouetteCustom::applyRules(Job* job, Body* body){
    double f_normal;
    KinematicVector delta_momentum_xyz = KinematicVector(job->JOB_TYPE);
    KinematicVector delta_momentum_rtz = KinematicVector(job->JOB_TYPE);
    KinematicVector tmpVec, tmpPoint;

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

        if (bcNodalMask(i) == -1 || bcNodalMask(i) == 1) {
            //then this is a boundary

            if (bcNodalMask(i) == -1 && inner_int_prop == NO_SLIP_WALL){
                for (int pos=0; pos < body->nodes->x_t.DIM; pos++){
                    //zero out velocity on boundary
                    bcReactionForce(i,pos) = -(body->nodes->mx_t(i,pos)/job->dt + body->nodes->f(i,pos));
                    body->nodes->x_t(i,pos) = 0;
                    body->nodes->mx_t(i,pos) = 0;
                    body->nodes->f(i,pos) = 0;
                }
            } else if (bcNodalMask(i) == 1 && outer_int_prop == NO_SLIP_WALL){
                for (int pos=0; pos < body->nodes->x_t.DIM; pos++){
                    //zero out velocity on boundary
                    bcReactionForce(i,pos) = -(body->nodes->mx_t(i,pos)/job->dt + body->nodes->f(i,pos));
                    body->nodes->x_t(i,pos) = 0;
                    body->nodes->mx_t(i,pos) = 0;
                    body->nodes->f(i,pos) = 0;
                }
            } else {
                tmpPoint = Regular2DTaylorCouetteCell::cPoint_to_rPoint(body->nodes->x(i));

                if (bcNodalMask(i) == -1) {
                    //zero out velocity on inner boundary
                    //only add friction for closing velocity
                    delta_momentum_rtz(0) = -std::min(0.0, Regular2DTaylorCouetteCell::cVector_to_rVector(
                            body->nodes->mx_t(i) + job->dt * bcNodalForce(i), body->nodes->x(i))[0]);

                    tmpVec = Regular2DTaylorCouetteCell::cVector_to_rVector(-(body->nodes->mx_t(i)/job->dt + body->nodes->f(i)), body->nodes->x(i));
                    tmpVec(1) = 0; // don't worry about theta component for now
                    bcReactionForce(i) = Regular2DTaylorCouetteCell::rVector_to_cVector(tmpVec, tmpPoint);

                    tmpVec = Regular2DTaylorCouetteCell::cVector_to_rVector(body->nodes->x_t(i), body->nodes->x(i));
                    tmpVec(0) = 0;
                    body->nodes->x_t(i) = Regular2DTaylorCouetteCell::rVector_to_cVector(tmpVec, tmpPoint);

                    tmpVec = Regular2DTaylorCouetteCell::cVector_to_rVector(body->nodes->mx_t(i), body->nodes->x(i));
                    tmpVec(0) = 0;
                    body->nodes->mx_t(i) = Regular2DTaylorCouetteCell::rVector_to_cVector(tmpVec, tmpPoint);

                    tmpVec = Regular2DTaylorCouetteCell::cVector_to_rVector(body->nodes->f(i), body->nodes->x(i));
                    tmpVec(0) = 0;
                    body->nodes->f(i) = Regular2DTaylorCouetteCell::rVector_to_cVector(tmpVec, tmpPoint);

                    //body->nodes->x_t(i,pos) = 0;
                    //body->nodes->mx_t(i,pos) = 0;
                    //body->nodes->f(i,pos) = 0;
                } else if (bcNodalMask(i) == 1) {
                    //zero out velocity on outer boundary
                    //only add friction for closing velocity
                    delta_momentum_rtz(0) = -std::max(0.0, Regular2DTaylorCouetteCell::cVector_to_rVector(
                            body->nodes->mx_t(i) + job->dt * bcNodalForce(i), body->nodes->x(i))[0]);

                    tmpVec = Regular2DTaylorCouetteCell::cVector_to_rVector(-(body->nodes->mx_t(i)/job->dt + body->nodes->f(i)), body->nodes->x(i));
                    tmpVec(1) = 0; // don't worry about theta component for now
                    bcReactionForce(i) = Regular2DTaylorCouetteCell::rVector_to_cVector(tmpVec, tmpPoint);

                    tmpVec = Regular2DTaylorCouetteCell::cVector_to_rVector(body->nodes->x_t(i), body->nodes->x(i));
                    tmpVec(0) = 0;
                    body->nodes->x_t(i) = Regular2DTaylorCouetteCell::rVector_to_cVector(tmpVec, tmpPoint);

                    tmpVec = Regular2DTaylorCouetteCell::cVector_to_rVector(body->nodes->mx_t(i), body->nodes->x(i));
                    tmpVec(0) = 0;
                    body->nodes->mx_t(i) = Regular2DTaylorCouetteCell::rVector_to_cVector(tmpVec, tmpPoint);

                    tmpVec = Regular2DTaylorCouetteCell::cVector_to_rVector(body->nodes->f(i), body->nodes->x(i));
                    tmpVec(0) = 0;
                    body->nodes->f(i) = Regular2DTaylorCouetteCell::rVector_to_cVector(tmpVec, tmpPoint);
                } else {
                    delta_momentum_rtz(0) = 0;
                }

                f_normal = delta_momentum_rtz.norm() / job->dt;

                //resulting motion should be tangential only...
                if (bcNodalMask(i) == -1 && inner_int_prop == FRICTION_LESS_WALL){
                    //do nothing
                } else if (bcNodalMask(i) == 1 && outer_int_prop == FRICTION_LESS_WALL){
                    //do nothing
                } else if (bcNodalMask(i) == -1 && inner_int_prop == FRICTIONAL_WALL){
                    //add friction to inner wall
                    delta_momentum_xyz = body->nodes->mx_t(i) + job->dt * body->nodes->f(i);
                    if (delta_momentum_xyz.norm() > 0) {
                        body->nodes->f(i) -=
                                std::min(delta_momentum_xyz.norm() / job->dt, f_normal * inner_fp64_prop) * delta_momentum_xyz /
                                delta_momentum_xyz.norm();

                        bcReactionForce(i) -= std::min(delta_momentum_xyz.norm() / job->dt, f_normal * inner_fp64_prop) * delta_momentum_xyz /
                                              delta_momentum_xyz.norm();
                    }
                } else if (bcNodalMask(i) == 1 && outer_int_prop == FRICTIONAL_WALL){
                    //add friction to outer wall
                    delta_momentum_xyz = body->nodes->mx_t(i) + job->dt * body->nodes->f(i);
                    if (delta_momentum_xyz.norm() > 0) {
                        body->nodes->f(i) -=
                                std::min(delta_momentum_xyz.norm() / job->dt, f_normal * outer_fp64_prop) * delta_momentum_xyz /
                                delta_momentum_xyz.norm();

                        bcReactionForce(i) -= std::min(delta_momentum_xyz.norm() / job->dt, f_normal * outer_fp64_prop) * delta_momentum_xyz /
                                              delta_momentum_xyz.norm();
                    }
                } else if (bcNodalMask(i) == -1 && inner_int_prop == DRIVEN_VELOCITY){
                    //add inner wall velocity
                    tmpVec = KinematicVector(job->JOB_TYPE);
                    tmpVec(1) = inner_fp64_prop*body->nodes->m(i); //radial velocity
                    delta_momentum_xyz = Regular2DTaylorCouetteCell::rVector_to_cVector(tmpVec,tmpPoint);
                    delta_momentum_xyz -= body->nodes->mx_t(i) + job->dt * body->nodes->f(i);
                    body->nodes->f(i) += delta_momentum_xyz / job->dt;

                    //only include tangential reaction force for this boundary condition
                    bcReactionForce(i) =  delta_momentum_xyz / job->dt;
                } else if (bcNodalMask(i) == 1 && outer_int_prop == DRIVEN_VELOCITY){
                    //add inner wall velocity
                    tmpVec = KinematicVector(job->JOB_TYPE);
                    tmpVec(1) = outer_fp64_prop*body->nodes->m(i); //radial velocity
                    delta_momentum_xyz = Regular2DTaylorCouetteCell::rVector_to_cVector(tmpVec,tmpPoint);
                    delta_momentum_xyz -= body->nodes->mx_t(i) + job->dt * body->nodes->f(i);
                    body->nodes->f(i) += delta_momentum_xyz / job->dt;

                    //only include tangential reaction force for this boundary condition
                    bcReactionForce(i) =  delta_momentum_xyz / job->dt;
                }


            }
        }
    }
    return;

    /*
    //for now, just implement no-slip
    for (int i=0;i<body->nodes->x_t.size();i++){
        for (int pos=0; pos < body->nodes->x_t.DIM; pos++){
            if (bcNodalMask(i) == 1 || bcNodalMask(i) == -1){
                //zero out velocity on boundary
                body->nodes->x_t(i,pos) = 0;
                body->nodes->mx_t(i,pos) = 0;
                body->nodes->f(i,pos) = 0;
            }
        }
    }
    return;
     */
}

/*----------------------------------------------------------------------------*/
//
void Regular2DTaylorCouetteCustom::writeFrame(Job* job, Body* body, Serializer* serializer){
    Eigen::VectorXd tmpVec = bcNodalMask.cast<double>();
    serializer->writeScalarArray(tmpVec, "bcNodalMask");
    serializer->writeVectorArray(bcReactionForce, "bcReactionForce");

    tmpVec.setOnes();
    KinematicVectorArray gradS = body->gradS.operate(tmpVec,MPMSparseMatrixBase::TRANSPOSED);
    serializer->writeVectorArray(gradS, "gradient");
    return;
}

/*----------------------------------------------------------------------------*/
//
std::string Regular2DTaylorCouetteCustom::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

/*----------------------------------------------------------------------------*/
//
int Regular2DTaylorCouetteCustom::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}
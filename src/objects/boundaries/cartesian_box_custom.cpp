//
// Created by aaron on 5/23/18.
// cartesian_box_custom.cpp
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

/* input to this file denotes the type of boundary for each domain face (-x, +x, -y, +y, -z, +z)
 * DO NOT USE WITH AXISYMMETRIC JOBS!
 * 0 -- NO-SLIP
 * 1 -- FRICTION-LESS
 * 2 -- FRICTIONAL (need mu_f defined)
 * 3 -- PERIODIC
 * 4 -- DRIVEN WALL (VELOCITY, need vel defined)
 */


void CartesianBoxCustom::init(Job* job, Body* body){
    if (job->grid->object_name.compare("CartesianCustom") != 0){
        std::cout << "\nBOUNDARY CONDITION WARNING!" << std::endl;
        std::cout << "\"CartesianBoxCustom\" boundary expects \"CartesianCustom\" grid NOT \"" << job->grid->object_name << "\"!\n" << std::endl;
    }

    //initialize velocity term
    v_set = KinematicVectorArray(2*job->grid->GRID_DIM, job->JOB_TYPE); //one vector for each face

    //check that boundary properties are set
    if ((int_props.size() < 2*job->grid->GRID_DIM)){
        //need 2 coefficients for direction
        std::cout << int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i properties defined ({limit_props}).\n",
                __FILE__, __func__, 2*job->grid->GRID_DIM);
        exit(0);
    } else {
        //define limit props
        limit_props = Eigen::VectorXi(2*job->grid->GRID_DIM);
        for (int i=0;i<limit_props.size();i++){
            limit_props[i] = int_props[i];
        }

        std::vector<int> PROPS_REF_LIST;
        PROPS_REF_LIST.push_back(0); //first position points to start of extra fields
        int tmpINT;
        for (int i=0;i<limit_props.size();i++){
            tmpINT = 0;
            if (limit_props[i] == FRICTIONAL_WALL){
                tmpINT = 1;
            } else if (limit_props[i] == DRIVEN_VELOCITY){
                tmpINT = 3;
            } else if (limit_props[i] == DRIVEN_TRACTION){
                tmpINT = 3;
            }
            tmpINT = PROPS_REF_LIST[i] + tmpINT;
            PROPS_REF_LIST.push_back(tmpINT);
        }

        //store last value:
        tmpINT = PROPS_REF_LIST[PROPS_REF_LIST.size() - 1];

        if ((fp64_props.size() < tmpINT)) {
            std::cout << fp64_props.size() << "\n";
            fprintf(stderr,
                    "%s:%s: Need at least %i properties defined ({limit_props <- OK},{other_props).\n",
                    __FILE__, __func__, tmpINT);
            exit(0);
        }

        std::cout << "Boundary properties (";
        for (int i=0;i<limit_props.size();i++){
            std::cout << " " << limit_props[i];

            if (limit_props[i] == FRICTIONAL_WALL){
                //check for friction and assign coefficient
                mu_f = fp64_props[PROPS_REF_LIST[i]];
                std::cout << "[" << mu_f << "]";
            } else if (limit_props[i] == PERIODIC){
                //check periodic and make sure that both walls are periodic
                if ((i == 0 || i == 1) && (limit_props[0] != PERIODIC || limit_props[1] != PERIODIC)){
                    //this is fatal
                    std::cerr << "ERROR: CartesianBoxCustom given mismatched limits: " << limit_props[0] << " " << limit_props[1] << std::endl;
                    exit(0);
                } else if ((i == 2 || i == 3) && (limit_props[2] != PERIODIC || limit_props[3] != PERIODIC)){
                    //this is fatal
                    std::cerr << "ERROR: CartesianBoxCustom given mismatched limits: " << limit_props[2] << " " << limit_props[3] << std::endl;
                    exit(0);
                } else if ((i == 4 || i == 5) && (limit_props[4] != PERIODIC || limit_props[5] != PERIODIC)){
                    //this is fatal
                    std::cerr << "ERROR: CartesianBoxCustom given mismatched limits: " << limit_props[4] << " " << limit_props[5] << std::endl;
                    exit(0);
                } else {
                    //wooh! you did it right!
                }
            } else if (limit_props[i] == DRIVEN_VELOCITY){
                std::cout << "[";
                for (int pos = 0; pos < v_set.DIM; pos++){
                    v_set(i,pos) = fp64_props[PROPS_REF_LIST[i]+pos];
                    std::cout << " " << v_set(i,pos);
                }
                std::cout << " ]";
            } else if (limit_props[i] == DRIVEN_TRACTION){
                std::cout << "[";
                for (int pos = 0; pos < v_set.DIM; pos++){
                    v_set(i,pos) = fp64_props[PROPS_REF_LIST[i]+pos];
                    std::cout << " " << v_set(i,pos);
                }
                std::cout << " ]";
            }
        }

        std::cout << ")." << std::endl;
    }

    //find bounds of box
    Lx = KinematicVector(job->JOB_TYPE);
    Lx.setZero();
    for (int i=0; i < body->nodes->x.size(); i++){
        for (int pos=0; pos < job->grid->GRID_DIM; pos++){
            if (body->nodes->x(i,pos) > Lx(pos)){
                Lx(pos) = body->nodes->x(i,pos);
            }
        }
    }
    for (int i=job->grid->GRID_DIM; i<body->nodes->x.DIM; i++){
        Lx(i) = 0;
    }

    //set bounding mask
    double len = body->nodes->x.size();
    bcNodalMask = KinematicVectorArray(len, job->JOB_TYPE);
    bcNodalMask.setZero();
    bool is_edge;

    for (int i=0;i<len;i++){
        is_edge = false;
        for (int pos=0;pos<job->grid->GRID_DIM;pos++){
            //check -x, -y, -z faces
            if (body->nodes->x(i,pos) == 0) {
                if (limit_props(2*pos) == FRICTION_LESS_WALL){
                    bcNodalMask(i,pos) = 1;
                } else if (limit_props(2*pos) == NO_SLIP_WALL) {
                    bcNodalMask(i).setOnes();
                } else if (limit_props(2*pos) == FRICTIONAL_WALL){
                    bcNodalMask(i,pos) = -2;
                } else if (limit_props(2*pos) == DRIVEN_VELOCITY && bcNodalMask(i,pos) != 1){
                    bcNodalMask(i).setOnes();
                    bcNodalMask(i,pos) = -4;
                } else if (limit_props(2*pos) == DRIVEN_TRACTION && bcNodalMask(i,pos) != 1){
                    bcNodalMask(i).setOnes();
                    bcNodalMask(i,pos) = -5;
                }
                is_edge = true;
            } else if (body->nodes->x(i,pos) == Lx(pos)) {
                if (limit_props(2*pos + 1) == FRICTION_LESS_WALL){
                    bcNodalMask(i,pos) = 1;
                } else if (limit_props(2*pos + 1) == NO_SLIP_WALL) {
                    bcNodalMask(i).setOnes();
                } else if (limit_props(2*pos + 1) == FRICTIONAL_WALL){
                    bcNodalMask(i,pos) = 2;
                } else if (limit_props(2*pos + 1) == DRIVEN_VELOCITY && bcNodalMask(i,pos) != 1){
                    bcNodalMask(i).setOnes();
                    bcNodalMask(i,pos) = 4;
                } else if (limit_props(2*pos + 1) == DRIVEN_TRACTION && bcNodalMask(i,pos) != 1){
                    bcNodalMask(i).setOnes();
                    bcNodalMask(i,pos) = 5;
                }
                is_edge = true;
            }
        }

        //assign frictional edge
        for (int pos=0;pos<body->nodes->x.DIM;pos++){
            if (is_edge && bcNodalMask(i,pos) != 1 &&
                    bcNodalMask(i,pos) != 2 && bcNodalMask(i,pos) != -2 &&
                    bcNodalMask(i,pos) != 4 && bcNodalMask(i,pos) != -4 &&
                    bcNodalMask(i,pos) != 5 && bcNodalMask(i,pos) != -5){
                bcNodalMask(i,pos) = 3;
            }
        }
    }

    std::cout << "Boundary Initialized: [" << body->name << "]." << std::endl;

    return;
}

void CartesianBoxCustom::generateRules(Job* job, Body* body){
    //setup friction
    tmp = body->points->T;
    for (int i=0;i<tmp.size();i++){
        tmp(i) *= body->points->v(i);
    }

    //setup traction
    if (job->JOB_TYPE == job->JOB_AXISYM) {
        Eigen::VectorXd A_tmp = Eigen::VectorXd(body->points->x.size());
        for (int i = 0; i < body->points->x.size(); i++){
            A_tmp(i) = body->points->v(i) / body->points->x(i,0); //in axisym simulations area of integration equals volume/r
        }
        v_n = body->S * A_tmp;
    } else {
        v_n = body->S * body->points->v;
    }

    //wrap particles about appropriate axes for periodic BCs
    for (int i=0;i<body->points->x.size();i++){
        for (int pos=0;pos<job->grid->GRID_DIM;pos++){
            if (limit_props(2*pos+1) == PERIODIC && body->points->x(i,pos) > Lx(pos)){
                body->points->x(i,pos) -= Lx(pos);
            } else if (limit_props(2*pos) == PERIODIC && body->points->x(i,pos) < 0){
                body->points->x(i,pos) += Lx(pos);
            }
        }
    }
    return;
}

void CartesianBoxCustom::applyRules(Job* job, Body* body){
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
            if (bcNodalMask(i,pos) == 1){
                //zero out boundary and do not adjust friction.
                body->nodes->x_t(i,pos) = 0;
                body->nodes->mx_t(i,pos) = 0;
                body->nodes->f(i,pos) = 0;
            } else if (bcNodalMask(i,pos) == -2){
                //zero out velocity on lower boundary
                //only add friction for closing velocity
                delta_momentum(pos) = -std::min(0.0, body->nodes->mx_t(i,pos) + job->dt * bcNodalForce(i,pos));
                body->nodes->x_t(i,pos) = 0;
                body->nodes->mx_t(i,pos) = 0;
                body->nodes->f(i,pos) = 0;
            } else if (bcNodalMask(i,pos) == 2){
                //zero out velocity on upper boundary
                //only add friction for closing velocity
                delta_momentum(pos) = -std::max(0.0, body->nodes->mx_t(i,pos) + job->dt * bcNodalForce(i,pos));
                body->nodes->x_t(i,pos) = 0;
                body->nodes->mx_t(i,pos) = 0;
                body->nodes->f(i,pos) = 0;
            } else if (bcNodalMask(i,pos) == -4){
                //zero out velocity closing velocity
                //do not add friction
                delta_momentum(pos) = 0; //-std::max(0.0, body->nodes->mx_t(i,pos) + job->dt * bcNodalForce(i,pos));
                body->nodes->f(i) = (body->nodes->m(i) * v_set(2*pos) - body->nodes->mx_t(i)) / job->dt;
                body->nodes->x_t(i,pos) = 0;
                body->nodes->mx_t(i,pos) = 0;
                body->nodes->f(i,pos) = 0;
                break;
            } else if (bcNodalMask(i,pos) == 4){
                //zero out velocity closing velocity
                //do not add friction
                delta_momentum(pos) = 0; //-std::max(0.0, body->nodes->mx_t(i,pos) + job->dt * bcNodalForce(i,pos));
                body->nodes->f(i) = (body->nodes->m(i) * v_set(2*pos + 1) - body->nodes->mx_t(i)) / job->dt;
                body->nodes->x_t(i,pos) = 0;
                body->nodes->mx_t(i,pos) = 0;
                body->nodes->f(i,pos) = 0;
                break;
            } else if (bcNodalMask(i,pos) == -5){
                //zero out velocity closing velocity
                //do not add friction
                delta_momentum(pos) = 0; //-std::max(0.0, body->nodes->mx_t(i,pos) + job->dt * bcNodalForce(i,pos));
                body->nodes->f(i) = v_set(2*pos)*job->grid->nodeSurfaceArea(job,i)*v_n(i)/job->grid->nodeVolume(job,i);
                break;
            } else if (bcNodalMask(i,pos) == 5){
                //zero out velocity closing velocity
                //do not add friction
                delta_momentum(pos) = 0; //-std::max(0.0, body->nodes->mx_t(i,pos) + job->dt * bcNodalForce(i,pos));
                body->nodes->f(i) = v_set(2*pos + 1)*job->grid->nodeSurfaceArea(job,i)*v_n(i)/job->grid->nodeVolume(job,i);
                break;
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

void CartesianBoxCustom::writeFrame(Job* job, Body* body, Serializer* serializer){
    serializer->writeVectorArray(bcNodalMask,"bc_nodal_mask");
    return;
}

std::string CartesianBoxCustom::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int CartesianBoxCustom::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}
//
// Created by aaron on 4/15/20.
// fish.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <dlfcn.h>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "materials.hpp"

/*----------------------------------------------------------------------------*/
//initialize assuming that general properties have been assigned correctly
//fp64_props etc. have been filled by configuration object
void Fish::init(Job* job, Body* body){
    if (fp64_props.size() < 11){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 11 properties defined (E, nu, alpha_h, alpha_a, alpha_t, beta..., gamma...).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        E = fp64_props[0];
        nu = fp64_props[1];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        printf("Material properties (E = %g, nu = %g, G = %g, K = %g).\n",
               E, nu, G, K);

        alpha_h = fp64_props[2];
        alpha_a = fp64_props[3];
        alpha_t = fp64_props[4];

        beta_h = fp64_props[5];
        beta_a = fp64_props[6];
        beta_t = fp64_props[7];

        gamma_h = fp64_props[8];
        gamma_a = fp64_props[9];
        gamma_t = fp64_props[10];

        std::cout << "         Alpha, Beta, Gamma" << std::endl;
        std::cout << "Head:    " << alpha_h << ", " << beta_h << ", " << gamma_h << std::endl;
        std::cout << "Abdomen: " << alpha_a << ", " << beta_a << ", " << gamma_a << std::endl;
        std::cout << "Tail:    " << alpha_t << ", " << beta_t << ", " << gamma_t << std::endl;
    }

    //initialize material matrices
    F = MaterialTensorArray(body->points->x.size());
    F_e = MaterialTensorArray(body->points->x.size());
    F_m = MaterialTensorArray(body->points->x.size());

    //initialize labels for parts of fish
    label = std::vector<int>(body->points->x.size());

    //determine center of mass, and max,min x-position
    double x_min = 0;
    double x_max = 0;
    KinematicVector x_cm = KinematicVector(job->JOB_TYPE);
    double m = 0;
    for (int i=0; i<body->points->x.size(); i++){
        if (i == 0){
            x_min = body->points->x(i,0);
            x_max = body->points->x(i,0);
        } else if (body->points->x(i,0) > x_max){
            x_max = body->points->x(i,0);
        } else if (body->points->x(i,0) < x_min){
            x_min = body->points->x(i,0);
        }

        x_cm += body->points->m(i) * body->points->x[i];
        m += body->points->m(i);
    }
    x_cm /= m;

    //initialize deformation gradients and point labels
    for (int i=0; i<body->points->x.size(); i++){
        //deformation gradients initially identity matrix
        F[i].setIdentity();
        F_e[i].setIdentity();
        F_m[i].setIdentity();

        if (body->points->x(i,1) > x_cm(1)){
            //port side of fish
            if (body->points->x(i,0) > (0.25 * x_min + 0.75 * x_max)){
                //head
                label[i] = PORT_HEAD;
            } else if (body->points->x(i,0) > (0.6 * x_min + 0.4 * x_max)){
                //abdomen
                label[i] = PORT_ABDOMEN;
            } else {
                //tail
                label[i] = PORT_TAIL;
            }
        } else {
            //starboard side of fish
            if (body->points->x(i,0) > (0.25 * x_min + 0.75 * x_max)){
                //head
                label[i] = STARBOARD_HEAD;
            } else if (body->points->x(i,0) > (0.6 * x_min + 0.4 * x_max)){
                //abdomen
                label[i] = STARBOARD_ABDOMEN;
            } else {
                //tail
                label[i] = STARBOARD_TAIL;
            }
        }
    }

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//calculate stress state given prior state and strain-rate
void Fish::calculateStress(Job* job, Body* body, int SPEC){
    //tensors for matrix calculation
    MaterialTensor F_m_inv, B, L;
    F_m_inv = MaterialTensor::Identity();
    double J;

    double lambda_head = std::exp(alpha_h * std::sin(beta_h*job->t + gamma_h));
    double lambda_abdomen = std::exp(alpha_a * std::sin(beta_a*job->t + gamma_a));
    double lambda_tail = std::exp(alpha_t * std::sin(beta_t*job->t + gamma_t));

    for (int i=0;i<body->points->x.size();i++){
        if (body->points->active[i] == 0){
            continue;
        }

        //deformation gradient
        L = body->points->L[i];

        //update F
        F[i] += job->dt*L*F[i];

        //define F_m
        switch (label[i]){
            case PORT_HEAD:
                F_m(i,0,0) = lambda_head;
                break;
            case STARBOARD_HEAD:
                F_m(i,0,0) = 1.0/lambda_head;
                break;
            case PORT_ABDOMEN:
                F_m(i,0,0) = lambda_abdomen;
                break;
            case STARBOARD_ABDOMEN:
                F_m(i,0,0) = 1.0/lambda_abdomen;
                break;
            case PORT_TAIL:
                F_m(i,0,0) = lambda_tail;
                break;
            case STARBOARD_TAIL:
                F_m(i,0,0) = 1.0/lambda_tail;
                break;
            default:
            std::cout << "u";
                break;
        }

        //set F_m_inverse
        F_m_inv(0,0) = 1.0/F_m(i,0,0);

        //determine F_e, J, and left Cauchy Green tensor
        F_e[i] = F[i] * F_m_inv;
        J = F_e[i].det();
        B = F_e[i]*F_e[i].transpose();

        //determine stress (sort of neo-hookean, but simpler to code)
        body->points->T[i] = G*(B - B.trace()/3.0*MaterialTensor::Identity()) + K*log(J)*MaterialTensor::Identity();
    }

    return;
}


/*----------------------------------------------------------------------------*/
//define stress assignement for consistency with history dependent materials
void Fish::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    std::cout << "ERROR: CANNOT ASSIGN STRESS!" << std::endl;
    return;
}


/*----------------------------------------------------------------------------*/
//define pressure assignement for consistency with history dependent materials
void Fish::assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    std::cout << "ERROR: CANNOT ASSIGN STRESS!" << std::endl;
    return;
}


/*----------------------------------------------------------------------------*/
//frame and state writing
void Fish::writeFrame(Job* job, Body* body, Serializer* serializer){
    serializer->writeTensorArray(F, "F");
    serializer->writeTensorArray(F_e, "F_e");
    serializer->writeTensorArray(F_m, "F_m");
    Eigen::VectorXd tmpArray = Eigen::VectorXd(body->points->x.size());
    for (int i=0; i<tmpArray.rows(); i++){
        tmpArray(i) = label[i];
    }
    serializer->writeScalarArray(tmpArray, "fish_label");
    return;
}

std::string Fish::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int Fish::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}
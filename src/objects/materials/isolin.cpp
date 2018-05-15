//
// Created by aaron on 5/15/18.
// isolin.cpp
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
void IsotropicLinearElasticity::init(Job* job, Body* body){
    if (body->material->fp64_props.size() < 2){
        std::cout << body->material->fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined (E, nu).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        E = body->material->fp64_props[0];
        nu = body->material->fp64_props[1];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        printf("Material properties (E = %g, nu = %g, G = %g, K = %g).\n",
               E, nu, G, K);
    }

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//calculate stress state given prior state and strain-rate
void IsotropicLinearElasticity::calculateStress(Job* job, Body* body, int SPEC){
    //tensors for matrix calculation
    MaterialTensor T, L, D, W, tmpMat;

    //double for trace of D
    double trD;

    for (size_t i=0;i<body->points->x.size();i++){
        if (body->points->active[i] == 0){
            continue;
        }

        L = body->points->L[i];
        T = body->points->T[i];

        D = 0.5*(L+L.transpose());
        W = 0.5*(L-L.transpose());

        trD = D.trace();

        //Jaumann objective rate
        tmpMat.setIdentity();
        tmpMat = (2*G*D) + (lambda*trD*MaterialTensor::Identity()) + (W*T) - (T*W);

        body->points->T[i] += job->dt * tmpMat;
    }

    return;
}


/*----------------------------------------------------------------------------*/
//define stress assignement for consistency with history dependent materials
void IsotropicLinearElasticity::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    body->points->T[idIN] = stressIN;
    return;
}


/*----------------------------------------------------------------------------*/
//define pressure assignement for consistency with history dependent materials
void IsotropicLinearElasticity::assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    body->points->T[idIN]
}


/*----------------------------------------------------------------------------*/
//frame and state writing
void IsotropicLinearElasticity::writeFrame(Job* job, Body* body, Serializer* serializer){
    //nothing to report
    return;
}

std::string IsotropicLinearElasticity::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int IsotropicLinearElasticity::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}
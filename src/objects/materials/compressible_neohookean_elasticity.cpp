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
void CompressibleNeohookeanElasticity::init(Job* job, Body* body){

    //call parent initialization
    IsotropicLinearElasticity::init(job, body);

    //size deformation gradient
    F = MaterialTensorArray(body->points->x.size());

    //initialize deformations to identity
    for (int i=0; i<F.size(); i++){
        F[i].setIdentity();
    }

    return;
}

/*----------------------------------------------------------------------------*/
//calculate stress state given prior state and strain-rate
void CompressibleNeohookeanElasticity::calculateStress(Job* job, Body* body, int SPEC){
    //tensors for matrix calculation
    MaterialTensor B, L, tmpMat;
    double J;

    for (int i=0;i<body->points->x.size();i++){
        if (body->points->active[i] == 0){
            continue;
        }

        //velocity gradient
        L = body->points->L[i];

        //evolution of F
        F[i] += L*F[i]*job->dt;

        //calculation of B
        B = F[i]*F[i].transpose();

        //calculation of J
        J = B.det();

        //T = G*B_0 + K*log(J)*I
        body->points->T[i] = G*B.deviator() + K*std::log(J)*MaterialTensor::Identity();
    }

    return;
}


/*----------------------------------------------------------------------------*/
//define stress assignement for consistency with history dependent materials
void CompressibleNeohookeanElasticity::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    std::cout << "WARNING: assignStress() not implemented in CompressibleNeohookeanElasticity!" << std::endl;
    return;
}


/*----------------------------------------------------------------------------*/
//define pressure assignement for consistency with history dependent materials
void CompressibleNeohookeanElasticity::assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    std::cout << "WARNING: assignPressure() not implemented in CompressibleNeohookeanElasticity!" << std::endl;
    return;
}


/*----------------------------------------------------------------------------*/
//frame and state writing
void CompressibleNeohookeanElasticity::writeFrame(Job* job, Body* body, Serializer* serializer){
    serializer->writeTensorArray(F,"F");
    return;
}
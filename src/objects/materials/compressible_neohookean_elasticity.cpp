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
        J = std::sqrt(B.det());

        //T = G*B_0/J + K*log(J)/J*I
        body->points->T[i] = G/J * (B - MaterialTensor::Identity())
                             + (K - 2.0*G/3.0)*std::log(J)/J*MaterialTensor::Identity();
    }

    return;
}


/*----------------------------------------------------------------------------*/
//define stress assignement for consistency with history dependent materials
void CompressibleNeohookeanElasticity::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    if (SPEC != 22766) {
        std::cout << "WARNING: assignStress() not implemented in CompressibleNeohookeanElasticity!" << std::endl;
    } else {
        //do not use this section of code UNLESS using 2D avoid-a-void algorithm
        //forgive me for the messiness...

        //compute volume ratio
        double J = body->points->v(idIN) / body->points->v0(idIN);

        //get s33
        double s33 = stressIN(2,2);

        //recompute J using fixed point iteration
        for (int n = 0; n < 5; n++){
            J = std::exp(J * s33 / (K - 2.0 * G / 3.0));
            //std::cout << idIN << ", " << n << " : " << J << std::endl;
        }

        //compute B
        MaterialTensor B = J / G * (stressIN - s33 * MaterialTensor::Identity())
                           + MaterialTensor::Identity();

        //compute F?
        double detF = std::sqrt(B(0, 0) * B(1, 1) - B(0, 1) * B(0, 1));
        double trF = std::sqrt(B(0, 0) + B(1, 1) + 2.0 * detF);
        double b = B(0, 1) / trF;
        double a = std::sqrt(B(0, 0) - b * b);
        double c = std::sqrt(B(1, 1) - b * b);

        MaterialTensor tmpF = MaterialTensor::Identity();
        tmpF(0, 0) = a;
        tmpF(0, 1) = b;
        tmpF(1, 0) = b;
        tmpF(1, 1) = c;

        //assign F
        F[idIN] = tmpF;

        /*
        if (!std::isfinite(a)) {
            std::cerr << idIN << " : " << "Uh oh! " << detF << " =?= " << J << std::endl;
            exit(0);
        }
        */

        //check that F*F^T = B?
        /*
        MaterialTensor FFt = tmpF * tmpF.transpose();
        std::cout << idIN << " : " << J << " =?= " << detF << std::endl;
        std::cout << idIN << " : " << J*J << " =?= " << B.det() << std::endl;
        std::cout << idIN << " : " << B(0,0) << " =?= " << FFt(0,0) << std::endl;
        std::cout << idIN << " : " << B(0,1) << " =?= " << FFt(0,1) << std::endl;
        std::cout << idIN << " : " << B(1,1) << " =?= " << FFt(1,1) << std::endl;
         */
    }

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

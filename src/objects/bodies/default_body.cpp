//
// Created by aaron on 5/16/18.
// default_body.cpp
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

#include "bodies.hpp"

/*----------------------------------------------------------------------------*/
//
void DefaultBody::init(Job* job){
    int cp = 1;
    for (int i=0;i<job->DIM;i++){
        cp *= 2; //square or cube
    }
    //initialize A matrix
    //for mapping position in cube to id
    //0 -> -1,-1,-1
    //1 -> +1,-1,-1
    //...
    //8 -> +1,+1,+1
    Eigen::VectorXi onoff = -1*Eigen::VectorXi::Ones(job->DIM);
    A = Eigen::MatrixXi(cp, job->DIM);
    for (int c=0; c<cp;c++){
        for (int i=0;i<onoff.rows();i++){
            A(c,i) = onoff(i);
        }
        for (int i=0;i<onoff.rows();i++) {
            if (onoff(i) == -1){
                onoff(i) = 1;
                break;
            } else {
                onoff(i) = -1;
            }
        }
    }

    points->init(job,this);
    nodes->init(job,this);
    material->init(job,this);
    boundary->init(job,this);


    //assign vector type
    S = MPMScalarSparseMatrix(nodes->x.size(), points->x.size());
    gradS = KinematicVectorSparseMatrix(nodes->x.size(), points->x.size(), job->JOB_TYPE);
    std::cout << "Body Initialized: [" << name << "]." << std::endl;
    return;
}


/*----------------------------------------------------------------------------*/
//
std::string DefaultBody::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}


/*----------------------------------------------------------------------------*/
//
int DefaultBody::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//
void DefaultBody::generateMap(Job* job, int SPEC) {
    S.clear();
    gradS.clear();

    //calculate phi and grad phi
    std::vector<int> nvec(0);
    std::vector<double> valvec(0);
    KinematicVectorArray gradvec(0,job->JOB_TYPE);
    KinematicVector tmpGrad(job->JOB_TYPE);
    KinematicVector tmpVec(job->JOB_TYPE);

    int ith_cpdi;

    for (int i = 0; i < points->x.size(); i++) {
        KinematicVector::Map x_i = points->x[i];
        ith_cpdi = SPEC;
        //check whether point is active:
        if (points->active(i) == 0) {
            continue;
        } else if (!job->grid->inDomain(job, x_i)) {
            points->active(i) = 0;
            continue;
        } else if (ith_cpdi == DefaultBody::CPDI_ON) {
            //check that corners are in domain
            for (int c = 0; c < A.rows(); c++) {
                for (int pos = 0; pos < tmpVec.DIM; pos++){
                    tmpVec[pos] = x_i[pos] + 0.5*points->extent(i)*A(c,pos);
                }
                if (!job->grid->inDomain(job, tmpVec)) {
                    //corner out of domain
                    ith_cpdi = DefaultBody::CPDI_OFF;
                    break;
                }
            }
        }

        //calculate map for ith point
        if (ith_cpdi == DefaultBody::CPDI_OFF) {
            nvec.resize(0);
            valvec.resize(0);
            job->grid->evaluateBasisFnValue(job, x_i, nvec, valvec);
            for (int j = 0; j < nvec.size(); j++) {
                S.push_back(nvec[j], i, valvec[j]); //node, point, value
            }

            nvec.resize(0);
            gradvec.resize(0);
            job->grid->evaluateBasisFnGradient(job, x_i, nvec, gradvec);
            for (int j = 0; j < nvec.size(); j++) {
                gradS.push_back(nvec[j], i, gradvec[j]); //node, point, value
            }
        } else if (ith_cpdi == DefaultBody::CPDI_ON) {
            nvec.resize(0);
            valvec.resize(0);
            for (int c = 0; c < A.rows(); c++) {
                //spread influence
                for (int pos = 0; pos < tmpVec.DIM; pos++){
                    tmpVec[pos] = x_i[pos] + 0.5*points->extent(i)*A(c,pos);
                }
                job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
            }
            for (int j = 0; j < nvec.size(); j++) {
                S.push_back(nvec[j], i, valvec[j] / A.rows()); //node, point, value
            }

            //average gradients along sides of extent
            nvec.resize(0);
            valvec.resize(0);
            gradvec.resize(0);
            for (int c = 0; c < A.rows(); c++) {
                //find shape function value at offset point location
                //add node ids to nodevec
                //add values to valvec
                valvec.resize(0);

                for (int pos = 0; pos < tmpVec.DIM; pos++){
                    tmpVec[pos] = x_i[pos] + 0.5*points->extent(i)*A(c,pos);
                }

                job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);

                for (int v = 0; v < valvec.size(); v++) {
                    //gradient contribution from corner
                    //G(x) = (S(x+a) - S(x-a))/(2a)
                    for (int pos = 0; pos < tmpGrad.DIM; pos++){
                        tmpGrad[pos] = valvec[v] / (A.rows() * 0.5 * points->extent(i)) * A(c,pos);
                    }
                    gradvec.push_back(tmpGrad);
                }
            }
            for (int j = 0; j < nvec.size(); j++) {
                gradS.push_back(nvec[j], i, gradvec[j]); //node, point, value
            }
        } else {
            std::cerr << "Unrecognized Argument for use_cpdi in DefaultBody::generateMap(): " << SPEC << std::endl;
        }
    }
    return;
}
//
// Created by aaron on 6/17/19.
// threadpool_points.cpp
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>
#include <sstream>
#include <Eigen/Core>
#include <ctime>

#include "mpm_objects.hpp"
#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "job.hpp"

#include "points.hpp"
#include "objects/bodies/bodies.hpp"

//method to fill in mapping matrix for points between i_begin and i_end
void ThreadPoolPoints::fillMap(Job* job, Body* body, ThreadPoolPoints* points,
             int SPEC,
             MPMScalarSparseMatrix& S,
             KinematicVectorSparseMatrix& gradS,
             int i_begin, int i_end, volatile bool &done){

    S.clear();
    gradS.clear();

    //calculate phi and grad phi
    std::vector<int> nvec(0);
    std::vector<double> valvec(0);
    KinematicVectorArray gradvec(0,job->JOB_TYPE);
    KinematicVector tmpGrad(job->JOB_TYPE);
    KinematicVector tmpVec(job->JOB_TYPE);

    int ith_cpdi;

    for (int i = i_begin; i <= i_end; i++) {
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
            for (int c = 0; c < points->A.rows(); c++) {
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                    tmpVec[pos] = x_i[pos] + 0.5*points->extent(i)*points->A(c,pos);
                }
                for (int pos = job->grid->GRID_DIM; pos < tmpVec.DIM; pos++){
                    tmpVec(pos) = 0;
                }
                if ((points->use_elem && !job->grid->inDomain(job, tmpVec, points->elem(i,c))) || !job->grid->inDomain(job, tmpVec)) {
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
            if (points->use_elem) {
                job->grid->evaluateBasisFnValue(job, x_i, nvec, valvec, points->elem(i,0));
                for (int c=1; c<points->A.rows(); c++){
                    //assign all corners to centroid value
                    points->elem(i,c) = points->elem(i,0);
                }
            } else {
                job->grid->evaluateBasisFnValue(job, x_i, nvec, valvec);
            }
            for (int j = 0; j < nvec.size(); j++) {
                //body->S.push_back(nvec[j], i, valvec[j]); //node, point, value
                S.push_back(nvec[j], i, valvec[j]);
            }

            nvec.resize(0);
            gradvec.resize(0);
            if (points->use_elem) {
                job->grid->evaluateBasisFnGradient(job, x_i, nvec, gradvec, points->elem(i,0));
            } else {
                job->grid->evaluateBasisFnGradient(job, x_i, nvec, gradvec);
            }
            for (int j = 0; j < nvec.size(); j++) {
                //body->gradS.push_back(nvec[j], i, gradvec[j]); //node, point, value
                gradS.push_back(nvec[j], i, gradvec[j]);
            }
        } else if (ith_cpdi == DefaultBody::CPDI_ON) {
            nvec.resize(0);
            valvec.resize(0);
            for (int c = 0; c < points->A.rows(); c++) {
                //spread influence
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                    tmpVec[pos] = x_i[pos] + 0.5*points->extent(i)*points->A(c,pos);
                }
                for (int pos = job->grid->GRID_DIM; pos < tmpVec.DIM; pos++){
                    tmpVec(pos) = 0;
                }
                //job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
                if (points->use_elem) {
                    job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec, points->elem(i,c));
                } else {
                    job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
                }
            }
            for (int j = 0; j < nvec.size(); j++) {
                //body->S.push_back(nvec[j], i, valvec[j] / points->A.rows()); //node, point, value
                S.push_back(nvec[j], i, valvec[j] / points->A.rows());
            }

            //average gradients along sides of extent
            nvec.resize(0);
            valvec.resize(0);
            gradvec.resize(0);
            for (int c = 0; c < points->A.rows(); c++) {
                //find shape function value at offset point location
                //add node ids to nodevec
                //add values to valvec
                valvec.resize(0);

                for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                    tmpVec[pos] = x_i[pos] + 0.5*points->extent(i)*points->A(c,pos);
                }
                for (int pos = job->grid->GRID_DIM; pos < tmpVec.DIM; pos++){
                    tmpVec(pos) = 0;
                }

                //job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
                if (points->use_elem) {
                    job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec, points->elem(i,c));
                } else {
                    job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
                }

                for (int v = 0; v < valvec.size(); v++) {
                    //gradient contribution from corner
                    //G(x) = (S(x+a) - S(x-a))/(2a)
                    for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                        tmpGrad[pos] = valvec[v] / (points->A.rows() * 0.5 * points->extent(i)) * points->A(c,pos);
                    }
                    for (int pos = job->grid->GRID_DIM; pos < tmpGrad.DIM; pos++){
                        tmpGrad[pos] = 0;
                    }
                    gradvec.push_back(tmpGrad);
                }
            }
            for (int j = 0; j < nvec.size(); j++) {
                //body->gradS.push_back(nvec[j], i, gradvec[j]); //node, point, value
                gradS.push_back(nvec[j], i, gradvec[j]);
            }
        } else {
            std::cerr << "Unrecognized Argument for use_cpdi in DefaultBody::generateMap(): " << SPEC << std::endl;
        }
    }

    done = true;

    return;
}

/*---------------------------------------------------------------------------*/

void ThreadPoolPoints::combineMaps(Job* job, Body* body, ThreadPoolPoints* points,
                        MPMScalarSparseMatrix& S,
                        KinematicVectorSparseMatrix& gradS,
                        int i_begin_S, int i_begin_gradS, volatile bool& done){
    //assume that body->S and body->gradS are correctly sized for parallel processing
    for (int i=0; i<S.size(); i++){
        //insert value into container
        body->S.insert(i_begin_S + i, S.i_vec[i], S.j_vec[i], S.buffer[i]);
    }

    for (int i=0; i<gradS.size(); i++){
        //insert value into container
        body->gradS.insert(i_begin_gradS + i, gradS.i_vec[i], gradS.j_vec[i], gradS.buffer[i]);
    }

    done = true;
    return;
}

/*----------------------------------------------------------------------------*/

void ThreadPoolPoints::generateMap(Job* job, Body* body, int SPEC) {
    //the default body will defer to points to generate the mapping
    //body->S.clear();
    //body->gradS.clear();

    //if thread_count <= 1 then run serial code
    if (job->thread_count <= 1){
        DefaultPoints::generateMap(job, body, SPEC);
        return;
    }

    //get maximum index for point calculations
    int i_max = x.size() - 1;

    //determine number of threads to use
    int thread_count = 0;
    if (x.size() >= job->thread_count){
        thread_count = job->thread_count;
    } else {
        thread_count = x.size();
    }

    //check size of memory units
    if (S_vec.size() != thread_count){
        S_vec.resize(thread_count);
    }
    if (gradS_vec.size() != thread_count){
        gradS_vec.resize(thread_count);
    }

    //task completion boolean
    volatile bool firstTaskComplete[thread_count] = {false};

    //choose interval size
    int i_interval = (x.size()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++) {
        //initialize container
        if (gradS_vec[t].VECTOR_TYPE != x.VECTOR_TYPE){
            gradS_vec[t] = KinematicVectorSparseMatrix(0,0,x.VECTOR_TYPE);
        }
        S_vec[t].clear();
        gradS_vec[t].clear();

        //set interval
        i_begin = t * i_interval;
        i_end = i_begin + i_interval - 1;
        if (i_end > i_max){
            i_end = i_max;
        }

        //begin threads
        job->threadPool.doJob(std::bind(fillMap,
                                        job,
                                        body,
                                        this,
                                        SPEC,
                                        std::ref(S_vec[t]),
                                        std::ref(gradS_vec[t]),
                                        i_begin,
                                        i_end,
                                        std::ref(firstTaskComplete[t])));
    }

    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            //std::cout << firstTaskComplete[t] << "?" << std::endl;
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //thread_count is the same, but the intervals need to change based on the size of containers
    volatile bool secondTaskComplete[thread_count] = {false};

    //determine size of S and gradS
    int sizeS = 0;
    int sizegradS = 0;
    for (int t=0; t<thread_count; t++){
        sizeS += S_vec[t].size();
        sizegradS += gradS_vec[t].size();
    }

    //resize body matrices
    if (body->S.size() != sizeS){
        body->S.resize(sizeS);
    }
    if (body->gradS.size() != sizegradS){
        body->gradS.resize(sizegradS);
    }
    //std::cout << sizeS << " ?= " << body->S.size() << std::endl;
    //std::cout << sizegradS << " ?= " << body->gradS.size() << std::endl;
    //exit(0);

    int begin_S = 0;
    int begin_gradS = 0;
    //start threads
    for (int t=0; t<thread_count; t++){
        //spin up job
        job->threadPool.doJob(std::bind(combineMaps,
                                        job,
                                        body,
                                        this,
                                        std::ref(S_vec[t]),
                                        std::ref(gradS_vec[t]),
                                        begin_S,
                                        begin_gradS,
                                        std::ref(secondTaskComplete[t])));
        //update beginning values
        begin_S += S_vec[t].size();
        begin_gradS += gradS_vec[t].size();
    }

    taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    return;
}
//
// Created by aaron on 6/10/19.
// threadpool_explicit_usl.cpp
//

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "solvers.hpp"
#include "objects/bodies/bodies.hpp"
#include "threadpool_explicit_usl.hpp"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <thread>
#include <time.h>
#include <functional>
#include <binders.h>


/*----------------------------------------------------------------------------*/
//
void ThreadPoolExplicitUSL::init(Job* job){
    //check that contact properties are set
    if (int_props.size() == 0){
        //do nothing
        cpdi_spec = DefaultBody::CPDI_ON;
        contact_spec = Contact::IMPLICIT;
        num_threads = job->thread_count;
    } else if (int_props.size() == 1){
        //cpdi_spec given as argument
        cpdi_spec = int_props[0];
        contact_spec = Contact::IMPLICIT;
        num_threads = job->thread_count;
    } else if (int_props.size() == 2){
        //cpdi_spec and contact_spec given
        cpdi_spec = int_props[0];
        contact_spec = int_props[1];
        num_threads = job->thread_count;
    } else if (int_props.size() >= 3){
        //cpdi_spec, contact_spec, num_threads and debug given
        cpdi_spec = int_props[0];
        contact_spec = int_props[1];
        num_threads = job->thread_count;
        if (int_props[2] == 1){
            debug = true;
        }
    }

    //get pointer to job->threadPool
    jobThreadPool = &job->threadPool;

    //initialize memory units
    for (int b=0; b<job->bodies.size(); b++){
        memoryUnits.push_back(parallelMemoryUnit());
    }

    printf("Solver properties (cpdi_spec = %i, contact_spec = %i, num_threads = %i).\n", cpdi_spec, contact_spec, num_threads);
    std::cout << "Solver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
std::string ThreadPoolExplicitUSL::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

/*----------------------------------------------------------------------------*/
//
int ThreadPoolExplicitUSL::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

/*----------------------------------------------------------------------------*/
void ThreadPoolExplicitUSL::mapPointsToNodes(Job* job){
    //check that threadpool exists and not serial
    if (job->thread_count <= 1){
        ExplicitUSL::mapPointsToNodes(job);
        return;
    }

    Body *body;
    Points *points;
    Nodes *nodes;
    Eigen::VectorXd pval;
    Eigen::VectorXd nval;
    KinematicVectorArray pvec;
    MaterialVectorArray nvec;
    MaterialTensorArray tmpMat;
    double tmpVAL;

    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = job->bodies[b].get();
        points = job->bodies[b]->points.get();
        nodes = job->bodies[b]->nodes.get();

        //map mass
        parallelMultiply(body->S, points->m, nodes->m, MPMScalarSparseMatrix::NORMAL, true, b);

        //map momentum
        parallelMultiply(points->x_t, points->m, 1.0, points->mx_t, true);

        parallelMultiply(body->S, points->mx_t, nodes->mx_t, MPMScalarSparseMatrix::NORMAL, true, b);

        parallelDivide(nodes->mx_t, nodes->m, 1.0, nodes->x_t, true);

        //map body force
        pvec = KinematicVectorArray(points->b.size(), points->b.VECTOR_TYPE);
        parallelMultiply(points->b, points->m, 1.0, pvec, true);
        parallelMultiply(body->S, pvec, nodes->f, MPMScalarSparseMatrix::NORMAL, true, b);

        //map divergence of stress
        tmpMat = MaterialTensorArray(points->T.size());
        parallelMultiply(points->T, points->v, 1.0, tmpMat, true);

        nvec = MaterialVectorArray(nodes->x.size());
        parallelMultiply(body->gradS, tmpMat, nvec, MPMSparseMatrixBase::NORMAL, true, b);
        for (int i = 0; i < nvec.size(); i++) {
            nodes->f[i] -= KinematicVector(nvec[i], nodes->f.VECTOR_TYPE);
        }

        if (job->JOB_TYPE == job->JOB_AXISYM){
            pval = Eigen::VectorXd(points->x.size());
            nval = Eigen::VectorXd(nodes->x.size());

            //scale s_tt by r and add contribution to f_r
            for (int i=0;i<pval.rows();i++){
                pval(i) = tmpMat(i,2,2) / points->x(i,0);
            }
            parallelMultiply(body->S, pval, nval, MPMScalarSparseMatrix::NORMAL, true, b);
            for (int i=0; i<nval.rows(); i++){
                nodes->f(i,0) -= nval(i);
            }

            //scale s_rt by r and add contribution to f_t
            for (int i=0; i<pval.rows(); i++){
                pval(i) = tmpMat(i,0,2) / points->x(i,0);
            }
            parallelMultiply(body->S, pval, nval, MPMScalarSparseMatrix::NORMAL, true, b);
            for (int i=0; i<nval.rows(); i++){
                nodes->f(i,2) += nval(i);
            }
        }
        //std::cout << points->m.sum() << " ?= " << nodes->m.sum() << std::endl;
    }
    return;
}

void ThreadPoolExplicitUSL::moveGrid(Job* job){
    //check that threadpool exists and not serial
    if (job->thread_count <= 1){
        ExplicitUSL::moveGrid(job);
        return;
    }

    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //update momentum
        job->bodies[b]->nodes->mx_t += job->dt * job->bodies[b]->nodes->f;

        //calculate velocity
        parallelDivide(job->bodies[b]->nodes->mx_t, job->bodies[b]->nodes->m, 1.0, job->bodies[b]->nodes->x_t, true);

        //set displacement
        job->bodies[b]->nodes->u = job->dt * job->bodies[b]->nodes->x_t;

        //calculate difference in velocity
        parallelDivide(job->bodies[b]->nodes->f, job->bodies[b]->nodes->m, job->dt, job->bodies[b]->nodes->diff_x_t, true);
    }
    return;
}

void ThreadPoolExplicitUSL::movePoints(Job* job){
    //check that threadpool exists and not serial
    if (job->thread_count <= 1){
        ExplicitUSL::movePoints(job);
        return;
    }

    Body* body;
    Points* points;
    Nodes* nodes;
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = job->bodies[b].get();
        points = job->bodies[b]->points.get();
        nodes = job->bodies[b]->nodes.get();

        //map nodal displacement to point positions
        parallelMultiply(body->S, nodes->u, points->x, MPMScalarSparseMatrix::TRANSPOSED, false, b);
        parallelMultiply(body->S, nodes->u, points->u, MPMScalarSparseMatrix::TRANSPOSED, false, b);

        //fix position for out of plane dimension
        if (job->grid->GRID_DIM < job->DIM){
            for (int i=0; i<points->x.size(); i++){
                for (int pos=job->grid->GRID_DIM; pos<job->DIM; pos++){
                    points->x(i,pos) = 0;
                }
            }
        }

        //map nodal velocity diff to points
        parallelMultiply(body->S, nodes->diff_x_t, points->x_t, MPMSparseMatrixBase::TRANSPOSED, false, b);

        //calculate momentum
        parallelMultiply(points->x_t, points->m, 1.0, points->mx_t, true);
    }
    return;
}

void ThreadPoolExplicitUSL::calculateStrainRate(Job* job){
    //check that threadpool exists and not serial
    if (job->thread_count <= 1){
        ExplicitUSL::calculateStrainRate(job);
        return;
    }

    Body* body;
    Points* points;
    Nodes* nodes;
    KinematicVectorArray pvec;

    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = job->bodies[b].get();
        points = job->bodies[b]->points.get();
        nodes = job->bodies[b]->nodes.get();

        //calculate gradient of nodal velocity at points

        if (debug) {
            struct timespec tStart, tFinish;
            clock_gettime(CLOCK_MONOTONIC, &tStart);
            parallelMultiply(body->gradS, nodes->x_t, points->L, MPMSparseMatrixBase::TRANSPOSED, true, b);
            clock_gettime(CLOCK_MONOTONIC, &tFinish);

            std::cout << b << " : " << (tFinish.tv_sec - tStart.tv_sec) +
                                       (tFinish.tv_nsec - tStart.tv_nsec) / 1000000000.0 << std::endl << std::endl;

        } else {
            parallelMultiply(body->gradS, nodes->x_t, points->L, MPMSparseMatrixBase::TRANSPOSED, true, b);
        }

        //correct L if axisymmetric
        if (job->JOB_TYPE == job->JOB_AXISYM){
            pvec = KinematicVectorArray(points->L.size(), points->L.TENSOR_TYPE);
            parallelMultiply(body->S, nodes->x_t, pvec, MPMSparseMatrixBase::TRANSPOSED, true, b);
            for (int i=0; i<points->L.size(); i++){
                points->L(i,0,2) = -pvec(i,2) / points->x(i,0);
                points->L(i,2,2) = pvec(i,0) / points->x(i,0);
            }
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/
void ThreadPoolExplicitUSL::parallelMultiply(const MPMScalarSparseMatrix &S,
                                           const Eigen::VectorXd &x,
                                           Eigen::VectorXd &lhs,
                                           int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = S.size() - 1;

    //get length of output vector
    int i_max = lhs.rows() - 1;

    //determine number of threads for matrix vector mult
    int thread_count;
    if (S.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = S.size();
    }

    //intermediate storage vector
    //std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>> lhs_vec(thread_count);
    if (memoryUnits[memUnitID].s.size() != thread_count){
        memoryUnits[memUnitID].s.resize(thread_count);
    }

    //boolean of completion status
    //std::vector<bool> taskComplete = std::vector<bool>(thread_count);
    /*for (int t=0; t<thread_count; t++){
        taskComplete[t] = false;
    }*/
    volatile bool firstTaskComplete[thread_count] = {false};

    //choose interval size
    int k_interval = (S.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        if (memoryUnits[memUnitID].s[t].rows() != lhs.rows()){
            memoryUnits[memUnitID].s[t] = Eigen::VectorXd(lhs.rows());
        }

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        //threads[t] = std::thread(ssmOperateStoS, std::ref(S), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end);
        //std::function<void (void)> tmpFunc = std::bind(ssmOperateStoSwithFlag, std::ref(S), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end, taskComplete[t]);
        jobThreadPool->doJob(std::bind(ssmOperateStoSwithFlag, std::ref(S), std::ref(x), std::ref(memoryUnits[memUnitID].s[t]), SPEC, k_begin, k_end, std::ref(firstTaskComplete[t])));
        //std::cout << "Started task: " << t << "." << std::endl;
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
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
    //std::cout << "Tasks completed." << std::endl;

    //determine number of threads for addition
    if (lhs.rows() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.rows();
    }

    //boolean for completion
    /*taskComplete = std::vector<bool>(thread_count);
    for (int t=0; t<thread_count; t++){
        taskComplete[t] = false;
    }*/
    volatile bool secondTaskComplete[thread_count] = {false};

    //choose interval size
    int i_interval = (lhs.rows()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        //threads[t] = std::thread(scalarAdd, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear);
        jobThreadPool->doJob(std::bind(scalarAddwithFlag, std::ref(memoryUnits[memUnitID].s), std::ref(lhs), i_begin, i_end, clear, std::ref(secondTaskComplete[t])));
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
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

    //should be all done
    return;
}

void ThreadPoolExplicitUSL::parallelMultiply(const MPMScalarSparseMatrix &S,
                                           const KinematicVectorArray &x,
                                           KinematicVectorArray &lhs,
                                           int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = S.size() - 1;

    //get length of output vector
    int i_max = lhs.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count;
    if (S.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = S.size();
    }

    //intermediate storage vector
    //std::vector<KinematicVectorArray> lhs_vec(thread_count);
    if (memoryUnits[memUnitID].kv.size() != thread_count){
        memoryUnits[memUnitID].kv.resize(thread_count);
    }

    //boolean of completion status
    volatile bool firstTaskComplete[thread_count] = {false};

    //choose interval size
    int k_interval = (S.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        //lhs_vec[t] = KinematicVectorArray(lhs.size(),lhs.VECTOR_TYPE);
        if (memoryUnits[memUnitID].kv[t].size() != lhs.size()){
            memoryUnits[memUnitID].kv[t] = KinematicVectorArray(lhs.size(),lhs.VECTOR_TYPE);
        }

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        //threads[t] = std::thread(ssmOperateVtoV, std::ref(S), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end);
        jobThreadPool->doJob(std::bind(ssmOperateVtoVwithFlag, std::ref(S), std::ref(x), std::ref(memoryUnits[memUnitID].kv[t]), SPEC, k_begin, k_end, std::ref(firstTaskComplete[t])));
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //determine number of threads for addition
    if (lhs.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.size();
    }

    //boolean of completion status
    volatile bool secondTaskComplete[thread_count] = {false};

    //choose interval size
    int i_interval = (lhs.size()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        //threads[t] = std::thread(vectorAddK, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear);
        jobThreadPool->doJob(std::bind(vectorAddKwithFlag, std::ref(memoryUnits[memUnitID].kv), std::ref(lhs), i_begin, i_end, clear, std::ref(secondTaskComplete[t])));
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
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

    //should be all done
    return;
}

void ThreadPoolExplicitUSL::parallelMultiply(const KinematicVectorSparseMatrix &gradS,
                                           const MaterialTensorArray &T,
                                           MaterialVectorArray &lhs,
                                           int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = gradS.size() - 1;

    //get length of output vector
    int i_max = lhs.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count;
    if (gradS.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = gradS.size();
    }

    //intermediate storage vector
    //std::vector<MaterialVectorArray> lhs_vec(thread_count);
    if (memoryUnits[memUnitID].mv.size() != thread_count){
        memoryUnits[memUnitID].mv.resize(thread_count);
    }

    //boolean of completion status
    /*std::vector<bool> taskComplete = std::vector<bool>(thread_count);
    for (int t=0; t<thread_count; t++){
        taskComplete[t] = false;
    }*/
    volatile bool firstTaskComplete[thread_count] = {false};

    //choose interval size
    int k_interval = (gradS.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        //lhs_vec[t] = MaterialVectorArray(lhs.size());
        if (memoryUnits[memUnitID].mv[t].size() != lhs.size()){
            memoryUnits[memUnitID].mv[t] = MaterialVectorArray(lhs.size());
        }

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        //threads[t] = std::thread(kvsmLeftMultiply, std::ref(gradS), std::ref(T), std::ref(lhs_vec[t]), SPEC, k_begin, k_end);
        jobThreadPool->doJob(std::bind(kvsmLeftMultiplywithFlag, std::ref(gradS), std::ref(T), std::ref(memoryUnits[memUnitID].mv[t]), SPEC, k_begin, k_end, std::ref(firstTaskComplete[t])));
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //determine number of threads for addition
    if (lhs.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.size();
    }

    //boolean for task completion
    /*taskComplete = std::vector<bool>(thread_count);
    for (int t=0; t<thread_count; t++){
        taskComplete[t] = false;
    }*/
    volatile bool secondTaskComplete[thread_count] = {false};

    //choose interval size
    int i_interval = (lhs.size()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        //threads[t] = std::thread(vectorAddM, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear);
        jobThreadPool->doJob(std::bind(vectorAddMwithFlag, std::ref(memoryUnits[memUnitID].mv), std::ref(lhs), i_begin, i_end, clear, std::ref(secondTaskComplete[t])));
    }

    //join threads
    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
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

    //should be all done
    return;
}

void ThreadPoolExplicitUSL::parallelMultiply(const KinematicVectorSparseMatrix &gradS,
                                           const KinematicVectorArray &x,
                                           KinematicTensorArray &L,
                                           int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = gradS.size() - 1;

    //get length of output vector
    int i_max = L.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count;
    if (gradS.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = gradS.size();
    }

    //intermediate storage vector
    //std::vector<KinematicTensorArray> lhs_vec(thread_count);
    if (memoryUnits[memUnitID].kt.size() != thread_count){
        memoryUnits[memUnitID].kt.resize(thread_count);
    }

    //boolean array to track task completion
    volatile bool firstTaskComplete[thread_count] = {false};

    //choose interval size
    int k_interval = (gradS.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        //lhs_vec[t] = KinematicTensorArray(L.size(),L.TENSOR_TYPE); //is this slowing things down?
        if (memoryUnits[memUnitID].kt[t].size() != L.size()){
            memoryUnits[memUnitID].kt[t] = KinematicTensorArray(L.size(),L.TENSOR_TYPE);
        }

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(kvsmTensorProductTwithFlag, std::ref(gradS), std::ref(x), std::ref(memoryUnits[memUnitID].kt[t]), SPEC, k_begin, k_end, std::ref(firstTaskComplete[t])));
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //determine number of threads for addition
    if (L.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = L.size();
    }

    //boolean array for second task
    volatile bool secondTaskComplete[thread_count] = {false};

    //choose interval size
    int i_interval = (L.size()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(tensorAddwithFlag, std::ref(memoryUnits[memUnitID].kt), std::ref(L), i_begin, i_end, clear, std::ref(secondTaskComplete[t])));
    }

    //join threads
    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
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

    //should be all done
    return;
}

/*----------------------------------------------------------------------------*/
//methods to manage threads during parallel operations

void ThreadPoolExplicitUSL::parallelMultiply(KinematicVectorArray &kv,
                                           Eigen::VectorXd &s,
                                           double scale,
                                           KinematicVectorArray &out,
                                           bool clear){
    //determine maximum index
    int thread_count;
    int i_max = kv.size() - 1;

    //determine number of threads for addition
    if (kv.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = kv.size();
    }

    //choose interval size
    int i_interval = (kv.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroKVwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        /*for (int t=0; t<thread_count; t++){
            threads[t].join();
        }*/
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //boolean array for second task
    volatile bool secondTaskComplete[thread_count] = {false};

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(multiplyKVbySwithFlag, std::ref(kv), std::ref(s), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
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

void ThreadPoolExplicitUSL::parallelDivide(KinematicVectorArray &kv,
                                         Eigen::VectorXd &s,
                                         double scale,
                                         KinematicVectorArray &out,
                                         bool clear){
    //determine maximum index
    int thread_count;
    int i_max = kv.size() - 1;

    //determine number of threads for addition
    if (kv.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = kv.size();
    }

    //choose interval size
    int i_interval = (kv.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroKVwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //boolean array for second task completion
    volatile bool secondTaskComplete[thread_count] = {false};

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(divideKVbySwithFlag, std::ref(kv), std::ref(s), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
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

void ThreadPoolExplicitUSL::parallelMultiply(Eigen::VectorXd &a,
                                           Eigen::VectorXd &b,
                                           double scale,
                                           Eigen::VectorXd &out,
                                           bool clear){
    //determine maximum index
    int thread_count;
    int i_max = a.size() - 1;

    //determine number of threads for addition
    if (a.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = a.size();
    }

    //choose interval size
    int i_interval = (a.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroSwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //boolean array for second task
    volatile bool secondTaskComplete[thread_count] = {false};

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(multiplySbySwithFlag, std::ref(a), std::ref(b), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
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

void ThreadPoolExplicitUSL::parallelMultiply(MaterialTensorArray &mt,
                                           Eigen::VectorXd &s,
                                           double scale,
                                           MaterialTensorArray &out,
                                           bool clear){
    //determine maximum index
    int thread_count;
    int i_max = mt.size() - 1;

    //determine number of threads for addition
    if (mt.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = mt.size();
    }

    //choose interval size
    int i_interval = (mt.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroMTwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //second boolean array
    volatile bool secondTaskComplete[thread_count] = {false};

    //multipy MaterialTensor by Scalar and add to output
    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(multiplyMTbySwithFlag, std::ref(mt), std::ref(s), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
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

void ThreadPoolExplicitUSL::parallelAdd(Eigen::VectorXd &a,
                                      Eigen::VectorXd &b,
                                      double scale,
                                      Eigen::VectorXd &out,
                                      bool clear){

    //determine maximum index
    int thread_count;
    int i_max = a.size() - 1;

    //determine number of threads for addition
    if (a.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = a.size();
    }

    //choose interval size
    int i_interval = (a.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroSwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //second task completion boolean array
    volatile bool secondTaskComplete[thread_count] = {false};

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(addStoSwithFlag, std::ref(a), std::ref(b), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
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

void ThreadPoolExplicitUSL::parallelSubtract(Eigen::VectorXd &a,
                                           Eigen::VectorXd &b,
                                           double scale,
                                           Eigen::VectorXd &out,
                                           bool clear){
    //determine maximum index
    int thread_count;
    int i_max = a.size() - 1;

    //determine number of threads for addition
    if (a.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = a.size();
    }

    //choose interval size
    int i_interval = (a.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroSwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //second boolean array for task completion
    volatile bool secondTaskComplete[thread_count] = {false};

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(subtractSfromSwithFlag, std::ref(a), std::ref(b), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
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
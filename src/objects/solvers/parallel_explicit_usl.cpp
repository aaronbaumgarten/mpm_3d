//
// Created by aaron on 5/27/19.
// parallel_explicit_usl.cpp
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
#include "parallel_explicit_usl.hpp"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <thread>
#include <time.h>

/*----------------------------------------------------------------------------*/
//
void ParallelExplicitUSL::step(Job* job){
    if (debug) {
        struct timespec timeStart, timeCreateMappings, timeMapPointsToNodes, timeLoads;
        struct timespec timeContacts, timeBCs, timeMoveGrid, timeMovePoints, timeStrainRate, timeDensity, timeGravity;
        struct timespec timeFinish;
        clock_gettime(CLOCK_MONOTONIC, &timeStart);

        //create map
        createMappings(job);
        clock_gettime(CLOCK_MONOTONIC, &timeCreateMappings);

        //map particles to grid
        mapPointsToNodes(job);
        clock_gettime(CLOCK_MONOTONIC, &timeMapPointsToNodes);

        //add arbitrary loading conditions
        generateLoads(job);
        applyLoads(job);
        clock_gettime(CLOCK_MONOTONIC, &timeLoads);

        //add contact forces
        generateContacts(job);
        addContacts(job);
        clock_gettime(CLOCK_MONOTONIC, &timeContacts);

        //enforce boundary conditions
        generateBoundaryConditions(job);
        addBoundaryConditions(job);
        clock_gettime(CLOCK_MONOTONIC, &timeBCs);

        //move grid
        moveGrid(job);
        clock_gettime(CLOCK_MONOTONIC, &timeMoveGrid);

        //move particles
        movePoints(job);
        clock_gettime(CLOCK_MONOTONIC, &timeMovePoints);

        //calculate strainrate
        calculateStrainRate(job);
        clock_gettime(CLOCK_MONOTONIC, &timeStrainRate);

        //update density
        updateDensity(job);
        clock_gettime(CLOCK_MONOTONIC, &timeDensity);

        //add body forces
        job->driver->generateGravity(job);
        job->driver->applyGravity(job);
        clock_gettime(CLOCK_MONOTONIC, &timeGravity);

        //update stress
        updateStress(job);
        clock_gettime(CLOCK_MONOTONIC, &timeFinish);

        std::cout << "createMappings(): " << (timeCreateMappings.tv_sec - timeStart.tv_sec) +
                                             (timeCreateMappings.tv_nsec - timeStart.tv_nsec) / 1000000000.0;
        std::cout << ", mapPointsToNodes(): " << (timeMapPointsToNodes.tv_sec - timeCreateMappings.tv_sec) +
                                                 (timeMapPointsToNodes.tv_nsec - timeCreateMappings.tv_nsec) /
                                                 1000000000.0;
        std::cout << ", applyLoads(): " << (timeLoads.tv_sec - timeMapPointsToNodes.tv_sec) +
                                           (timeLoads.tv_nsec - timeMapPointsToNodes.tv_nsec) / 1000000000.0;
        std::cout << ", addContacts(): " << (timeContacts.tv_sec - timeLoads.tv_sec) +
                                            (timeContacts.tv_nsec - timeLoads.tv_nsec) / 1000000000.0;
        std::cout << ", addBoundaryConditions(): "
                  << (timeBCs.tv_sec - timeContacts.tv_sec) + (timeBCs.tv_nsec - timeContacts.tv_nsec) / 1000000000.0;
        std::cout << ", moveGrid(): "
                  << (timeMoveGrid.tv_sec - timeBCs.tv_sec) + (timeMoveGrid.tv_nsec - timeBCs.tv_nsec) / 1000000000.0;
        std::cout << ", movePoints(): " << (timeMovePoints.tv_sec - timeMoveGrid.tv_sec) +
                                           (timeMovePoints.tv_nsec - timeMoveGrid.tv_nsec) / 1000000000.0;
        std::cout << ", calculateStrainRate(): " << (timeStrainRate.tv_sec - timeMovePoints.tv_sec) +
                                                    (timeStrainRate.tv_nsec - timeMovePoints.tv_nsec) / 1000000000.0;
        std::cout << ", updateDensity(): " << (timeDensity.tv_sec - timeStrainRate.tv_sec) +
                                              (timeDensity.tv_nsec - timeStrainRate.tv_nsec) / 1000000000.0;
        std::cout << ", applyGravity(): " << (timeGravity.tv_sec - timeDensity.tv_sec) +
                                             (timeGravity.tv_nsec - timeDensity.tv_nsec) / 1000000000.0;
        std::cout << ", updateStress(): " << (timeFinish.tv_sec - timeGravity.tv_sec) +
                                             (timeFinish.tv_nsec - timeGravity.tv_nsec) / 1000000000.0;
        std::cout << std::endl;
    } else {
        //run standard job
        ExplicitUSL::step(job);
    }

    return;
}

/*----------------------------------------------------------------------------*/
//
void ParallelExplicitUSL::init(Job* job){
    //check that contact properties are set
    if (int_props.size() == 0){
        //do nothing
        cpdi_spec = DefaultBody::CPDI_ON;
        contact_spec = Contact::IMPLICIT;
        num_threads = 1;
    } else if (int_props.size() == 1){
        //cpdi_spec given as argument
        cpdi_spec = int_props[0];
        contact_spec = Contact::IMPLICIT;
        num_threads = 1;
    } else if (int_props.size() == 2){
        //cpdi_spec and contact_spec given
        cpdi_spec = int_props[0];
        contact_spec = int_props[1];
        num_threads = 1;
    } else if (int_props.size() == 3){
        //cpdi_spec, contact_spec, and num_threads given
        cpdi_spec = int_props[0];
        contact_spec = int_props[1];
        num_threads = int_props[2];
    } else if (int_props.size() >= 4){
        //cpdi_spec, contact_spec, num_threads and debug given
        cpdi_spec = int_props[0];
        contact_spec = int_props[1];
        num_threads = int_props[2];
        if (int_props[3] == 1){
            debug = true;
        }
    }

    //declare threads
    threads = std::vector<std::thread>(num_threads);

    printf("Solver properties (cpdi_spec = %i, contact_spec = %i, num_threads = %i).\n", cpdi_spec, contact_spec, num_threads);
    std::cout << "Solver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
std::string ParallelExplicitUSL::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

/*----------------------------------------------------------------------------*/
//
int ParallelExplicitUSL::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//

void ParallelExplicitUSL::mapPointsToNodes(Job* job){
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
        //nodes->m = body->S * points->m; //m_i = S_ip * m_p
        parallelMultiply(body->S, points->m, nodes->m, MPMScalarSparseMatrix::NORMAL);

        //map momentum
        /*for (int i = 0; i < points->mx_t.size(); i++) {
            points->mx_t(i) = points->m(i) * points->x_t(i);
        }*/
        parallelMultiply(points->x_t, points->m, 1.0, points->mx_t, true);
        //nodes->mx_t = body->S * points->mx_t;
        parallelMultiply(body->S, points->mx_t, nodes->mx_t, MPMScalarSparseMatrix::NORMAL);

        //calculate velocity
        /*
        for (int i = 0; i < nodes->x_t.size(); i++) {
            if (nodes->m(i) > 0) {
                nodes->x_t(i) = nodes->mx_t(i) / nodes->m(i);
            } else {
                nodes->x_t(i).setZero();
            }
        }*/
        parallelDivide(nodes->mx_t, nodes->m, 1.0, nodes->x_t, true);

        //map body force
        pvec = KinematicVectorArray(points->b.size(), points->b.VECTOR_TYPE);
        /*for (int i = 0; i < points->b.size(); i++) {
            pvec(i) = points->m(i) * points->b(i);
        }*/
        parallelMultiply(points->b, points->m, 1.0, pvec, true);
        //nodes->f = body->S * pvec;
        parallelMultiply(body->S, pvec, nodes->f, MPMScalarSparseMatrix::NORMAL);

        //map divergence of stress
        /*
        tmpMat = points->T;
        for (int i = 0; i < tmpMat.size(); i++) {
            tmpMat[i] *= points->v[i];
        }
        //nvec = body->gradS.left_multiply_by_tensor(tmpMat); //f_i = v_p T_p * gradS_ip
         */
        tmpMat = MaterialTensorArray(points->T.size());
        parallelMultiply(points->T, points->v, 1.0, tmpMat, true);

        nvec = MaterialVectorArray(nodes->x.size());
        parallelMultiply(body->gradS, tmpMat, nvec, MPMSparseMatrixBase::NORMAL);
        for (int i = 0; i < nvec.size(); i++) {
            nodes->f[i] -= KinematicVector(nvec[i], nodes->f.VECTOR_TYPE);
        }

        if (job->JOB_TYPE == job->JOB_AXISYM){
            pval = Eigen::VectorXd(points->x.size());

            //scale s_tt by r and add contribution to f_r
            for (int i=0;i<pval.rows();i++){
                pval(i) = tmpMat(i,2,2) / points->x(i,0);
            }
            //nval = body->S * pval;
            parallelMultiply(body->S, pval, nval, MPMScalarSparseMatrix::NORMAL);
            for (int i=0; i<nval.rows(); i++){
                nodes->f(i,0) -= nval(i);
            }

            //scale s_rt by r and add contribution to f_t
            for (int i=0; i<pval.rows(); i++){
                pval(i) = tmpMat(i,0,2) / points->x(i,0);
            }
            //nval = body->S * pval;
            parallelMultiply(body->S, pval, nval, MPMScalarSparseMatrix::NORMAL);
            for (int i=0; i<nval.rows(); i++){
                nodes->f(i,2) += nval(i);
            }
        }
        //std::cout << points->m.sum() << " ?= " << nodes->m.sum() << std::endl;
    }
    return;
}

void ParallelExplicitUSL::moveGrid(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //update momentum
        job->bodies[b]->nodes->mx_t += job->dt * job->bodies[b]->nodes->f;

        //calculate velocity
        /*
        for (int i=0;i<job->bodies[b]->nodes->x_t.size();i++){
            if (job->bodies[b]->nodes->m(i) > 0) {
                job->bodies[b]->nodes->x_t(i) = job->bodies[b]->nodes->mx_t(i) / job->bodies[b]->nodes->m(i);
            } else {
                job->bodies[b]->nodes->x_t(i).setZero();
            }
        }*/
        parallelDivide(job->bodies[b]->nodes->mx_t, job->bodies[b]->nodes->m, 1.0, job->bodies[b]->nodes->x_t, true);

        //set displacement
        job->bodies[b]->nodes->u = job->dt * job->bodies[b]->nodes->x_t;

        //calculate difference in velocity
        /*
        for (int i=0;i<job->bodies[b]->nodes->diff_x_t.size();i++){
            if (job->bodies[b]->nodes->m(i) > 0) {
                job->bodies[b]->nodes->diff_x_t(i) = job->dt * job->bodies[b]->nodes->f(i) / job->bodies[b]->nodes->m(i);
            } else {
                job->bodies[b]->nodes->diff_x_t(i).setZero();
            }
        }*/
        parallelDivide(job->bodies[b]->nodes->f, job->bodies[b]->nodes->m, job->dt, job->bodies[b]->nodes->diff_x_t, true);
    }
    return;
}

void ParallelExplicitUSL::movePoints(Job* job){
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
        //points->x += body->S.operate(nodes->u, MPMSparseMatrixBase::TRANSPOSED);
        //points->u += body->S.operate(nodes->u, MPMSparseMatrixBase::TRANSPOSED);
        parallelMultiply(body->S, nodes->u, points->x, MPMScalarSparseMatrix::TRANSPOSED, false);
        parallelMultiply(body->S, nodes->u, points->u, MPMScalarSparseMatrix::TRANSPOSED, false);

        //fix position for out of plane dimension
        if (job->grid->GRID_DIM < job->DIM){
            for (int i=0; i<points->x.size(); i++){
                for (int pos=job->grid->GRID_DIM; pos<job->DIM; pos++){
                    points->x(i,pos) = 0;
                }
            }
        }

        //map nodal velocity diff to points
        //points->x_t += body->S.operate(nodes->diff_x_t, MPMSparseMatrixBase::TRANSPOSED);
        parallelMultiply(body->S, nodes->diff_x_t, points->x_t, MPMSparseMatrixBase::TRANSPOSED, false);

        //calculate momentum
        /*for (int i=0;i<points->mx_t.size();i++){
            points->mx_t(i) = points->m(i) * points->x_t(i);
        }*/
        parallelMultiply(points->x_t, points->m, 1.0, points->mx_t, true);
    }
    return;
}

void ParallelExplicitUSL::calculateStrainRate(Job* job){
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
        //body->bodyCalcPointGradient(job,points->L,nodes->x_t,Body::SET);
        //points->L = body->gradS.tensor_product_transpose(nodes->x_t, MPMSparseMatrixBase::TRANSPOSED);
        parallelMultiply(body->gradS, nodes->x_t, points->L, MPMSparseMatrixBase::TRANSPOSED);

        //correct L if axisymmetric
        if (job->JOB_TYPE == job->JOB_AXISYM){
            //pvec = body->S.operate(nodes->x_t, MPMSparseMatrixBase::TRANSPOSED);
            pvec = KinematicVectorArray(points->L.size(), points->L.TENSOR_TYPE);
            parallelMultiply(body->S, nodes->x_t, pvec, MPMSparseMatrixBase::TRANSPOSED);
            for (int i=0; i<points->L.size(); i++){
                points->L(i,0,2) = -pvec(i,2) / points->x(i,0);
                points->L(i,2,2) = pvec(i,0) / points->x(i,0);
            }
        }
    }
    return;
}

void ParallelExplicitUSL::updateDensity(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (int i=0;i<job->bodies[b]->points->v.rows();i++) {
            job->bodies[b]->points->v(i) *= std::exp(job->dt * job->bodies[b]->points->L(i).trace());
        }

        //this is new, but maybe useful
        job->bodies[b]->points->updateIntegrators(job,job->bodies[b].get());
    }
    return;
}
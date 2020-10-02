//
// Created by aaron on 12/23/19.
// fvm_default_solver.cpp
//


#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <fstream>
#include <job.hpp>
#include <time.h>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "mpm_objects.hpp"
#include "fvm_objects.hpp"
#include "fvm_solvers.hpp"
#include "objects/bodies/bodies.hpp"
#include "objects/solvers/solvers.hpp"

void ParallelMixtureSolverRK4::init(Job* job, FiniteVolumeDriver* driver){

    //call parent initialization f'ns
    FVMMixtureSolverRK4::init(job, driver);
    ThreadPoolExplicitUSL::init(job);

    return;
}

/*----------------------------------------------------------------------------*/
//one forward mixed time-stpe

void ParallelMixtureSolverRK4::step(Job* job, FiniteVolumeDriver* driver){

    //tell mixture solver base class to take step
    FVMMixtureSolverRK4::step(job, driver);

    return;
}

/*---------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ParallelMixtureSolverRK4::mapPointsToNodes(Job* job){
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

void ParallelMixtureSolverRK4::moveGrid(Job* job){
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

void ParallelMixtureSolverRK4::movePoints(Job* job){
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

void ParallelMixtureSolverRK4::calculateStrainRate(Job* job){
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
        parallelMultiply(body->gradS, nodes->x_t, points->L, MPMSparseMatrixBase::TRANSPOSED, true, b);

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
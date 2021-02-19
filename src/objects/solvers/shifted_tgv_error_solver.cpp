//
// Created by aaron on 5/16/18.
// explicit_usl.cpp
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

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

KinematicVector ShiftedTGVErrorSolver::getVelocity(Job* job, KinematicVector const &x){
    KinematicVector result = KinematicVector(job->JOB_TYPE);
    result[0] = u_max * std::sin(a*x[0]) * std::cos(a*x[1]);
    result[1] = -u_max * std::cos(a*x[0]) * std::sin(a*x[1]);
    return result;
}

KinematicVector ShiftedTGVErrorSolver::getAcceleration(Job* job, KinematicVector const &x){
    KinematicVector result = KinematicVector(job->JOB_TYPE);
    result[0] = u_max * u_max * a * std::sin(a*x[0]) * std::cos(a*x[0]);
    result[1] = u_max * u_max * a * std::sin(a*x[1]) * std::cos(a*x[1]);
    return result;
}


double ShiftedTGVErrorSolver::getPressure(Job* job, KinematicVector const &x){
    //this time, with gravity!
    return density * u_max * u_max * 0.25 * (std::cos(2.0*a*x[0]) + std::cos(2.0*a*x[1])) - density * 9.81 * x[1];
}

void ShiftedTGVErrorSolver::calculateAcceleration(Job* job){
    Body *body;
    Points *points;
    Nodes *nodes;
    Eigen::VectorXd pval;
    Eigen::VectorXd nval;
    KinematicVectorArray pvec;
    MaterialVectorArray nvec;
    MaterialTensorArray tmpMat;
    KinematicVector g = KinematicVector(job->JOB_TYPE);
    g.setZero(); g[1] = -9.81;
    double tmpVAL;

    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = job->bodies[b].get();
        points = job->bodies[b]->points.get();
        nodes = job->bodies[b]->nodes.get();

        //map mass
        nodes->m = body->S * points->m; //m_i = S_ip * m_p

        //map divergence of stress
        nodes->f.setZero();
        tmpMat = points->T;
        for (int i = 0; i < tmpMat.size(); i++) {
            tmpMat[i] *= points->v[i];
        }
        nvec = body->gradS.left_multiply_by_tensor(tmpMat); //f_i = v_p T_p * gradS_ip
        for (int i = 0; i < nvec.size(); i++) {
            nodes->f[i] -= KinematicVector(nvec[i], nodes->f.VECTOR_TYPE);
        }

        //add gravity
        pvec = KinematicVectorArray(points->x.size(), points->x.VECTOR_TYPE);
        for (int i = 0; i < points->x.size(); i++) {
            pvec(i) = points->m(i) * g;
        }
        nodes->f += body->S * pvec;

        //zero out forces at +/- y,z boundaries (slow unfortunately, but useful for simplicity)
        KinematicVector x_min = nodes->x[0];
        KinematicVector x_max = x_min;
        for (int i=1; i<nodes->x.size(); i++){
            for (int dir=0; dir<job->grid->GRID_DIM; dir++){
                if (nodes->x(i,dir) > x_max[dir]){
                    x_max[dir] = nodes->x(i,dir);
                } else if (nodes->x(i,dir) < x_min[dir]){
                    x_min[dir] = nodes->x(i,dir);
                }
            }
        }
        for (int i=0; i<nodes->x.size(); i++){
            for (int dir=0; dir<job->grid->GRID_DIM; dir++){
                if (nodes->x(i,dir) == x_min[dir] || nodes->x(i,dir) == x_max[dir]){
                    //set dir velocity to zero
                    nodes->f(i,dir) = 0;
                    a_M(i,dir) = 0; //boundary condition
                }
            }
        }
    }
    return;
}

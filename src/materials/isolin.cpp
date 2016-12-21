//
// Created by aaron on 9/8/16.
// isolin.cpp
//

/**
    \file isolin.c
    \author Sachith Dunatunga
    \date 23.10.13

    The isotropic linear elastic material model.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "process.hpp"
#include "tensor.hpp"

#define MAT_VERSION_STRING "1.0" __DATE__ " " __TIME__

/* timestep passed as argument */
double dt = 1e-6;

/* Material Constants (set by configuration file). */
double E, nu, G, K;
double lambda;

void material_init(Body *body);

void calculate_stress_threaded(threadtask_t *task, Body *body, double dt);

void calculate_stress(Body *body, double dt);

#define SZZ_STATE 0

/*----------------------------------------------------------------------------*/
void material_init(Body *body) {
    size_t i, j;

    for (i = 0; i < body->p; i++) {
        for (j = 0; j < DEPVAR; j++) {
            body->particles.state(i,j) = 0;
        }
    }

    if (body->material.num_fp64_props < 2) {
        std::cout << body->material.num_fp64_props << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined (E, nu).\n",
                __FILE__, __func__);
        exit(0); //replace error code handler from sachith's work with 0
    } else {
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        printf("%s:%s: properties (E = %g, nu = %g, G = %g, K = %g).\n",
               __FILE__, __func__, E, nu, G, K);
    }

    printf("%s:%s: (material version %s) done initializing material.\n",
           __FILE__, __func__, MAT_VERSION_STRING);
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress(Body *body, double dtIn) {
    dt = dtIn;
    threadtask_t t;
    t.offset = 0;
    t.blocksize = body->p;
    calculate_stress_threaded(&t,body, dtIn);
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_implicit(Body *body, double dtIn) {

    dt = dtIn;

    for (size_t i = 0; i < body->p; i++) {
        if (body->particles.active[i] == 0) {
            continue;
        }

        Eigen::VectorXd tmpVec(9);

        tmpVec << body->particles.L.row(i).transpose();
        Eigen::Matrix3d L(tmpVec.data());

        tmpVec << body->particles.T.row(i).transpose();
        Eigen::Matrix3d T(tmpVec.data());

        Eigen::Matrix3d D = 0.5*(L+L.transpose());
        Eigen::Matrix3d W = 0.5*(L-L.transpose());

        double trD = D.trace();

        Eigen::Matrix3d gleft = W*T;
        Eigen::Matrix3d gright = T*W;

        Eigen::Matrix3d tmpMat = gleft*gright;

        Eigen::Matrix3d CD = 2*G*D + lambda*trD*Eigen::Matrix3d::Identity();

        Eigen::Matrix3d dsj = CD + tmpMat;

        for (size_t pos=0;pos<9;pos++){
            body->particles.Ttrial(i,pos) = body->particles.T(i,pos) + dt*dsj(pos);
        }
    }

}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_threaded(threadtask_t *task, Body *body, double dtIn) {
    //printf("t = %g\nx = %g\ny = %g\n", job->t, job->particles[0].x, job->particles[0].y);
    /*printf("L: \n%g %g %g\n%g %g %g\n%g %g %g\n\n",
           job->particles[0].L[0],
           job->particles[0].L[1],
           job->particles[0].L[2],
           job->particles[0].L[3],
           job->particles[0].L[4],
           job->particles[0].L[5],
           job->particles[0].L[6],
           job->particles[0].L[7],
           job->particles[0].L[8]);*/
    /* Since this is local, we can split the particles among the threads. */
    dt = dtIn;
    size_t p_start = task->offset;
    size_t p_stop = task->offset + task->blocksize;

    for (size_t i = 0; i < body->p; i++) {
        if (body->particles.active[i] == 0) {
            continue;
        }

        Eigen::VectorXd tmpVec(9);

        tmpVec << body->particles.L.row(i).transpose();
        Eigen::Matrix3d L(tmpVec.data());

        tmpVec << body->particles.T.row(i).transpose();
        Eigen::Matrix3d T(tmpVec.data());

        Eigen::Matrix3d D = 0.5*(L+L.transpose());
        Eigen::Matrix3d W = 0.5*(L-L.transpose());

        double trD = D.trace();

        Eigen::Matrix3d gleft = W*T;
        Eigen::Matrix3d gright = T*W;

        Eigen::Matrix3d tmp = gleft*gright;

        Eigen::Matrix3d CD = 2*G*D + lambda*trD*Eigen::Matrix3d::Identity();

        Eigen::Matrix3d dsj = CD + tmp;

        for (size_t pos=0;pos<9;pos++){
            body->particles.T(i,pos) += dt*dsj(pos);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/


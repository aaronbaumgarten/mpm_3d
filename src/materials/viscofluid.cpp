//
// Created by aaron on 2/6/17.
// viscofluid.cpp
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "process.hpp"
#include "tensor.hpp"
#include <Eigen/Dense>

#define MAT_VERSION_STRING "1.0" __DATE__ " " __TIME__

/* timestep passed as argument */
double dt = 1e-6;

/* Material Constants (set by configuration file). */
double mu, K;
double lambda;

extern "C" void material_init(Body *body);

extern "C" void calculate_stress_threaded(threadtask_t *task, Body *body, double dt);

extern "C" void calculate_stress(Body *body, double dt);

extern "C" void calculate_stress_implicit(Body *body, double dt);

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
                "%s:%s: Need at least 2 properties defined (K, mu).\n",
                __FILE__, __func__);
        exit(0); //replace error code handler from sachith's work with 0
    } else {
        K = body->material.fp64_props[0];
        mu = body->material.fp64_props[1];
        printf("Material properties (K = %g, mu = %g).\n",
               K,mu);
    }

    std::cout << "Done initializing material (" << body->id << ").\n";
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
        //Eigen::Matrix3d L(tmpVec.data());
        Eigen::Matrix<double, 3, 3, Eigen::RowMajor> L(tmpVec.data());

        tmpVec << body->particles.T.row(i).transpose();
        //Eigen::Matrix3d T(tmpVec.data());
        Eigen::Matrix<double, 3, 3, Eigen::RowMajor> T(tmpVec.data());

        Eigen::Matrix3d D = 0.5*(L+L.transpose());

        double trT = T.trace();
        double trD = D.trace();

        //assume Je = J and F(n+1) ~ (1+L*dt)F(n)
        double dJ = (dt*L + Eigen::Matrix3d::Identity()).determinant();

        Eigen::Matrix3d Sv = 2 * mu * (D - trD/3.0 * Eigen::Matrix3d::Identity());
        Eigen::Matrix3d Se = (trT/3.0 + K*std::log(dJ))*Eigen::Matrix3d::Identity();

        for (size_t pos=0;pos<9;pos++){
            body->particles.Ttrial(i,pos) = Se(pos) + Sv(pos);
        }
    }

    return;

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

        body->material.calculate_stress_implicit(body,dtIn);

        for (size_t pos=0;pos<9;pos++){
            body->particles.T(i,pos) = body->particles.Ttrial(i,pos);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/


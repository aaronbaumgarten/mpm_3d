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
            body->particles[i].state[j] = 0;
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
    calculate_stress_threaded(&t,body);
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_threaded(threadtask_t *task, Body *body) {
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
    size_t p_start = task->offset;
    size_t p_stop = task->offset + task->blocksize;

    for (size_t i = p_start; i < p_stop; i++) {
        if (body->particle_active[i] == 0) {
            continue;
        }

        double e_t[9];
        double w_t[9];
        tensor_sym3(e_t,body->particles[i].L);
        tensor_skw3(w_t,body->particles[i].L);

        double trD;
        tensor_trace3(&trD,body->particles[i].L);

        double dsj[9]; //this formulation may be wrong. Second two terms are confusing
        dsj[XX] = lambda * trD + 2.0 * G * e_t[XX]
                  + 2.0 * w_t[XY] * body->particles[i].T[XY]
                  - 2.0 * w_t[ZX] * body->particles[i].T[XZ];

        dsj[YY] = lambda * trD + 2.0 * G * e_t[YY]
                  - 2.0 * w_t[XY] * body->particles[i].T[XY]
                  + 2.0 * w_t[YZ] * body->particles[i].T[YZ];

        dsj[ZZ] = lambda * trD + 2.0 * G * e_t[ZZ]
                  - 2.0 * w_t[YZ] * body->particles[i].T[YZ]
                  + 2.0 * w_t[ZX] * body->particles[i].T[XZ];

        dsj[XY] = 2.0 * G * e_t[XY]
                  - w_t[XY] * (body->particles[i].T[XX] - body->particles[i].T[YY]);

        dsj[XZ] = 2.0 * G * e_t[XZ]
                  - w_t[ZX] * (body->particles[i].T[ZZ] - body->particles[i].T[XX]);

        dsj[YZ] = 2.0 * G * e_t[YZ]
                  - w_t[YZ] * (body->particles[i].T[YY] - body->particles[i].T[ZZ]);

        body->particles[i].T[XX] += dt * dsj[XX];
        body->particles[i].T[XY] += dt * dsj[XY];
        body->particles[i].T[XZ] += dt * dsj[XZ];
        body->particles[i].T[YX] += dt * dsj[XY];//*
        body->particles[i].T[YY] += dt * dsj[YY];
        body->particles[i].T[YZ] += dt * dsj[YZ];
        body->particles[i].T[ZX] += dt * dsj[XZ];//*
        body->particles[i].T[ZY] += dt * dsj[YZ];//*
        body->particles[i].T[ZZ] += dt * dsj[ZZ];

        /*const double exx_t = body->particles[i].L[XX];
        const double exy_t = 0.5 * (job->particles[i].L[XY] + job->particles[i].L[YX]);
        const double wxy_t = 0.5 * (job->particles[i].L[XY] - job->particles[i].L[YX]);
        const double eyy_t = job->particles[i].L[YY];

        const double trD = exx_t + eyy_t;
        double dsjxx = lambda * trD + 2.0 * G * exx_t;
        double dsjxy = 2.0 * G * exy_t;
        double dsjyy = lambda * trD + 2.0 * G * eyy_t;
        dsjxx += 2 * wxy_t * job->particles[i].T[XY];
        dsjxy -= wxy_t * (job->particles[i].T[XX] - job->particles[i].T[YY]);
        dsjyy -= 2 * wxy_t * job->particles[i].T[XY];
        const double dsjzz = lambda * trD;

        job->particles[i].T[XX] += job->dt * dsjxx;
        job->particles[i].T[XY] += job->dt * dsjxy;
        job->particles[i].T[YY] += job->dt * dsjyy;
        job->particles[i].T[ZZ] += job->dt * dsjzz;*/
    }

    return;
}
/*----------------------------------------------------------------------------*/
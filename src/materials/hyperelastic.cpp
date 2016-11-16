//
// Created by aaron on 11/4/16.
// hyperelastic.cpp
//

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
        for (j = 0; j < 9; j++){
            body->particles[i].F[j] = 0;
        }

        body->particles[i].state[XX] = 1;
        body->particles[i].state[YY] = 1;
        body->particles[i].state[ZZ] = 1;

        body->particles[i].F[XX] = 1;
        body->particles[i].F[YY] = 1;
        body->particles[i].F[ZZ] = 1;
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
    dt = dtIn; //this is dangerous
    threadtask_t t;
    t.offset = 0;
    t.blocksize = body->p;
    calculate_stress_threaded(&t,body, dtIn);
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

    for (size_t i = p_start; i < p_stop; i++) {
        if (body->particles[i].active[0] == 0) {
            continue;
        }

        double f[9]; //finite deformation
        double fbar[9]; //isochoric finite deformation
        double fbarT[9];
        double bebar[9]; //left cauchy green tensor (stored as state)
        double s[9]; //deviatoric stress
        double tmp[9]; //temporary tensor
        double one[9] = {1,0,0,0,1,0,0,0,1};

        double J; //detF
        double detf;
        double p; //pressure

        tensor_copy3(tmp,body->particles[i].L);
        tensor_scale3(tmp,dt);
        tensor_add3(f,one,tmp); //calculate f

        tensor_copy3(tmp,body->particles[i].F);
        tensor_multiply3(body->particles[i].F,f,tmp); //update F

        tensor_det3(&detf,f);
        tensor_copy3(fbar,f);
        tensor_scale3(fbar,1/cbrt(detf)); //calculate f isochoric

        tensor_transpose3(fbarT,fbar);
        tensor_multiply3(tmp,body->particles[i].state,fbarT);
        tensor_multiply3(bebar,fbar,tmp); //calculate be isochoric

        tensor_dev3(s,bebar);
        tensor_scale3(s,G); //calculate deviatoric stress

        tensor_det3(&J,body->particles[i].F);
        p = K/2.0 * (J*J-1) / J;
        tensor_copy3(tmp,one);
        tensor_scale3(tmp,J*p);
        tensor_add3(body->particles[i].T,tmp,s); //cauchy stress

        tensor_copy3(body->particles[i].state,bebar); //isochoric left cauchy-green tensor

    }

    return;
}
/*----------------------------------------------------------------------------*/
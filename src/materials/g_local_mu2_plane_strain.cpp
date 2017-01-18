//
// Created by aaron on 1/18/17.
// g_local_mu2_plane_strain.cpp
//

/**
    \file g_local_mu2.c
    \author Sachith Dunatunga
    \date 04.12.13

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <Eigen/Core>
#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "process.hpp"
#include "tensor.hpp"

/* Estimated (see Kamrin and Koval 2014 paper) */
#define MU_S 0.280
#define GRAINS_RHO 2450
#define RHO_CRITICAL 1500
#define MU_2 (I_0+MU_S)
#define DELTA_MU (MU_2 - MU_S)
#define I_0 0.278

#define ELEMENT_DISTANCE_STATE 4
#define GAMMAP_STATE 9
#define GAMMADOTP_STATE 10

#define MAT_VERSION_STRING "1.0 " __DATE__ " " __TIME__

extern "C" void material_init(Body *body);

extern "C" void calculate_stress_threaded(threadtask_t *task, Body *body, double dt);

extern "C" void calculate_stress(Body *body, double dt);

extern "C" void calculate_stress_implicit(Body *body, double dt);

/*
    The Young's Modulus (E) and Poisson ratio (nu), set in the material init
    procedure. These are Ccpies of the properties given in the
    configuration file. Shear modulus (G) and bulk modulus (K) are derived from
    these, as is Lame Coefficient (lambda).
*/
static double E;
static double nu;
static double G;
static double K;
static double lambda;
static double grains_d;

void quadratic_roots(double *x1, double *x2, double a, double b, double c)
{
    if (a == 0) {
        /* not a quadratic... */
        if (b == 0) {
            /* what */
            *x1 = 0;
            *x2 = 0;
        } else {
            *x1 = c / b;
            *x2 = 0;
        }
    } else {
        *x1 = (-b - copysign(sqrt(b * b - 4 * a * c), b)) / (2 * a);
        *x2 = c / (a * *x1);
    }
    return;
}

double negative_root(double a, double b, double c)
{
    double x;
    if (b > 0) {
        x = (-b - sqrt(b*b - 4*a*c)) / (2*a);
    } else {
        x = (2*c) / (-b + sqrt(b*b - 4*a*c));
    }
    return x;
}

/*----------------------------------------------------------------------------*/
void material_init(Body *body)
{
    for (size_t i = 0; i < body->p; i++) {
        for (size_t j = 0; j < DEPVAR; j++) {
            body->particles.state(i,j) = 0;
        }
    }

    if (job->material.num_fp64_props < 3) {
        // Bit of a hack, but it's okay for now. Just close your eyes and code it anyways.
        fprintf(stderr,
                "%s:%s: Need at least 3 properties defined (E, nu, grain diameter).\n",
                __FILE__, __func__);
        exit(EXIT_FAILURE);
    } else {
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        grains_d = body->material.fp64_props[2];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2*nu));
        lambda = K - 2.0 * G / 3.0;

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, grain diameter = %g).\n",
               E, nu, G, K, grains_d);
    }

    double total_particle_mass = body->particles.m.sum();

    double average_particle_mass = total_particle_mass / job->num_particles;

    for (size_t i = 0; i < body->p; i++) {
        body->particles.state(i, GAMMAP_STATE) = 0;
        body->particles.state(i, GAMMADOTP_STATE) = 0;
    }

    std::cout << "Done initializing material (" << body->id << ").\n";

    return;
}
/*----------------------------------------------------------------------------*/

/* Local granular fluidity model. */
void calculate_stress(Body *body, double dtIn)
{
    threadtask_t t;
    t.offset = 0;
    t.blocksize = body->p;
    calculate_stress_threaded(&t, body, dtIn);
    return;
}

/*----------------------------------------------------------------------------*/
void calculate_stress_threaded(threadtask_t *task, Body *body, double dtIn)
{
    dt = dtIn;

    /* Since this is local, we can split the particles among the threads. */
    size_t p_start = task->offset;
    size_t p_stop = task->offset + task->blocksize;

    /*
        cohesion, but other values should probably be changed to make this
        nonzero.
    *///////////////////////////////////////////////////////////////////////
    double c = 0;

    /* fprintf(stderr, "processing particle ids [%zu %zu].\n", p_start, p_stop); */

    for (size_t i = p_start; i < p_stop; i++) {
        if (job->active[i] == 0) {
            continue;
        }

        /* more convenient to use these forms for jaumann rate integration */
        const double exx_t = job->particles[i].L[XX];
        const double exy_t = 0.5 * (job->particles[i].L[XY] + job->particles[i].L[YX]);
        const double wxy_t = 0.5 * (job->particles[i].L[XY] - job->particles[i].L[YX]);
        const double eyy_t = job->particles[i].L[YY];

        /* trial elastic increment using jaumann rate */
        const double trD = exx_t + eyy_t;
        double dsjxx = lambda * trD + 2.0 * G * exx_t;
        double dsjxy = 2.0 * G * exy_t;
        double dsjyy = lambda * trD + 2.0 * G * eyy_t;
        double dsjzz = lambda * trD;
        dsjxx += 2 * wxy_t * job->particles[i].T[XY];
        dsjxy -= wxy_t * (job->particles[i].T[XX] - job->particles[i].T[YY]);
        dsjyy -= 2 * wxy_t * job->particles[i].T[XY];

        /* trial stress tensor */
        const double sxx_tr = job->particles[i].T[XX] + job->dt * dsjxx;
        const double sxy_tr = job->particles[i].T[XY] + job->dt * dsjxy;
        const double syy_tr = job->particles[i].T[YY] + job->dt * dsjyy;
        const double szz_tr = job->particles[i].T[ZZ] + job->dt * dsjzz;

        /* Calculate tau and p trial values. */

        /* Treat the intruder as a linear-elastic material. */
        if (job->particles[i].material == kIntruderMaterial) {
            const double scaling = 10;
            const double lambda_block = scaling * lambda;
            const double G_block = scaling * G;
            const double trD = exx_t + eyy_t;
            dsjxx = lambda_block * trD + 2.0 * G_block * exx_t;
            dsjxy = 2.0 * G_block * exy_t;
            dsjyy = lambda_block * trD + 2.0 * G_block * eyy_t;
            dsjxx += 2 * wxy_t * job->particles[i].T[XY];
            dsjxy -= wxy_t * (job->particles[i].T[XX] - job->particles[i].T[YY]);
            dsjyy -= 2 * wxy_t * job->particles[i].T[XY];
            double dsjzz = lambda_block * trD;

            const double sxx_tr = job->particles[i].T[XX] + job->dt * dsjxx;
            const double sxy_tr = job->particles[i].T[XY] + job->dt * dsjxy;
            const double syy_tr = job->particles[i].T[YY] + job->dt * dsjyy;
            const double szz_tr = job->particles[i].T[ZZ] + job->dt * dsjzz;
            job->particles[i].T[XX] = sxx_tr;
            job->particles[i].T[XY] = sxy_tr;
            job->particles[i].T[YY] = syy_tr;
            job->particles[i].T[ZZ] = szz_tr;

            continue;
        }

        // make everything linear elastic
        /*
            job->particles[i].T[XX] = sxx_tr;
            job->particles[i].T[XY] = sxy_tr;
            job->particles[i].T[YY] = syy_tr;
            job->particles[i].T[ZZ] = szz_tr;
            continue;
        */

        /* trial deviator values */
        const double p_tr = -(sxx_tr + syy_tr + szz_tr) / 3.0;
        const double t0xx_tr = sxx_tr + p_tr;
        const double t0xy_tr = sxy_tr;
        const double t0yy_tr = syy_tr + p_tr;
        const double t0zz_tr = szz_tr + p_tr;
        const double tau_tr = sqrt(0.5*(t0xx_tr*t0xx_tr + 2*t0xy_tr*t0xy_tr + t0yy_tr*t0yy_tr + t0zz_tr*t0zz_tr));

        bool density_flag = false;
        if ((job->particles[i].m / job->particles[i].v) < RHO_CRITICAL) {
            density_flag = true;
            /* printf("%4d: density %lf\n", i, (job->particles[i].m / job->particles[i].v));*/
        }

        double nup_tau;
        if (density_flag || p_tr <= c) {
            nup_tau = (tau_tr) / (G * job->dt);
            // beta = -p_tr / (K * job->dt);

            job->particles[i].T[XX] = 0;
            job->particles[i].T[XY] = 0;
            job->particles[i].T[YY] = 0;
            job->particles[i].T[ZZ] = 0;
        } else if (p_tr > c) {
            const double mu_scaling = 1;
            const double S0 = mu_scaling * MU_S * p_tr;
            double tau_tau;
            double scale_factor;
            if (tau_tr <= S0) {
                tau_tau = tau_tr;
                scale_factor = 1.0;
            } else {
                const double S2 = mu_scaling * MU_2 * p_tr;
                const double alpha = G * I_0 * job->dt * sqrt(p_tr / GRAINS_RHO) / grains_d;
                const double B = -(S2 + tau_tr + alpha);
                const double H = S2 * tau_tr + S0 * alpha;
                tau_tau = negative_root(1.0, B, H);
                scale_factor = (tau_tau / tau_tr);
            }

            nup_tau = ((tau_tr - tau_tau) / G) / job->dt;
            // beta = 0;

            job->particles[i].T[XX] = scale_factor * t0xx_tr - p_tr;
            job->particles[i].T[XY] = scale_factor * t0xy_tr;
            job->particles[i].T[YY] = scale_factor * t0yy_tr - p_tr;
            job->particles[i].T[ZZ] = scale_factor * t0zz_tr - p_tr;
        } else {
            /* fprintf(stderr, "u %zu %3.3g %3.3g %d ", i, f, p_tr, density_flag);*/
            fprintf(stderr, "u");
            nup_tau = 0;
        }

        /* use strain rate to calculate stress increment */
        job->particles[i].state[GAMMAP_STATE] += nup_tau * job->dt;
        job->particles[i].state[GAMMADOTP_STATE] = nup_tau;
    }

    return;
}
/*----------------------------------------------------------------------------*/


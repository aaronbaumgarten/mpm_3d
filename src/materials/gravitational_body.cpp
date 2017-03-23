//
// Created by aaron on 3/9/17.
// gravitational_body.cpp
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
#include <Eigen/Dense>
#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "process.hpp"
#include "tensor.hpp"

/* Estimated (see Kamrin and Koval 2014 paper) */
#define MU_S 0.6745//0.3819//0.280
#define GRAINS_RHO 2500//2450
#define RHO_CRITICAL (GRAINS_RHO*0.6)//1500
#define MU_2 MU_S//(I_0+MU_S)
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

extern "C" void volumetric_smoothing(Body *body, Eigen::VectorXd trE, Eigen::VectorXd trT);

extern "C" void volumetric_smoothing_implicit(Body *body, Eigen::VectorXd trE, Eigen::VectorXd trT);

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
static double rate_in;

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

void gravity_field(Body *body)
{
    //sum potentials over nodes with mass
    double Gr = 6.674e-11;
    double dx, dy, dz;
    Eigen::VectorXd g_potential(body->n);
    g_potential.setZero();
    /*for (size_t i = 0; i < body->n; i++){
        if (body->nodes.m[i] > 0){
            for (size_t p = 0; p < body->p; p++) {
                if (body->particles.active[p] == 1) {
                    dx = body->particles.x[p] - body->nodes.x[i];
                    dy = body->particles.y[p] - body->nodes.y[i];
                    dz = body->particles.z[p] - body->nodes.z[i];
                    g_potential[i] -= Gr * body->particles.m[p] / std::sqrt(dx * dx + dy * dy + dz * dz);
                }
            }
        }
    }*/
    body->particles.bx.setZero();// = -body->gradPhiX.transpose()*g_potential;
    body->particles.by.setZero();// = -body->gradPhiY.transpose()*g_potential;
    body->particles.bz.setZero();// = -body->gradPhiZ.transpose()*g_potential;

    //add coriolis and centrifugal terms
    double x_cm;
    double y_cm;
    double z_cm;
    x_cm = (body->particles.x.array()*body->particles.m.array()).sum()/body->particles.m.sum();
    y_cm = (body->particles.y.array()*body->particles.m.array()).sum()/body->particles.m.sum();
    z_cm = (body->particles.z.array()*body->particles.m.array()).sum()/body->particles.m.sum();

    Eigen::Vector3d v;
    Eigen::Vector3d r;
    Eigen::Vector3d w;
    Eigen::Vector3d b_co;
    Eigen::Vector3d b_ce;
    for (size_t p=0;p<body->p;p++){
        r[0] = body->particles.x[p] - x_cm;
        r[1] = body->particles.y[p] - y_cm;
        r[2] = body->particles.z[p] - z_cm;
        v[0] = body->particles.x_t[p];
        v[1] = body->particles.y_t[p];
        v[2] = body->particles.z_t[p];
        w.setZero();
        w[2] = rate_in;
        b_co = -2.0*w.cross(v);
        b_ce = -w.cross(w.cross(r));//-w.dot(r)*w + w.dot(w)*r;
        body->particles.bx[p] += b_co[0] + b_ce[0];
        body->particles.by[p] += b_co[1] + b_ce[1];
        body->particles.bz[p] += b_co[2] + b_ce[2];
    }
    return;
}

/*----------------------------------------------------------------------------*/
void material_init(Body *body)
{
    for (size_t i = 0; i < body->p; i++) {
        for (size_t j = 0; j < DEPVAR; j++) {
            body->particles.state(i,j) = 0;
        }
    }

    if (body->material.num_fp64_props < 4) {
        // Bit of a hack, but it's okay for now. Just close your eyes and code it anyways.
        fprintf(stderr,
                "%s:%s: Need at least 4 properties defined (E, nu, grain diameter, rate).\n",
                __FILE__, __func__);
        exit(EXIT_FAILURE);
    } else {
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        grains_d = body->material.fp64_props[2];
        rate_in = body->material.fp64_props[3];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2*nu));
        lambda = K - 2.0 * G / 3.0;

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, grain diameter = %g).\n",
               E, nu, G, K, grains_d);
    }

    for (size_t i = 0; i < body->p; i++) {
        body->particles.state(i, GAMMAP_STATE) = 0;
        body->particles.state(i, GAMMADOTP_STATE) = 0;
    }

    //z-axis rotation
    /*double x_cm;
    double y_cm;
    x_cm = (body->particles.x.array()*body->particles.m.array()).sum()/body->particles.m.sum();
    y_cm = (body->particles.y.array()*body->particles.m.array()).sum()/body->particles.m.sum();
    body->particles.x_t = -(body->particles.y.array() - y_cm)*rate_in;
    body->particles.y_t = (body->particles.x.array() - x_cm)*rate_in;*/

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
    double dt = dtIn;

    /* Since this is local, we can split the particles among the threads. */
    size_t p_start = task->offset;
    size_t p_stop = task->offset + task->blocksize;

    /*
        cohesion, but other values should probably be changed to make this
        nonzero.
    *///////////////////////////////////////////////////////////////////////
    double c = 0;

    /* fprintf(stderr, "processing particle ids [%zu %zu].\n", p_start, p_stop); */

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
        Eigen::Matrix3d W = 0.5*(L-L.transpose());

        double trD = D.trace();

        Eigen::Matrix3d gleft = W*T;
        Eigen::Matrix3d gright = T*W;

        Eigen::Matrix3d tmp = gleft-gright;

        Eigen::Matrix3d CD = 2*G*D + lambda*trD*Eigen::Matrix3d::Identity();

        Eigen::Matrix3d dsj = CD + tmp;

        //trial stress
        Eigen::Matrix3d s_tr = T + dt * dsj;
        double p_tr = -s_tr.trace() / 3.0;

        //trial deviator
        Eigen::Matrix3d t0_tr = s_tr + p_tr * Eigen::Matrix3d::Identity();
        double tau_tr = t0_tr.norm()/std::sqrt(2.0);

        bool density_flag = false;
        if ((body->particles.m[i] / body->particles.v[i]) < RHO_CRITICAL){
            density_flag = true;
        }

        double nup_tau;
        if (density_flag || p_tr <= c){
            nup_tau = tau_tr / (G * dt);
            for (size_t pos=0;pos<9;pos++){
                body->particles.T(i,pos) = 0;
            }
        } else if (p_tr > c) {
            double S0 = MU_S * p_tr;
            double tau_tau;
            double scale_factor;
            if (tau_tr <= S0) {
                tau_tau = tau_tr;
                scale_factor = 1.0;
            } else {
                double S2 = MU_2 * p_tr;
                double alpha = G * I_0 * dt * std::sqrt(p_tr / GRAINS_RHO) / grains_d;
                double B = -(S2 + tau_tr + alpha);
                double H = S2 * tau_tr + S0 * alpha;
                tau_tau = negative_root(1.0,B,H);
                scale_factor = (tau_tau/tau_tr);
            }

            nup_tau = ((tau_tr - tau_tau) / G) / dt;

            for (size_t pos=0;pos<9;pos++){
                body->particles.T(i,pos) = scale_factor * t0_tr(pos);
                if (pos == XX || pos == YY || pos == ZZ){
                    body->particles.T(i,pos) -= p_tr;
                }
            }
        } else {
            std::cerr << "u\n";
            nup_tau = 0;
        }

        body->particles.state(i,GAMMAP_STATE) += nup_tau * dt;
        body->particles.state(i,GAMMADOTP_STATE) = nup_tau;
    }

    gravity_field(body);

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_implicit(Body *body, double dtIn)
{
    double dt = dtIn;

    /*
        cohesion, but other values should probably be changed to make this
        nonzero.
    *///////////////////////////////////////////////////////////////////////
    double c = 0;

    /* fprintf(stderr, "processing particle ids [%zu %zu].\n", p_start, p_stop); */

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
        Eigen::Matrix3d W = 0.5*(L-L.transpose());

        double trD = D.trace();

        Eigen::Matrix3d gleft = W*T;
        Eigen::Matrix3d gright = T*W;

        Eigen::Matrix3d tmp = gleft-gright;

        Eigen::Matrix3d CD = 2*G*D + lambda*trD*Eigen::Matrix3d::Identity();

        Eigen::Matrix3d dsj = CD + tmp;

        //trial stress
        Eigen::Matrix3d s_tr = T + dt * dsj;
        double p_tr = -s_tr.trace() / 3.0;

        //trial deviator
        Eigen::Matrix3d t0_tr = s_tr + p_tr * Eigen::Matrix3d::Identity();
        double tau_tr = t0_tr.norm()/std::sqrt(2.0);

        bool density_flag = false;
        if ((body->particles.m[i] / body->particles.v[i]) < RHO_CRITICAL){
            density_flag = true;
        }

        double nup_tau;
        if (density_flag || p_tr <= c){
            nup_tau = tau_tr / (G * dt);
            for (size_t pos=0;pos<9;pos++){
                body->particles.Ttrial(i,pos) = 0;
            }
        } else if (p_tr > c) {
            double S0 = MU_S * p_tr;
            double tau_tau;
            double scale_factor;
            if (tau_tr <= S0) {
                tau_tau = tau_tr;
                scale_factor = 1.0;
            } else {
                double S2 = MU_2 * p_tr;
                double alpha = G * I_0 * dt * std::sqrt(p_tr / GRAINS_RHO) / grains_d;
                double B = -(S2 + tau_tr + alpha);
                double H = S2 * tau_tr + S0 * alpha;
                tau_tau = negative_root(1.0,B,H);
                scale_factor = (tau_tau/tau_tr);
            }

            nup_tau = ((tau_tr - tau_tau) / G) / dt;

            for (size_t pos=0;pos<9;pos++){
                body->particles.Ttrial(i,pos) = scale_factor * t0_tr(pos);
                if (pos == XX || pos == YY || pos == ZZ){
                    body->particles.Ttrial(i,pos) -= p_tr;
                }
            }
        } else {
            std::cerr << "u\n";
            nup_tau = 0;
        }
    }

    gravity_field(body);

    return;
}
/*----------------------------------------------------------------------------*/

void volumetric_smoothing(Body *body, Eigen::VectorXd trE, Eigen::VectorXd trT) {
    assert(trE.size() == body->p);
    assert(trT.size() == body->p);
    /*for (size_t i = 0; i < body->p; i++) {
        body->particles.T(i,XX) = trT[i]/3.0;
        body->particles.T(i,YY) = trT[i]/3.0;
        body->particles.T(i,ZZ) = trT[i]/3.0;
    }*/
    Eigen::VectorXd tmpVec(9);
    for (size_t i=0;i<body->p;i++) {
        tmpVec << body->particles.T.row(i).transpose();
        //Eigen::Matrix3d T(tmpVec.data());
        Eigen::Matrix<double, 3, 3, Eigen::RowMajor> T(tmpVec.data());
        T = T - (T.trace() - trT[i])/3.0*Eigen::Matrix3d::Identity();
        body->particles.T(i,XX) = T(XX);
        body->particles.T(i,YY) = T(YY);
        body->particles.T(i,ZZ) = T(ZZ);
    }

    //volume smoothing
    body->particles.v = body->particles.v0.array() * trE.array().exp();

    return;
}

void volumetric_smoothing_implicit(Body *body, Eigen::VectorXd trE, Eigen::VectorXd trT) {
    assert(trE.size() == body->p);
    assert(trT.size() == body->p);
    /*for (size_t i = 0; i < body->p; i++) {
        body->particles.Ttrial(i,XX) = trT[i]/3.0;
        body->particles.Ttrial(i,YY) = trT[i]/3.0;
        body->particles.Ttrial(i,ZZ) = trT[i]/3.0;
    }*/
    Eigen::VectorXd tmpVec(9);
    for (size_t i=0;i<body->p;i++) {
        tmpVec << body->particles.Ttrial.row(i).transpose();
        //Eigen::Matrix3d T(tmpVec.data());
        Eigen::Matrix<double, 3, 3, Eigen::RowMajor> T(tmpVec.data());
        T = T - (T.trace() - trT[i])/3.0*Eigen::Matrix3d::Identity();
        body->particles.Ttrial(i,XX) = T(XX);
        body->particles.Ttrial(i,YY) = T(YY);
        body->particles.Ttrial(i,ZZ) = T(ZZ);
    }

    //volume smoothing
    body->particles.v_trial = body->particles.v0.array() * trE.array().exp();

    return;
}
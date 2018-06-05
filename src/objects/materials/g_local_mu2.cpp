//
// Created by aaron on 5/25/18.
// g_local_mu2.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <dlfcn.h>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "materials.hpp"

/*----------------------------------------------------------------------------*/

void Sand_SachithLocal::quadratic_roots(double *x1, double *x2, double a, double b, double c)
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

double Sand_SachithLocal::negative_root(double a, double b, double c)
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

void Sand_SachithLocal::init(Job* job, Body* body){
    if (fp64_props.size() < 6) {
        // Bit of a hack, but it's okay for now. Just close your eyes and code it anyways.
        fprintf(stderr,
                "%s:%s: Need at least 6 properties defined (E, nu, grain diameter, mu_s, grains_rho, rho_critical).\n",
                __FILE__, __func__);
        exit(EXIT_FAILURE);
    } else if (fp64_props.size() < 8){
        E = fp64_props[0];
        nu = fp64_props[1];
        grains_d = fp64_props[2];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2*nu));
        lambda = K - 2.0 * G / 3.0;

        MU_S = fp64_props[3];
        GRAINS_RHO = fp64_props[4];
        RHO_CRITICAL = fp64_props[5];

        MU_2 = MU_S;
        DELTA_MU = MU_2 - MU_S;
        I_0 = 1;

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, grain diameter = %g, mu_s = %g, grains_rho = %g, rho_critical = %g).\n",
               E, nu, G, K, grains_d, MU_S, GRAINS_RHO, RHO_CRITICAL);
    } else {
        E = fp64_props[0];
        nu = fp64_props[1];
        grains_d = fp64_props[2];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2*nu));
        lambda = K - 2.0 * G / 3.0;

        MU_S = fp64_props[3];
        GRAINS_RHO = fp64_props[4];
        RHO_CRITICAL = fp64_props[5];

        MU_2 = fp64_props[6];
        DELTA_MU = MU_2 - MU_S;
        I_0 = fp64_props[7];

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, grain diameter = %g, mu_s = %g, grains_rho = %g, rho_critical = %g, mu_2 = %g, i_0 = %g).\n",
               E, nu, G, K, grains_d, MU_S, GRAINS_RHO, RHO_CRITICAL, MU_2, I_0);
    }

    gammap_dot.resize(body->points->x.size());
    gammap_dot.setZero();
    gammap.resize(body->points->x.size());
    gammap.setZero();

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void Sand_SachithLocal::writeFrame(Job* job, Body* body, Serializer* serializer) {
    serializer->writeScalarArray(gammap,"gammap");
    serializer->writeScalarArray(gammap_dot,"gammap_dot");
    return;
}

/*----------------------------------------------------------------------------*/

std::string Sand_SachithLocal::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

/*----------------------------------------------------------------------------*/

int Sand_SachithLocal::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}

/*----------------------------------------------------------------------------*/

void Sand_SachithLocal::calculateStress(Job* job, Body* body, int SPEC){
    /*
        cohesion, but other values should probably be changed to make this
        nonzero.
    *///////////////////////////////////////////////////////////////////////
    double c = 0;

    MaterialTensor T, L, D, W;
    MaterialTensor s_tr, t0_tr;

    double trD, p_tr, tau_tr, nup_tau;
    double S0, tau_tau, scale_factor;
    double S2, alpha, B, H;

    bool density_flag;

    MaterialTensor tmpMat;

    for (int i=0;i<body->points->x.size();i++){
        if (body->points->active[i] == 0){
            continue;
        }

        L = body->points->L(i);
        T = body->points->T(i);

        D = 0.5*(L+L.transpose());
        W = 0.5*(L-L.transpose());

        trD = D.trace();

        tmpMat = (2*G*D) + (lambda*trD*MaterialTensor::Identity()) + (W*T) - (T*W);

        //trial stress
        s_tr = T + job->dt * tmpMat;
        p_tr = -s_tr.trace()/3.0;

        //trial deviator
        t0_tr = s_tr + p_tr * MaterialTensor::Identity();
        tau_tr = t0_tr.norm()/std::sqrt(2.0);

        density_flag = false;
        if ((body->points->m(i) / body->points->v(i)) < RHO_CRITICAL){
            density_flag = true;
        }

        if (density_flag || p_tr <= c){
            nup_tau = tau_tr / (G * job->dt);
            body->points->T(i).setZero();
        } else if (p_tr > c) {
            S0 = MU_S * p_tr;
            if (tau_tr <= S0) {
                tau_tau = tau_tr;
                scale_factor = 1.0;
            } else {
                S2 = MU_2 * p_tr;
                alpha = G * I_0 * job->dt * std::sqrt(p_tr / GRAINS_RHO) / grains_d;
                B = -(S2 + tau_tr + alpha);
                H = S2 * tau_tr + S0 * alpha;
                tau_tau = negative_root(1.0,B,H);
                scale_factor = (tau_tau/tau_tr);
            }

            nup_tau = ((tau_tr - tau_tau) / G) / job->dt;

            tmpMat = scale_factor * t0_tr - p_tr * MaterialTensor::Identity();
            body->points->T(i) = tmpMat;

        } else {
            if (SPEC == Material::UPDATE) {
                std::cerr << "u\n";
            }
            nup_tau = 0;
        }

        if (SPEC == Material::UPDATE) {
            gammap(i) += nup_tau * job->dt;
            gammap_dot(i) = nup_tau;
        }
    }

    return;
}

/*----------------------------------------------------------------------------*/

void Sand_SachithLocal::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    for (int pos=0;pos<stressIN.size();pos++){
        body->points->T(idIN,pos) = stressIN(pos);
    }
    return;
}

/*----------------------------------------------------------------------------*/

void Sand_SachithLocal::assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    MaterialTensor tmpMat;
    Eigen::VectorXd tmpVec;
    double trT;

    tmpMat = body->points->T(idIN);

    trT = tmpMat.trace();
    tmpMat = tmpMat - (trT/3.0 + pressureIN) * MaterialTensor::Identity();
    body->points->T(idIN) = tmpMat;
    return;
}

/*----------------------------------------------------------------------------*/
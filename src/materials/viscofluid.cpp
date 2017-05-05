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
#include <unsupported/Eigen/MatrixFunctions>

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

extern "C" void volumetric_smoothing(Body *body, Eigen::VectorXd trE, Eigen::VectorXd trT);

extern "C" void volumetric_smoothing_implicit(Body *body, Eigen::VectorXd trE, Eigen::VectorXd trT);

#define J_STATE 8
#define J_TRIAL_STATE 7

/*----------------------------------------------------------------------------*/
void material_init(Body *body) {
    size_t i, j;

    for (i = 0; i < body->p; i++) {
        for (j = 0; j < DEPVAR; j++) {
            body->particles.state(i,j) = 0;
        }
        body->particles.state(i,J_STATE) = 1;
        body->particles.state(i,J_TRIAL_STATE) = 1;
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
        Eigen::Matrix<double, 3, 3, Eigen::RowMajor> L(tmpVec.data());

        Eigen::Matrix3d D = 0.5*(L+L.transpose());

        //double trT = T.trace();
        double trD = D.trace();

        //if in slurry mixture, state9 will be 1
        if (body->particles.state(i,9) != 0){
            body->particles.state(i,J_TRIAL_STATE) *= std::exp(dt*body->particles.state(i,10));
        } else {
            body->particles.state(i,J_TRIAL_STATE) *= std::exp(dt*L.trace());
        }

        double J = body->particles.state(i,J_TRIAL_STATE);
        //J = body->particles.v[i] / body->particles.v0[i];

        Eigen::Matrix3d Sv = 2 * mu * (D - trD/3.0 * Eigen::Matrix3d::Identity());
        Eigen::Matrix3d Se = K*std::log(J)*Eigen::Matrix3d::Identity();

        for (size_t pos=0;pos<9;pos++){
            body->particles.Ttrial(i,pos) = Se(pos) + Sv(pos);
        }
    }

    return;

}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_threaded(threadtask_t *task, Body *body, double dtIn) {
    dt = dtIn;
    size_t p_start = task->offset;
    size_t p_stop = task->offset + task->blocksize;

    body->particles.state.col(J_TRIAL_STATE) = body->particles.state.col(J_STATE);
    body->material.calculate_stress_implicit(body,dtIn);
    body->particles.T = body->particles.Ttrial;
    //adjust state
    body->particles.state.col(J_STATE) = body->particles.state.col(J_TRIAL_STATE);

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
        body->particles.state(i,J_STATE) = std::exp(trT[i]/(3*K));
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
        body->particles.state(i,J_TRIAL_STATE) = std::exp(trT[i]/(3*K));
    }

    //volume smoothing
    body->particles.v_trial = body->particles.v0.array() * trE.array().exp();

    return;
}


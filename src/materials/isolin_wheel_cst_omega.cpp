//
// Created by aaron on 2/28/17.
// isolin_wheel_cst_omega.cpp
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
double t = 0;
double t_start = 0;
double w_min = 1e0;

/* Material Constants (set by configuration file). */
double E, nu, omega_in, axle_x, axle_y, axle_z;
double G, K, lambda;

extern "C" void material_init(Body *body);

extern "C" void calculate_stress_threaded(threadtask_t *task, Body *body, double dt);

extern "C" void calculate_stress(Body *body, double dt);

extern "C" void calculate_stress_implicit(Body *body, double dt);

extern "C" void volumetric_smoothing(Body* body, Eigen::VectorXd trE, Eigen::VectorXd trT);

extern "C" void volumetric_smoothing_implicit(Body * body, Eigen::VectorXd trE, Eigen::VectorXd trT);

#define SZZ_STATE 0

/*----------------------------------------------------------------------------*/
void material_init(Body *body) {
    size_t i, j;

    for (i = 0; i < body->p; i++) {
        for (j = 0; j < DEPVAR; j++) {
            body->particles.state(i,j) = 0;
        }
    }

    if (body->material.num_fp64_props < 5) {
        std::cout << body->material.num_fp64_props << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 5 properties defined (E, nu, omega, [t_start], [axle_x, axle_y, axle_z]).\n",
                __FILE__, __func__);
        exit(0); //replace error code handler from sachith's work with 0
    } else if (body->material.num_fp64_props == 7) {
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        omega_in = body->material.fp64_props[2];
        t_start = body->material.fp64_props[3];
        axle_x = body->material.fp64_props[4];
        axle_y = body->material.fp64_props[5];
        axle_z = body->material.fp64_props[6];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, omega = %g, Axle = < %g , %g , %g >).\n",
               E, nu, G, K, omega_in, axle_x, axle_y, axle_z);
    } else if (body->material.num_fp64_props == 4) {
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        omega_in = body->material.fp64_props[2];
        t_start = body->material.fp64_props[3];
        axle_x = 0;
        axle_y = 0;
        axle_z = 1;
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, omega = %g, Axle = < %g , %g , %g >).\n",
               E, nu, G, K, omega_in, axle_x, axle_y, axle_z);
    } else {
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        omega_in = body->material.fp64_props[2];
        t_start = 0;
        axle_x = 0;
        axle_y = 0;
        axle_z = 1;
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, omega = %g, Axle = < %g , %g , %g >).\n",
               E, nu, G, K, omega_in, axle_x, axle_y, axle_z);
    }

    std::cout << "Done initializing material (" << body->id << ").\n";
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress(Body *body, double dtIn) {
    dt = dtIn;
    t += dtIn;
    threadtask_t t_kk;
    t_kk.offset = 0;
    t_kk.blocksize = body->p;
    calculate_stress_threaded(&t_kk,body, dtIn);
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_implicit(Body *body, double dtIn) {

    dt = dtIn;

    //find properties of center of wheel
    Eigen::Vector3d r_cp;
    Eigen::Vector3d v_cp;
    Eigen::Vector3d mr_cp;
    Eigen::Vector3d mv_cp;
    double m_cp = 0;

    mr_cp.setZero();
    mv_cp.setZero();

    for (size_t i = 0; i < body->p; i++) {
        mr_cp[0] += body->particles.m[i]*body->particles.x[i];
        mr_cp[1] += body->particles.m[i]*body->particles.y[i];
        mr_cp[2] += body->particles.m[i]*body->particles.z[i];

        mv_cp[0] += body->particles.m[i]*body->particles.x_t[i];
        mv_cp[1] += body->particles.m[i]*body->particles.y_t[i];
        mv_cp[2] += body->particles.m[i]*body->particles.z_t[i];

        m_cp += body->particles.m[i];
    }

    r_cp = mr_cp/m_cp;
    v_cp = mv_cp/m_cp;

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

        Eigen::Matrix3d tmpMat = gleft-gright;

        Eigen::Matrix3d CD = 2*G*D + lambda*trD*Eigen::Matrix3d::Identity();

        Eigen::Matrix3d dsj = CD + tmpMat;

        for (size_t pos=0;pos<9;pos++){
            body->particles.Ttrial(i,pos) = body->particles.T(i,pos) + dt*dsj(pos);
        }

        /* add in body force for rolling */
        if (t >= t_start) {
            Eigen::Vector3d r_p;
            r_p << body->particles.x[i], body->particles.y[i], body->particles.z[i];

            Eigen::Vector3d b = r_p - r_cp;
            Eigen::Vector3d a;
            a << axle_x, axle_y, axle_z;

            //radius from axle
            Eigen::Vector3d r = b - a * (a.dot(b)) / a.squaredNorm();

            //direction of rotation
            Eigen::Vector3d c = a.cross(r)/(a.norm() * r.norm());

            //relative velocity
            Eigen::Vector3d dv_cp;
            dv_cp << (body->particles.x_t[i] - v_cp[0]), (body->particles.y_t[i] - v_cp[1]), (body->particles.z_t[i] - v_cp[2]);

            //normal velocity
            Eigen::Vector3d v_n = dv_cp.dot(r)*r/r.squaredNorm();

            //theta velocity
            Eigen::Vector3d v_t = omega_in*r.norm()*c;

            body->particles.x_t[i] = v_n[0] + v_t[0] + v_cp[0];
            body->particles.y_t[i] = v_n[1] + v_t[1] + v_cp[1];
            body->particles.z_t[i] = v_n[2] + v_t[2] + v_cp[2];

        }

    }
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_threaded(threadtask_t *task, Body *body, double dtIn) {
    //printf("t = %g\nx = %g\ny = %g\n", body->t, body->particles[0].x, body->particles[0].y);
    /*printf("L: \n%g %g %g\n%g %g %g\n%g %g %g\n\n",
           body->particles[0].L[0],
           body->particles[0].L[1],
           body->particles[0].L[2],
           body->particles[0].L[3],
           body->particles[0].L[4],
           body->particles[0].L[5],
           body->particles[0].L[6],
           body->particles[0].L[7],
           body->particles[0].L[8]);*/
    /* Since this is local, we can split the particles among the threads. */
    dt = dtIn;
    calculate_stress_implicit(body,dt);

    body->particles.T = body->particles.Ttrial;

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
    /*Eigen::VectorXd tmpVec(9);
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
    body->particles.v = body->particles.v0.array() * trE.array().exp();*/

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
    /*Eigen::VectorXd tmpVec(9);
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
    body->particles.v_trial = body->particles.v0.array() * trE.array().exp();*/
    return;
}


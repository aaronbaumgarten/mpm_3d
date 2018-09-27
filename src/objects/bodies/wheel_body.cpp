//
// Created by aaron on 6/18/18.
// wheel_body.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "bodies.hpp"

/*----------------------------------------------------------------------------*/
//
void WheelBody::init(Job* job){
    if (fp64_props.size() < 2){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined ({omega, t_start}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        omega = fp64_props[0];
        t_start = fp64_props[1];

        if (job->DIM == 1){
            std::cerr << "CANNOT RUN WHEELS IN 1D! Exiting." << std::endl;
            exit(0);
        } else if (job->DIM == 2){
            a.setZero();
            a[2] = 1; //positive z
        } else if (job->DIM == 3 && fp64_props.size() < 5){
            std::cout << fp64_props.size() << "\n";
            fprintf(stderr,
                    "%s:%s: Need at least 5 properties defined ({omega, t_start, <axle_vector>}).\n",
                    __FILE__, __func__);
            exit(0);
        } else {
            a[0] = fp64_props[2];
            a[1] = fp64_props[3];
            a[2] = fp64_props[4];
        }
    }

    r = KinematicVector(job->JOB_TYPE);
    mv_cp = KinematicVector(job->JOB_TYPE);
    v_cp = KinematicVector(job->JOB_TYPE);
    mx_cp = KinematicVector(job->JOB_TYPE);
    x_cp = KinematicVector(job->JOB_TYPE);
    dv_cp = KinematicVector(job->JOB_TYPE);
    v_n = KinematicVector(job->JOB_TYPE);
    v_t = KinematicVector(job->JOB_TYPE);

    points->init(job,this);
    nodes->init(job,this);

    if (activeMaterial != 0) {
        material->init(job, this);
    }
    if (activeBoundary != 0) {
        boundary->init(job, this);
    }

    adjusted_x = KinematicVectorArray(points->x.size(), job->JOB_TYPE);

    //find bounds of box
    Lx = KinematicVector(job->JOB_TYPE);
    Lx.setZero();
    for (int i=0; i < nodes->x.size(); i++){
        for (int pos=0; pos < nodes->x.DIM; pos++){
            if (nodes->x(i,pos) > Lx(pos)){
                Lx(pos) = nodes->x(i,pos);
            }
        }
    }

    //determine initial center of mass
    m_cp = 0;
    mx_cp.setZero();
    for (int i=0; i<points->x.size(); i++){
        if (points->active[i] != 0) {
            m_cp += points->m[i];
            mx_cp += points->m[i] * points->x[i];
        }
    }
    x_cp = mx_cp/m_cp;

    //assign vector type
    S = MPMScalarSparseMatrix(nodes->x.size(), points->x.size());
    gradS = KinematicVectorSparseMatrix(nodes->x.size(), points->x.size(), job->JOB_TYPE);


    printf("Body properties (omega = %g, t_start = %g, <axle_vector> = <%g,%g,%g>)\n" ,omega, t_start, a[0], a[1], a[2]);
    std::cout << "Body Initialized: [" << name << "]." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
/*
void WheelBody::generateLoads(Job* job){
    points->generateLoads(job, this);
    nodes->generateLoads(job, this);

    if (job->t >= t_start){
        //determine properties of wheel
        m_cp = 0;
        mx_cp.setZero();
        mv_cp.setZero();
        for (int i=0; i<nodes->x.size(); i++){
            if (nodes->m[i] > 0) {
                m_cp += nodes->m[i];
                mx_cp += nodes->m[i] * nodes->x[i];
                //mv_cp += nodes->mx_t[i] + nodes->f[i] * job->dt;
                mv_cp += nodes->mx_t[i];
            }
        }
        x_cp = mx_cp/m_cp;
        v_cp = mv_cp/m_cp;
    }

    return;
}

void WheelBody::applyLoads(Job* job){
    points->applyLoads(job, this);
    nodes->applyLoads(job, this);

    double tmp;
    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    tmpVec.setZero();
    if (job->t >= t_start){
        for (int i=0; i < nodes->x.size(); i++){
            if (nodes->m[i] > 0) {
                //relative position
                b = nodes->x[i] - x_cp;

                //radius from axle
                tmp = a.norm();
                r = KinematicVector(b - a * (a.dot(b)) / (tmp * tmp), job->JOB_TYPE);

                //direction of rotation
                c = a.cross(r) / (tmp * r.norm());

                //relative velocity
                //dv_cp = nodes->x_t[i] + nodes->f[i]/(nodes->m[i]) * job->dt - v_cp;
                dv_cp = nodes->x_t[i] - v_cp;

                //normal velocity
                tmp = r.norm();
                v_n = dv_cp.dot(r) * r / (tmp * tmp);

                //theta velocity
                v_t = omega * tmp * KinematicVector(c, job->JOB_TYPE);

                //update forces on nodes
                nodes->f[i] += 1.0/job->dt * (nodes->m[i]*(v_cp + v_n + v_t) - nodes->mx_t[i]);
                tmpVec += 1.0/job->dt * (nodes->m[i]*(v_cp + v_n + v_t) - nodes->mx_t[i]);
            }
        }

        std::cout << tmpVec.norm() << " ?= 0" << std::endl;

    }

    return;
}
*/

void WheelBody::generateLoads(Job* job){
    points->generateLoads(job, this);
    nodes->generateLoads(job, this);

    if (job->t >= t_start){
        //fix old x_cp
        for (int pos=0; pos<job->DIM; pos++){
            if (x_cp(pos) > Lx(pos)){
                //center wrapped around axis
                x_cp(pos) -= Lx(pos);
            } else if (x_cp(pos) < 0){
                //center wrapped around axis
                x_cp(pos) += Lx(pos);
            }
        }

        //find 'adjusted' point positions
        for (int i=0; i<points->x.size(); i++){
            for (int pos=0; pos<job->DIM; pos++){
                if (points->x(i,pos) - x_cp(pos) > Lx(pos)/2.0){
                    //point wrapped around axis
                    adjusted_x(i,pos) = points->x(i,pos) - Lx(pos);
                } else if (x_cp(pos) - points->x(i,pos) > Lx(pos)/2.0){
                    //point wrapped around axis
                    adjusted_x(i,pos) = points->x(i,pos) + Lx(pos);
                } else {
                    adjusted_x(i,pos) = points->x(i,pos);
                }
            }
        }

        //determine properties of wheel
        m_cp = 0;
        mx_cp.setZero();
        mv_cp.setZero();
        for (int i=0; i<points->x.size(); i++){
            if (points->active[i] != 0) {
                m_cp += points->m[i];
                //mx_cp += points->m[i] * points->x[i];
                mx_cp += points->m[i] * adjusted_x[i];
                mv_cp += points->mx_t[i];
            }
        }
        x_cp = mx_cp/m_cp;
        v_cp = mv_cp/m_cp;

        I_cp = 0;
        L_cp.setZero();
        for (int i=0; i<points->x.size(); i++){
            if (points->active[i] != 0) {
                //calculate moment of interia
                //I_cp += points->m[i] * (points->x[i] - x_cp).dot(points->x[i] - x_cp);
                I_cp += points->m[i] * (adjusted_x[i] - x_cp).dot(adjusted_x[i] - x_cp);
                //calculate angular momentum
                //L_cp += points->m[i] * (points->x[i] - x_cp).cross(points->x_t[i]);
                L_cp += points->m[i] * (adjusted_x[i] - x_cp).cross(points->x_t[i]);
            }
        }

        //calculate required delta omega (vector)
        d_omega = omega*a/a.norm() - L_cp/I_cp;
    }

    return;
}

void WheelBody::applyLoads(Job* job){
    points->applyLoads(job, this);
    nodes->applyLoads(job, this);

    double tmp;
    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    tmpVec.setZero();
    if (job->t >= t_start){
        for (int i=0; i < points->x.size(); i++){
            if (points->active[i] != 0) {
                //relative position
                //b = points->x[i] - x_cp;
                b = adjusted_x[i] - x_cp;

                //radius from d_omega
                tmp = d_omega.norm();
                r = KinematicVector(b - (d_omega * (d_omega.dot(b)) / (tmp * tmp)), job->JOB_TYPE);

                //direction of rotation
                c = d_omega.cross(r) / (tmp * r.norm());

                //relative velocity
                //dv_cp = nodes->x_t[i] + nodes->f[i]/(nodes->m[i]) * job->dt - v_cp;
                dv_cp = points->x_t[i] - v_cp;

                //normal velocity
                //tmp = r.norm();
                //v_n = dv_cp.dot(r) * r / (tmp * tmp);
                //v_n.setZero();

                //theta velocity
                v_t = d_omega.norm() * r.norm() * KinematicVector(c, job->JOB_TYPE);

                //update forces on nodes
                points->x_t[i] += v_t;
                points->mx_t[i] = points->m[i] * points->x_t[i];
                //tmpVec += points->m[i] * (v_t);

                //tmpVec += points->m[i]*(v_n + v_t);
            }
        }

        //std::cout << tmpVec.norm()/m_cp << " ?= 0" << std::endl;

    }

    return;
}
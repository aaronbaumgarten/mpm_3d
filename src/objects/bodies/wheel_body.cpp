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

    points->init(job,this);
    nodes->init(job,this);
    material->init(job,this);
    boundary->init(job,this);


    //assign vector type
    S = MPMScalarSparseMatrix(nodes->x.size(), points->x.size());
    gradS = KinematicVectorSparseMatrix(nodes->x.size(), points->x.size(), job->JOB_TYPE);


    printf("Body properties (omega = %g, t_start = %g, <axle_vector> = <%g,%g,%g>" ,omega, t_start, a[0], a[1], a[2]);
    std::cout << "Body Initialized: [" << name << "]." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
void WheelBody::generateLoads(Job* job){
    points->generateLoads(job, this);
    nodes->generateLoads(job, this);

    double tmp;
    if (job->t >= t_start){
        //determine properties of wheel
        m_cp = 0;
        mx_cp.setZero();
        mv_cp.setZero();
        for (int i=0; i<nodes->x.size(); i++){
            if (nodes->m[i] > 0) {
                m_cp += nodes->m[i];
                mx_cp += nodes->m[i] * nodes->x[i];
                mv_cp += nodes->mx_t[i] + nodes->f[i] * job->dt;
            }
        }
        x_cp = mx_cp/m_cp;
        v_cp = mv_cp/m_cp;

        for (int i=0; i < nodes->x.size(); i++){
            if (nodes->m[i] > 0) {
                //relative position
                b = nodes->x[i] - x_cp + (nodes->f[i] * job->dt)/nodes->m[i];

                //radius from axle
                tmp = a.norm();
                r = KinematicVector(b - a * (a.dot(b)) / (tmp * tmp), job->JOB_TYPE);

                //direction of rotation
                c = a.cross(r) / (tmp * r.norm());

                //relative velocity
                dv_cp = nodes->x_t[i] - v_cp;

                //normal velocity
                tmp = r.norm();
                v_n = dv_cp.dot(r) / (tmp * tmp);

                //theta velocity
                v_t = omega * tmp * KinematicVector(c, job->JOB_TYPE);

                //update forces on nodes
                nodes->f[i] += 1.0/job->dt * (nodes->m[i]*(v_cp + v_n + v_t) - nodes->mx_t[i]);
            }
        }
    }

    return;
}

void WheelBody::applyLoads(Job* job){
    points->applyLoads(job, this);
    nodes->applyLoads(job, this);
    return;
}
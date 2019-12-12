//
// Created by aaron on 5/25/18.
// viscousfluid.cpp
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

void BarotropicViscousFluid::init(Job* job, Body* body){
    if (fp64_props.size() < 2){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined (K, eta).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        K = fp64_props[0];
        eta = fp64_props[1];
        printf("Material properties (K = %g, mu = %g)",
               K, eta);
    }

    if (fp64_props.size() == 4){
        h = fp64_props[2];
        alpha = fp64_props[3];
        use_delta_correction = true;
        std::cout << "[" << h << ", " << alpha << "]." << std::endl;
    } else {
        h = 0;
        alpha = 0;
        use_delta_correction = false;
        std::cout << std::endl;
    }

    V_i.resize(job->grid->node_count,1);
    v_i.resize(job->grid->node_count,1);
    e.resize(job->grid->node_count,1);
    grad_e = KinematicVectorArray(body->points->x.size(),job->JOB_TYPE);
    del_pos = KinematicVectorArray(body->points->x.size(),job->JOB_TYPE);
    del_pos.setZero();

    for (int i=0; i<job->grid->node_count;i++){
        V_i(i) = job->grid->nodeVolume(job,i);
    }

    Lx = KinematicVector(job->JOB_TYPE);
    Lx.setZero();

    for (int pos=0;pos<Lx.rows();pos++){
        for (int i=0;i<body->nodes->x.size();i++){
            if (body->nodes->x(i,pos) > Lx(pos)){
                Lx(pos) = body->nodes->x(i,pos);
            }
        }
    }

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void BarotropicViscousFluid::calculateStress(Job* job, Body* body, int SPEC){
    MaterialTensor T, L, D, tmpMat;

    double trD;
    double alpha = 1.0;
    double w = 0;
    int max_index = -1;
    int rnd_num;
    double tmpNum;
    KinematicVector delta = KinematicVector(job->JOB_TYPE);

    if (use_delta_correction) {
        //calculate nodal volume integral
        if (job->JOB_TYPE == job->JOB_AXISYM) {
            //need to adjust point integration area
            Eigen::VectorXd A_tmp = body->points->v;
            for (int i = 0; i < A_tmp.rows(); i++) {
                A_tmp(i) /= body->points->x(i, 0);
            }
            v_i = body->S * A_tmp;
        } else {
            //otherwise integration area and volume are the same
            v_i = body->S * body->points->v;
        }

        for (int i = 0; i < V_i.rows(); i++) {
            tmpNum = (v_i(i) - V_i(i)) / (V_i(i));
            e(i) = std::max(0.0, tmpNum);
        }

        //calculate gradient of nodal overshoot
        grad_e = body->gradS.operate(e, MPMSparseMatrixBase::TRANSPOSED);
    }

    //strain smoothing
    Eigen::VectorXd J = Eigen::VectorXd(body->points->x.size());
    for (int i=0;i<J.size();i++){
        J(i) = body->points->v(i) / body->points->v0(i);
    }

    for (int i=0;i<body->points->x.size();i++){
        if (body->points->active[i] == 0){
            continue;
        }

        L = body->points->L(i);

        D = 0.5*(L+L.transpose());

        trD = D.trace();

        //stress update after density correction
        if (use_delta_correction) {

            delta = -alpha * job->dt * (L - (trD / 3.0) * MaterialTensor::Identity()).norm() * grad_e(i) * h * h;

            for (int pos = 0; pos < delta.rows(); pos++) {
                if (//((body->points->x(i, pos) + delta(pos)) >= Lx(pos)) ||
                    //((body->points->x(i, pos) + delta(pos)) <= 0) ||
                    !std::isfinite(delta(pos))) {
                    delta(pos) = 0;
                }
            }

            body->points->x(i) += delta;
            body->points->u(i) += delta;
            del_pos(i) += delta;

            body->points->x_t(i) += body->points->L(i) * delta;
            body->points->mx_t(i) = body->points->m(i) * body->points->x_t(i);
        }

        //T = 2*mu*D_0 + K*log(J)*I
        T = 2*eta*(D - (trD/3.0)*MaterialTensor::Identity()) + K*std::log(J(i))*MaterialTensor::Identity();

        body->points->T(i) = T;
    }

    return;
}

void BarotropicViscousFluid::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    body->points->T(idIN) = stressIN;
    double pressureIN = -stressIN.trace() / 3.0;

    body->points->v(idIN) = body->points->v0(idIN)*std::exp(-pressureIN/K); //p = -K*log(J)

    return;
}

void BarotropicViscousFluid::assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    MaterialTensor tmpMat;
    double trT;

    tmpMat = body->points->T(idIN);

    trT = tmpMat.trace();
    tmpMat = tmpMat - (trT/3.0 + pressureIN) * MaterialTensor::Identity();
    body->points->T(idIN) = tmpMat;

    body->points->v(idIN) = body->points->v0(idIN)*std::exp(-pressureIN/K); //p = -K*log(J)

    return;
}


/*----------------------------------------------------------------------------*/

void BarotropicViscousFluid::writeFrame(Job* job, Body* body, Serializer* serializer){
    if (use_delta_correction) {
        serializer->writeVectorArray(grad_e, "grad_err");
        serializer->writeVectorArray(del_pos, "del_pos");
        serializer->writeScalarArray(V_i, "grid_volume");
        serializer->writeScalarArray(e, "err");
    }

    Eigen::VectorXd nvec(body->nodes->x.size());
    Eigen::VectorXd pvec(body->points->x.size());
    Eigen::VectorXd tmpVec;
    MaterialTensor T;

    for (int i=0;i<body->points->x.size();i++) {
        if (body->points->active[i] == 0) {
            pvec(i) = 0;
            continue;
        }

        T = body->points->T(i);
        pvec(i) = -T.trace() * body->points->v(i) / 3.0;
    }

    nvec = body->S * pvec;
    v_i = body->S * body->points->v;

    for (int i=0;i<nvec.size();i++){
        nvec(i) = nvec(i) / v_i(i);
    }

    pvec = body->S.operate(nvec, MPMSparseMatrixBase::TRANSPOSED);

    serializer->writeScalarArray(pvec, "pressure_smooth");
    //nothing to report
    return;
}

std::string BarotropicViscousFluid::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int BarotropicViscousFluid::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}
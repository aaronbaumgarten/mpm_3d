//
// Created by aaron on 5/25/18.
// slurry_viscousfluid.cpp
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

void SlurryFluidPhase::init(Job* job, Body* body){
    if (fp64_props.size() < 3 || (str_props.size() < 1 && int_props.size() < 1)){
        std::cout << fp64_props.size() <<", " << str_props.size() << ", " << int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 4 properties defined ({K, eta, solid_rho}, {solid_body}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        K = fp64_props[0];
        eta_0 = fp64_props[1];
        solid_rho = fp64_props[2];

        //set body ids by name
        if (str_props.size() == 1){
            for (int b = 0; b < job->bodies.size(); b++) {
                if (str_props[0].compare(job->bodies[b]->name) == 0){
                    solid_body_id = b;
                    break;
                }
            }
        }

        // or set body ids by int
        if (solid_body_id < 0){
            if (int_props.size() == 1) {
                solid_body_id = int_props[0];
            } else {
                std::cout << fp64_props.size() <<", " << str_props.size() << ", " << int_props.size() << "\n";
                fprintf(stderr,
                        "%s:%s: Need at least 3 properties defined ({K, eta, solid_rho}, {solid_body}).\n",
                        __FILE__, __func__);
                exit(0);
            }
        }

        printf("Material properties ({K = %g, eta = %g, solid_rho = %g}, {solid_body: %i}).",
               K, eta_0, solid_rho, solid_body_id);
    }

    if (fp64_props.size() > 3){
        h = fp64_props[3];
        alpha = fp64_props[4];std::cout << "[" << h << ", " << alpha << "]." << std::endl;
    } else {
        h = 0;
        alpha = 0;
        std::cout << std::endl;
    }

    J.resize(body->points->x.size());
    J.setOnes();
    n_p.resize(body->points->x.size());
    n_p.setOnes();

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

void SlurryFluidPhase::calculateStress(Job* job, Body* body, int SPEC){
    MaterialTensor T, L, D;

    double m1, m2;
    KinematicVector mv1i, mv2i;

    Eigen::VectorXd n(body->nodes->x.size());
    Eigen::VectorXd pvec(body->points->x.size());
    Eigen::VectorXd m1_axi(job->bodies[solid_body_id]->nodes->x.size());
    //Eigen::VectorXd n_p(body->points.x.rows());
    KinematicVectorArray nMat = KinematicVectorArray(body->nodes->x.size(), job->JOB_TYPE);

    double trD, tmpNum;

    //perturbation variables
    KinematicVector delta = KinematicVector(job->JOB_TYPE);

    //calculate nodal volume integral
    if (job->JOB_TYPE == job->JOB_AXISYM){
        //need to adjust point integration area
        Eigen::VectorXd A_tmp = body->points->v;
        Eigen::VectorXd B_tmp = job->bodies[solid_body_id]->points->m;
        for (int i=0; i<A_tmp.rows(); i++){
            A_tmp(i) /= body->points->x(i,0);
        }
        v_i = body->S * A_tmp;

        //need to adjust 2D integral of mass
        for (int i=0; i<B_tmp.rows(); i++){
            B_tmp(i) /= job->bodies[solid_body_id]->points->x(i,0);
        }
        m1_axi = job->bodies[solid_body_id]->S * B_tmp;

    } else {
        //otherwise integration area and volume are the same
        v_i = body->S * body->points->v;
    }

    for (int i=0; i<V_i.rows(); i++){
        tmpNum = (v_i(i) - V_i(i))/(V_i(i));
        e(i) = std::max(0.0,tmpNum);
    }

    //calculate gradient of nodal overshoot
    grad_e = body->gradS.operate(e, MPMSparseMatrixBase::TRANSPOSED);

    Eigen::VectorXd J_tr = J;

    /*----------------------------------------------------------------------------*/
    for (int i=0;i<n.rows();i++){

        m1 = job->bodies[solid_body_id]->nodes->m[i];
        m2 = body->nodes->m[i];

        // set ((1-n)*vs + n*vw)
        if (m1 > 0 && m2 > 0) {
            if (job->JOB_TYPE == job->JOB_AXISYM){
                n(i) = 1 - (m1_axi(i) / (job->grid->nodeVolume(job,i) * solid_rho));
            } else {
                n(i) = 1 - (m1 / (job->grid->nodeVolume(job,i)*solid_rho));
            }
            mv1i = job->bodies[solid_body_id]->nodes->mx_t(i);
            mv2i = body->nodes->mx_t(i);

            if (n(i) < 0.2){
                n(i) = 0.2; //keep packing from overestimates...
            }

            nMat(i) = (1-n(i))/m1 * mv1i + n(i)/m2 * mv2i;
        } else if (m2 > 0) {
            n(i) = 1.0;
            mv2i = body->nodes->mx_t(i);

            nMat(i) = 1.0/m2 * mv2i;
        } else {
            n(i) = 1.0;
            nMat(i).setZero();
        }
    }

    //interpolate n to points
    n_p  = body->S.operate(n,MPMSparseMatrixBase::TRANSPOSED);

    //calculate divergence of ((1-n)*vs + n*vw)
    pvec = body->gradS.dot(nMat, MPMSparseMatrixBase::TRANSPOSED);

    //adjust pvec for axisym case
    if (job->JOB_TYPE == job->JOB_AXISYM){
        KinematicVectorArray tmpVec(body->points->x.size(),body->points->x.VECTOR_TYPE);
        tmpVec = body->S.operate(nMat, MPMSparseMatrixBase::TRANSPOSED);
        for (int p=0; p<pvec.rows(); p++){
            pvec(p) += tmpVec(p,0)/body->points->x(p,0);
        }
    }

    MaterialTensor tmpMat;
    double eta; //fluid viscosity change w/ phi

    for (int i=0;i<body->points->x.size();i++){
        if (body->points->active[i] == 0){
            continue;
        }

        L = body->points->L(i);

        //n*rho_dot/rho = -div((1-n)*v_s + n*v_w)
        //n*v_dot/v = div((1-n)*v_s + n*v_w)
        if (n_p(i) > 0 && n_p(i) < 1.0){
            J_tr(i) *= std::exp(job->dt * pvec(i) / n_p(i));
        } else {
            J_tr(i) *= std::exp(job->dt * L.trace());
        }

        D = 0.5*(L+L.transpose());

        trD = D.trace();

        //T = 2*mu*D_0 + K*log(J)*I
        eta = eta_0 * (1 + 5.0/2.0 * (1 - n_p(i)));
        T = 2*eta*(D - (trD/3.0)*MaterialTensor::Identity()) + K*std::log(J_tr(i))*MaterialTensor::Identity();

        body->points->T(i) = T;

        //position correction
        if (SPEC == Material::UPDATE && alpha > 0) {
            delta = -alpha * job->dt * (L - (trD / 3.0) * MaterialTensor::Identity()).norm() *
                    grad_e(i) * h * h;

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
    }

    if (SPEC == Material::UPDATE){
        J = J_tr;
    }

    return;
}

void SlurryFluidPhase::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    body->points->T(idIN) = stressIN;
    J(idIN) = std::exp(stressIN.trace()/(K*3.0)); //p = -K*log(J)
    return;
}

void SlurryFluidPhase::assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    MaterialTensor tmpMat;
    double trT;

    tmpMat = body->points->T(idIN);

    trT = tmpMat.trace();
    tmpMat = tmpMat - (trT/3.0 + pressureIN) * MaterialTensor::Identity();
    body->points->T(idIN) = tmpMat;

    J(idIN) = std::exp(-pressureIN/K); //p = -K*log(J)

    return;
}

/*----------------------------------------------------------------------------*/

void SlurryFluidPhase::writeFrame(Job* job, Body* body, Serializer* serializer){
    serializer->writeScalarArray(J,"J");
    serializer->writeScalarArray(n_p,"porosity");
    serializer->writeVectorArray(grad_e,"grad_err");
    serializer->writeVectorArray(del_pos,"del_pos");
    serializer->writeScalarArray(V_i, "grid_volume");
    serializer->writeScalarArray(e, "err");


    Eigen::VectorXd nvec(body->nodes->x.size());
    Eigen::VectorXd pvec(body->points->x.size());
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
    serializer->writeScalarArray(nvec, "nodal_pressure");
    return;
}

std::string SlurryFluidPhase::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int SlurryFluidPhase::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}
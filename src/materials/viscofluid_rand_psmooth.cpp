//
// Created by aaron on 12/8/17.
// viscofluid_rand_psmooth.cpp
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigenvalues>

#include "job.hpp"
#include "serializer.hpp"

#include "body.hpp"
#include "points.hpp"
#include "material.hpp"

double alpha = 1.0;
double mu, K;
Eigen::VectorXd Lx;
double h; //characteristic length of problem
Eigen::VectorXd delV;
Eigen::MatrixXd grad_delV;
Eigen::VectorXd V_i;
Eigen::VectorXd v_i;
Eigen::MatrixXd del_pos;

extern "C" void materialWriteFrame(Job* job, Body* body, Serializer* serializer);
extern "C" std::string materialSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
extern "C" int materialLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);

extern "C" void materialInit(Job* job, Body* body);
extern "C" void materialCalculateStress(Job* job, Body* body, int SPEC);
extern "C" void materialAssignStress(Job* job, Body* body, Eigen::MatrixXd stressIN, int idIN, int SPEC);
extern "C" void materialAssignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

/*----------------------------------------------------------------------------*/

void materialInit(Job* job, Body* body){
    if (body->material.fp64_props.size() < 3){
        std::cout << body->material.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 3 properties defined (K, mu, h).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        K = body->material.fp64_props[0];
        mu = body->material.fp64_props[1];
        h = body->material.fp64_props[2];
        printf("Material properties (K = %g, mu = %g, h = %g).\n",
               K, mu, h);
    }

    if (body->material.fp64_props.size() > 3){
        alpha = body->material.fp64_props[3];
    }

    V_i.resize(job->grid.node_count,1);
    v_i.resize(job->grid.node_count,1);
    delV.resize(job->grid.node_count,1);
    grad_delV = job->jobVectorArray<double>(body->points.x.rows());
    del_pos = job->jobVectorArray<double>(body->points.x.rows());
    del_pos.setZero();

    for (size_t i=0; i<job->grid.node_count;i++){
        V_i(i) = job->grid.gridNodalVolume(job,i);
    }

    Lx = body->nodes.x.colwise().maxCoeff();

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void materialWriteFrame(Job* job, Body* body, Serializer* serializer) {
    serializer->serializerWriteVectorArray(grad_delV,"grad_delV");
    serializer->serializerWriteVectorArray(del_pos,"del_pos");

    Eigen::VectorXd nvec(body->nodes.x.rows());
    Eigen::VectorXd pvec(body->points.x.rows());
    Eigen::VectorXd tmpVec;
    Eigen::MatrixXd T;
    for (size_t i=0;i<body->points.x.rows();i++) {
        if (body->points.active[i] == 0) {
            pvec(i) = 0;
            continue;
        }

        tmpVec = body->points.T.row(i).transpose();
        T = job->jobTensor<double>(tmpVec.data());

        pvec(i) = -T.trace() * body->points.v(i) / job->DIM;
    }
    body->bodyCalcNodalValues(job, nvec, pvec, Body::SET);
    body->bodyCalcNodalValues(job, v_i, body->points.v, Body::SET);
    nvec = nvec.array() / v_i.array();
    body->bodyCalcPointValues(job, pvec, nvec, Body::SET);


    serializer->serializerWriteScalarArray(pvec, "p_smooth");
    //nothing to report
    return;
}

/*----------------------------------------------------------------------------*/

std::string materialSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename
    std::stringstream s;
    s << "mpm_v2."  << body->name << "." << body->id << ".material." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    if (ffile.is_open()){
        ffile << "# mpm_v2 materials/isolin.so\n";
        ffile << K << "\n" << mu << "\n" << h << "\n" << alpha << "\n";
        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filepath+filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Material Saved: [" << body->name << "]." << std::endl;

    return filename;
}

/*----------------------------------------------------------------------------*/

int materialLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    std::string line;
    std::ifstream fin(fullpath);

    if(fin.is_open()){
        std::getline(fin,line); //first line
        std::getline(fin,line); //K
        K = std::stod(line);
        std::getline(fin,line); //mu
        mu = std::stod(line);
        std::getline(fin,line); //h
        h = std::stod(line);
        std::getline(fin,line);
        alpha = std::stod(line);

        fin.close();

        printf("Material properties (K = %g, mu = %g).\n",
               K, mu);
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    V_i.resize(job->grid.node_count,1);
    v_i.resize(job->grid.node_count,1);
    delV.resize(job->grid.node_count,1);
    grad_delV = job->jobVectorArray<double>(body->points.x.rows());
    del_pos = job->jobVectorArray<double>(body->points.x.rows());
    del_pos.setZero();

    for (size_t i=0; i<job->grid.node_count;i++){
        V_i(i) = job->grid.gridNodalVolume(job,i);
    }

    Lx = body->nodes.x.colwise().maxCoeff();

    std::cout << "Material Loaded: [" << body->name << "]." << std::endl;
    return 1;
}

/*----------------------------------------------------------------------------*/

void materialCalculateStress(Job* job, Body* body, int SPEC){
    Eigen::MatrixXd T = job->jobTensor<double>();
    Eigen::MatrixXd L = job->jobTensor<double>();
    Eigen::MatrixXd D = job->jobTensor<double>();
    Eigen::MatrixXd tmp = job->jobTensor<double>();

    Eigen::EigenSolver<Eigen::MatrixXd> es;

    double trD;
    double alpha = 1.0;
    double w = 0;
    int max_index = -1;
    int rnd_num;
    double tmpNum;
    Eigen::VectorXd e = job->jobVector<double>();

    //Eigen::VectorXd J = body->points.v.array() / body->points.v0.array();
    //double J;

    Eigen::MatrixXd tmpMat = job->jobTensor<double>();
    Eigen::VectorXd tmpVec;

    //calculate nodal volume integral
    body->bodyCalcNodalValues(job, v_i, body->points.v, Body::SET);
    for (size_t i=0; i<V_i.rows(); i++){
        tmpNum = (v_i(i) - V_i(i))/V_i(i);
        delV(i) = std::max(0.0,tmpNum);
    }

    //calculate gradient of nodal overshoot
    body->bodyCalcPointGradient(job, grad_delV, delV, Body::SET);

    //calculate nodal density
    //rho_i = body->nodes.m.array() / v_i.array();

    //calculate gradient of nodal density
    //body->bodyCalcPointGradient(job, grad_rho, rho_i, Body::SET);

    //strain smoothing
    Eigen::VectorXd J = body->points.v.array() / body->points.v0.array();
    /*
    Eigen::VectorXd pvec(body->points.x.rows());
    pvec = J.array().log();
    pvec = pvec.array()*body->points.v.array();
    Eigen::VectorXd nvec(body->nodes.x.rows());
    body->bodyCalcNodalValues(job, nvec, pvec, Body::SET);
    nvec = nvec.array() / v_i.array();
    body->bodyCalcPointValues(job, pvec, nvec, Body::SET);
    J = pvec.array().exp();
    body->points.v = body->points.v0.array() * J.array();
     */

    //p_i = -K * nvec;

    for (size_t i=0;i<body->points.x.rows();i++){
        if (body->points.active[i] == 0){
            continue;
        }

        tmpVec = body->points.L.row(i).transpose();
        L = job->jobTensor<double>(tmpVec.data());

        D = 0.5*(L+L.transpose());

        trD = D.trace();

        //stress update after density correction

        e = -alpha * job->dt * (L - (trD/D.rows())*job->jobTensor<double>(Job::IDENTITY)).norm() * grad_delV.row(i).transpose() * h * h;

        for (size_t pos=0;pos<e.rows();pos++){
            if (((body->points.x(i,pos) + e(pos)) >= Lx(pos)) || ((body->points.x(i,pos) + e(pos)) <= 0)){
                e(pos) = 0;
            }
        }

        body->points.x.row(i) += e.transpose();
        body->points.u.row(i) += e.transpose();
        del_pos.row(i) += e.transpose();

        body->points.x_t.row(i) += (L*e).transpose();
        body->points.mx_t.row(i) = body->points.m(i) * body->points.x_t.row(i);

        //rho = rho + e*grad_rho
        //m/vnew = m/v + e*grad_rho
        //vnew = 1/(1/v + e*grad_rho/m)
        //body->points.v(i) = body->points.m(i) / (body->points.m(i)/body->points.v(i) + e.dot(grad_rho.row(i).transpose()));
        //J = body->points.v(i)/body->points.v0(i);

        //T = 2*mu*D_0 + K*log(J)*I
        T = 2*mu*(D - (trD/D.rows())*job->jobTensor<double>(Job::IDENTITY)) + K*std::log(J(i))*job->jobTensor<double>(Job::IDENTITY);

        for (size_t i=0;i<tmpVec.size();i++){
            tmpVec(i) = T(i);
        }

        for (size_t pos=0;pos<T.size();pos++){
            body->points.T(i,pos) = T(pos);
        }

        /*
        es.compute(D - (trD/D.rows())*job->jobTensor<double>(Job::IDENTITY),true);

        //find max eigenvalue of D_0
        max_index = -1;
        w = 0;
        for (int pos=0;pos<D.cols();pos++){
            if (es.pseudoEigenvalueMatrix()(pos,pos) > w){
                w = es.pseudoEigenvalueMatrix()(pos,pos);
                max_index = pos;
            }
        }

        //find associated eigenvector
        if (max_index >= 0){
            tmp = es.pseudoEigenvectors();
            e = tmp.col(max_index);

            rnd_num = 2*(rand()%2) - 1;

            body->points.x.row(i) += w * job->dt * body->points.extent(i) * alpha * rnd_num * e.transpose();
            body->points.u.row(i) += w * job->dt * body->points.extent(i) * alpha * rnd_num * e.transpose();
        }
        */
    }

    return;
}

/*----------------------------------------------------------------------------*/

void materialAssignStress(Job* job, Body* body, Eigen::MatrixXd stressIN, int idIN, int SPEC){
    for (size_t pos=0;pos<stressIN.size();pos++){
        body->points.T(idIN,pos) = stressIN(pos);
    }
    return;
}

/*----------------------------------------------------------------------------*/

void materialAssignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    Eigen::MatrixXd tmpMat;
    Eigen::VectorXd tmpVec;
    double trT;

    tmpVec = body->points.T.row(idIN).transpose();
    tmpMat = job->jobTensor<double>(tmpVec.data());

    trT = tmpMat.trace();
    tmpMat = tmpMat - (trT/tmpMat.rows() + pressureIN) * job->jobTensor<double>(Job::IDENTITY);
    for (size_t pos=0;pos<tmpMat.size();pos++){
        body->points.T(idIN,pos) = tmpMat(pos);
    }

    body->points.v(idIN) = body->points.v0(idIN)*std::exp(-pressureIN/K); //p = -K*log(J)

    return;
}

/*----------------------------------------------------------------------------*/
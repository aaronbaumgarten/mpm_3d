//
// Created by aaron on 12/5/17.
// slurry_viscofluid_rand.cpp
//

//
// Created by aaron on 6/10/17.
// slurry_viscofluid.cpp
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

#include "job.hpp"
#include "serializer.hpp"

#include "body.hpp"
#include "points.hpp"
#include "material.hpp"

//double mu, K, solid_rho;
double mu = 8.9e-4;
double K = 1e9;
double solid_rho = 2500;
double h = 0;
double alpha = 1.0;

Eigen::VectorXd Lx;
Eigen::VectorXd e;
Eigen::MatrixXd grad_e;
Eigen::VectorXd V_i;
Eigen::VectorXd v_i;
Eigen::MatrixXd del_pos;

int solid_body_id = -1;
Eigen::VectorXd J(0,1);
Eigen::VectorXd n_p(0,1);

extern "C" void materialWriteFrame(Job* job, Body* body, Serializer* serializer);
extern "C" std::string materialSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
extern "C" int materialLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);

extern "C" void materialInit(Job* job, Body* body);
extern "C" void materialCalculateStress(Job* job, Body* body, int SPEC);
extern "C" void materialAssignStress(Job* job, Body* body, Eigen::MatrixXd stressIN, int idIN, int SPEC);
extern "C" void materialAssignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

/*----------------------------------------------------------------------------*/

void materialInit(Job* job, Body* body){
    if (body->material.fp64_props.size() < 4 || (body->material.str_props.size() < 1 && body->material.int_props.size() < 1)){
        std::cout << body->material.fp64_props.size() <<", " << body->material.str_props.size() << ", " << body->material.int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 5 properties defined ({K, mu, solid_rho, h}, {solid_body}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        K = body->material.fp64_props[0];
        mu = body->material.fp64_props[1];
        solid_rho = body->material.fp64_props[2];
        h = body->material.fp64_props[3];

        //set body ids by name
        if (body->material.str_props.size() == 1){
            for (size_t b = 0; b < job->bodies.size(); b++) {
                if (body->material.str_props[0].compare(job->bodies[b].name) == 0){
                    solid_body_id = b;
                    break;
                }
            }
        }

        // or set body ids by int
        if (solid_body_id < 0){
            if (body->material.int_props.size() == 1) {
                solid_body_id = body->material.int_props[0];
            } else {
                std::cout << body->material.fp64_props.size() <<", " << body->material.str_props.size() << ", " << body->material.int_props.size() << "\n";
                fprintf(stderr,
                        "%s:%s: Need at least 4 properties defined ({K, mu, solid_rho}, {solid_body}).\n",
                        __FILE__, __func__);
                exit(0);
            }
        }

        printf("Material properties ({K = %g, mu = %g, solid_rho = %g}, {solid_body: %i}).\n",
               K, mu, solid_rho, solid_body_id);
    }

    if (body->material.fp64_props.size() > 4){
        alpha = body->material.fp64_props[4];
    }

    J.resize(body->points.x.rows());
    J.setOnes();
    n_p.resize(body->points.x.rows());
    n_p.setOnes();

    V_i.resize(job->grid.node_count,1);
    v_i.resize(job->grid.node_count,1);
    e.resize(job->grid.node_count,1);
    grad_e = job->jobVectorArray<double>(body->points.x.rows());
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
    serializer->serializerWriteScalarArray(J,"J");
    serializer->serializerWriteScalarArray(n_p,"porosity");
    serializer->serializerWriteVectorArray(grad_e,"grad_err");
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

    serializer->serializerWriteScalarArray(nvec, "nodal_pressure");
    serializer->serializerWriteScalarArray(pvec, "p_smooth");
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
        ffile << K << "\n" << mu << "\n" << solid_rho << "\n" << h << "\n" << solid_body_id << "\n" << alpha << "\n";
        ffile << J.rows() << "\n";
        ffile << "J\n";
        ffile << "{\n";
        job->jobScalarArrayToFile(J, ffile);
        ffile << "}";
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
    int len = 0;

    if(fin.is_open()){
        std::getline(fin,line); //first line
        std::getline(fin,line); //K
        K = std::stod(line);
        std::getline(fin,line); //mu
        mu = std::stod(line);
        std::getline(fin,line); //solid_rho
        solid_rho = std::stod(line);
        std::getline(fin,line);//h
        h = std::stod(line);
        std::getline(fin,line); //solid_body_id
        solid_body_id = std::stoi(line);
        std::getline(fin,line);
        alpha = std::stod(line);

        std::getline(fin,line); //len
        len = std::stoi(line);
        J.resize(len);

        std::getline(fin, line); //J
        std::getline(fin, line); //{
        job->jobScalarArrayFromFile(J, fin);
        std::getline(fin, line); //}

        fin.close();

        printf("Material properties ({K = %g, mu = %g, solid_rho = %g}, {solid_body: %i}).\n",
               K, mu, solid_rho, solid_body_id);
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    n_p.resize(body->points.x.rows());
    n_p.setOnes();

    V_i.resize(job->grid.node_count,1);
    v_i.resize(job->grid.node_count,1);
    e.resize(job->grid.node_count,1);
    grad_e = job->jobVectorArray<double>(body->points.x.rows());
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

    double m1, m2;
    Eigen::VectorXd mv1i = job->jobVector<double>();
    Eigen::VectorXd mv2i = job->jobVector<double>();

    Eigen::VectorXd n(body->nodes.x.rows());
    Eigen::VectorXd pvec(body->points.x.rows());
    //Eigen::VectorXd n_p(body->points.x.rows());
    Eigen::MatrixXd nMat = job->jobVectorArray<double>(body->nodes.x.rows());

    double trD, tmpNum;

    //perturbation variables
    Eigen::VectorXd delta = job->jobVector<double>();

    //calculate nodal volume integral
    body->bodyCalcNodalValues(job, v_i, body->points.v, Body::SET);
    for (size_t i=0; i<V_i.rows(); i++){
        tmpNum = (v_i(i) - V_i(i))/(V_i(i));
        e(i) = std::max(0.0,tmpNum);
    }

    //calculate gradient of nodal overshoot
    body->bodyCalcPointGradient(job, grad_e, e, Body::SET);

    Eigen::VectorXd J_tr = J;

    /*----------------------------------------------------------------------------*/
    for (size_t i=0;i<n.rows();i++){

        m1 = job->bodies[solid_body_id].nodes.m[i];
        m2 = body->nodes.m[i];

        // set ((1-n)*vs + n*vw)
        if (m1 > 0 && m2 > 0) {
            n(i) = 1 - (m1 / (job->grid.gridNodalVolume(job,i)*solid_rho));
            mv1i = job->bodies[solid_body_id].nodes.mx_t.row(i).transpose();
            mv2i = body->nodes.mx_t.row(i).transpose();

            if (n(i) < 0.2){
                n(i) = 0.2; //keep packing from overestimates...
            }

            nMat.row(i) = (1-n(i))/m1 * mv1i.transpose() + n(i)/m2 * mv2i.transpose();
        } else if (m2 > 0) {
            n(i) = 1.0;
            mv2i = body->nodes.mx_t.row(i).transpose();

            nMat.row(i) = 1.0/m2 * mv2i.transpose();
        } else {
            n(i) = 1.0;
            nMat.row(i).setZero();
        }
    }

    //interpolate n to points
    body->bodyCalcPointValues<Eigen::VectorXd,Eigen::VectorXd>(job,n_p,n,Body::SET);

    //calculate divergence of ((1-n)*vs + n*vw)
    body->bodyCalcPointDivergence(job,pvec,nMat,Body::SET);

    Eigen::MatrixXd tmpMat = job->jobTensor<double>();
    Eigen::VectorXd tmpVec;
    double eta; //fluid viscosity change w/ phi

    for (size_t i=0;i<body->points.x.rows();i++){
        if (body->points.active[i] == 0){
            continue;
        }

        tmpVec = body->points.L.row(i).transpose();
        L = job->jobTensor<double>(tmpVec.data());

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
        eta = mu * (1 + 5.0/2.0 * (1 - n_p(i)));
        T = 2*eta*(D - (trD/D.rows())*job->jobTensor<double>(Job::IDENTITY)) + K*std::log(J_tr(i))*job->jobTensor<double>(Job::IDENTITY);

        for (size_t i=0;i<tmpVec.size();i++){
            tmpVec(i) = T(i);
        }

        for (size_t pos=0;pos<T.size();pos++){
            body->points.T(i,pos) = T(pos);
        }

        //position correction
        if (SPEC == Material::UPDATE) {
            delta = -alpha * job->dt * (L - (trD / D.rows()) * job->jobTensor<double>(Job::IDENTITY)).norm() *
                grad_e.row(i).transpose() * h * h;

            for (size_t pos = 0; pos < delta.rows(); pos++) {
                if (((body->points.x(i, pos) + delta(pos)) >= Lx(pos)) || ((body->points.x(i, pos) + delta(pos)) <= 0)) {
                    delta(pos) = 0;
                }
            }

            body->points.x.row(i) += delta.transpose();
            body->points.u.row(i) += delta.transpose();
            del_pos.row(i) += delta.transpose();

            body->points.x_t.row(i) += (L * delta).transpose();
            body->points.mx_t.row(i) = body->points.m(i) * body->points.x_t.row(i);
        }
    }

    if (SPEC == Material::UPDATE){
        J = J_tr;
    }

    return;
}

/*----------------------------------------------------------------------------*/

void materialAssignStress(Job* job, Body* body, Eigen::MatrixXd stressIN, int idIN, int SPEC){
    for (size_t pos=0;pos<stressIN.size();pos++){
        body->points.T(idIN,pos) = stressIN(pos);
    }
    J(idIN) = std::exp(-stressIN.trace()/(K*stressIN.rows())); //p = K*log(J)
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

    J(idIN) = std::exp(-pressureIN/K); //p = -K*log(J)

    return;
}

/*----------------------------------------------------------------------------*/
//
// Created by aaron on 6/6/17.
// g_local_mu2.cpp
//


//
// Created by aaron on 5/24/17.
// isolin.cpp
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "job.hpp"
#include "serializer.hpp"

#include "body.hpp"
#include "points.hpp"
#include "material.hpp"

Eigen::VectorXd gammap(0,1);
Eigen::VectorXd gammap_dot(0,1);

double E, nu, G, K;
double lambda;
double grains_d;

double MU_S = 0.280;
double GRAINS_RHO = 2500;//2450
double RHO_CRITICAL = (GRAINS_RHO*0.6);//1500
double MU_2 = MU_S;//MU_S//(I_0+MU_S)
double DELTA_MU = (MU_2 - MU_S);
double I_0 = 0.278;

extern "C" void materialWriteFrame(Job* job, Body* body, Serializer* serializer);
extern "C" std::string materialSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
extern "C" int materialLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);

extern "C" void materialInit(Job* job, Body* body);
extern "C" void materialCalculateStress(Job* job, Body* body, int SPEC);
extern "C" void materialAssignStress(Job* job, Body* body, Eigen::MatrixXd stressIN, int idIN, int SPEC);
extern "C" void materialAssignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

/*----------------------------------------------------------------------------*/

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

/*----------------------------------------------------------------------------*/

void materialInit(Job* job, Body* body){
    if (body->material.fp64_props.size() < 6) {
        // Bit of a hack, but it's okay for now. Just close your eyes and code it anyways.
        fprintf(stderr,
                "%s:%s: Need at least 6 properties defined (E, nu, grain diameter, mu_s, grains_rho, rho_critical).\n",
                __FILE__, __func__);
        exit(EXIT_FAILURE);
    } else if (body->material.fp64_props.size() < 8){
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        grains_d = body->material.fp64_props[2];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2*nu));
        lambda = K - 2.0 * G / 3.0;

        MU_S = body->material.fp64_props[3];
        GRAINS_RHO = body->material.fp64_props[4];
        RHO_CRITICAL = body->material.fp64_props[5];

        MU_2 = MU_S;
        DELTA_MU = MU_2 - MU_S;

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, grain diameter = %g, mu_s = %g, grains_rho = %g, rho_critical = %g).\n",
               E, nu, G, K, grains_d, MU_S, GRAINS_RHO, RHO_CRITICAL);
    } else {
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        grains_d = body->material.fp64_props[2];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2*nu));
        lambda = K - 2.0 * G / 3.0;

        MU_S = body->material.fp64_props[3];
        GRAINS_RHO = body->material.fp64_props[4];
        RHO_CRITICAL = body->material.fp64_props[5];

        MU_2 = body->material.fp64_props[6];
        DELTA_MU = MU_2 - MU_S;
        I_0 = body->material.fp64_props[7];

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, grain diameter = %g, mu_s = %g, grains_rho = %g, rho_critical = %g, mu_2 = %g, i_0 = %g).\n",
               E, nu, G, K, grains_d, MU_S, GRAINS_RHO, RHO_CRITICAL, MU_2, I_0);
    }

    gammap_dot.resize(body->points.x.rows());
    gammap_dot.setZero();
    gammap.resize(body->points.x.rows());
    gammap.setZero();

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void materialWriteFrame(Job* job, Body* body, Serializer* serializer) {
    serializer->serializerWriteScalarArray(gammap,"gammap");
    serializer->serializerWriteScalarArray(gammap_dot,"gammap_dot");
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
    std::ostringstream s;
    s << "mpm_v2."  << body->name << "." << body->id << ".material." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    if (ffile.is_open()){
        ffile << "# mpm_v2 materials/isolin.so\n";
        ffile << E << "\n" << nu << "\n" << gammap.rows() << "\n";
        ffile << "gammap_dot\n{\n";
        job->jobScalarArrayToFile(gammap_dot,ffile);
        ffile << "}\n";
        ffile << "gammap\n{\n";
        job->jobScalarArrayToFile(gammap,ffile);
        ffile << "}\n";
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
        std::getline(fin,line); //E
        E = std::stod(line);
        std::getline(fin,line); //nu
        nu = std::stod(line);
        std::getline(fin,line); //len
        int len = std::stoi(line);

        gammap.resize(len);
        gammap.setZero();
        gammap_dot.resize(len);
        gammap_dot.setZero();

        while(std::getline(fin,line)){
            if (line.compare("gammap_dot") == 0){
                std::getline(fin,line); //{
                job->jobScalarArrayFromFile(gammap_dot,fin);
                std::getline(fin,line); //}
            } else if (line.compare("gammap") == 0){
                std::getline(fin,line);
                job->jobScalarArrayFromFile(gammap,fin);
                std::getline(fin,line);
            } else {
                std::cerr << "Unknown field title: " << line << std::endl;
            }
        }

        fin.close();

        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g).\n",
               E, nu, G, K);
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Material Loaded: [" << body->name << "]." << std::endl;
    return 1;
}

/*----------------------------------------------------------------------------*/

void materialCalculateStress(Job* job, Body* body, int SPEC){
    /*
        cohesion, but other values should probably be changed to make this
        nonzero.
    *///////////////////////////////////////////////////////////////////////
    double c = 0;

    Eigen::MatrixXd T = job->jobTensor<double>();
    Eigen::MatrixXd L = job->jobTensor<double>();
    Eigen::MatrixXd D = job->jobTensor<double>();
    Eigen::MatrixXd W = job->jobTensor<double>();

    Eigen::MatrixXd s_tr = job->jobTensor<double>();
    Eigen::MatrixXd t0_tr = job->jobTensor<double>();

    double trD, p_tr, tau_tr, nup_tau;
    double S0, tau_tau, scale_factor;
    double S2, alpha, B, H;

    bool density_flag;

    Eigen::MatrixXd tmpMat = job->jobTensor<double>();
    Eigen::VectorXd tmpVec;

    for (size_t i=0;i<body->points.x.rows();i++){
        if (body->points.active[i] == 0){
            continue;
        }

        tmpVec = body->points.L.row(i).transpose();
        L = job->jobTensor<double>(tmpVec.data());

        tmpVec = body->points.T.row(i).transpose();
        T = job->jobTensor<double>(tmpVec.data());

        D = 0.5*(L+L.transpose());
        W = 0.5*(L-L.transpose());

        trD = D.trace();

        tmpMat = (2*G*D) + (lambda*trD*job->jobTensor<double>(Job::IDENTITY)) + (W*T) - (T*W);

        //trial stress
        s_tr = T + job->dt * tmpMat;
        p_tr = -s_tr.trace()/s_tr.rows();

        //trial deviator
        t0_tr = s_tr + p_tr * job->jobTensor<double>(Job::IDENTITY);
        tau_tr = t0_tr.norm()/std::sqrt(2.0);

        density_flag = false;
        if ((body->points.m(i) / body->points.v(i)) < RHO_CRITICAL){
            density_flag = true;
        }

        if (density_flag || p_tr <= c){
            nup_tau = tau_tr / (G * job->dt);
            body->points.T.row(i).setZero();
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

            tmpMat = scale_factor * t0_tr - p_tr * job->jobTensor<double>(Job::IDENTITY);
            for (size_t i=0;i<tmpVec.size();i++){
                tmpVec(i) = tmpMat(i);
            }
            body->points.T.row(i) = tmpVec.transpose();

            /*for (size_t pos=0;pos<9;pos++){
                body->particles.T(i,pos) = scale_factor * t0_tr(pos);
                if (pos == XX || pos == YY || pos == ZZ){
                    body->particles.T(i,pos) -= p_tr;
                }
            }*/

        } else {
            if (SPEC == Material::UPDATE) {
                std::cerr << "u\n";
            }
            nup_tau = 0;
        }

        gammap(i) += nup_tau * job->dt;
        gammap_dot(i) = nup_tau;
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
    return;
}

/*----------------------------------------------------------------------------*/
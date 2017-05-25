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
#include <eigen3/Eigen/Core>

#include "job.hpp"
#include "serializer.hpp"

#include "body.hpp"
#include "points.hpp"
#include "material.hpp"

double E, nu, G, K;
double lambda;

extern "C" void materialWriteFrame(Job* job, Body* body, Serializer* serializer);
extern "C" std::string materialSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
extern "C" int materialLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);

extern "C" void materialInit(Job* job, Body* body);
extern "C" void materialCalculateStress(Job* job, Body* body, int update = 1);
extern "C" void materialAssignStress(Job* job, Body* body, Eigen::MatrixXd stressIN, int idIN);
extern "C" void materialAssignPressure(Job* job, Body* body, double pressureIN, int idIN);

/*----------------------------------------------------------------------------*/

void materialWriteFrame(Job* job, Body* body, Serializer* serializer) {
    //nothing to report
    return;
}

/*----------------------------------------------------------------------------*/

std::string materialSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);

    //create filename
    std::ostringstream s;
    s << "mpm_v2."  << body->name << "." << body->id << ".material." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    std::string filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    if (ffile.is_open()){
        ffile << "# mpm_v2 materials/isolin.so\n";
        ffile << E << "\n" << nu << "\n";
        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filepath+filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Material Saved: [" << body->id << "]." << std::endl;

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

    std::cout << "Material Loaded: [" << body->id << "]." << std::endl;
    return 1;
}

/*----------------------------------------------------------------------------*/

void materialInit(Job* job, Body* body){
    if (body->material.fp64_props.size() < 2){
        std::cout << body->material.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined (E, nu).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        printf("Material properties (E = %g, nu = %g, G = %g, K = %g).\n",
               E, nu, G, K);
    }

    std::cout << "Material Initialized: [" << body->id << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void materialCalculateStress(Job* job, Body* body, int update=1){
    Eigen::MatrixXd T = job->jobTensor<double>();
    Eigen::MatrixXd L = job->jobTensor<double>();
    Eigen::MatrixXd D = job->jobTensor<double>();
    Eigen::MatrixXd W = job->jobTensor<double>();

    double trD;

    Eigen::Matrix tmpMat = job->jobTensor<double>();
    Eigen::Matrix tmpVec = job->jobVector<double>();

    for (size_t i=0;i<body->points.x.rows();i++){
        if (body->points.active[i] == 0){
            continue;
        }

        tmpVec << body->points.L.row(i).transpose();
        L = job->jobTensor<double>(tmpVec.data());

        tmpVec << body->points.T.row(i).transpose();
        T = job->jobTensor<double>(tmpVec.data());

        D = 0.5*(L+L.transpose());
        W = 0.5*(L-L.transpose());

        trD = D.trace();

        tmpMat = (2*G*D) + (lambda*trD*job->jobTensor<double>(Job::IDENTITY)) + (W*T) - (T*W);
        tmpVec = job->jobVector<double>(tmpMat.data());

        body->points.T.row(i) = body->points.T.row(i) + job->dt * tmpVec.transpose();
    }

    return;
}

/*----------------------------------------------------------------------------*/

void materialAssignStress(Job* job, Body* body, Eigen::MatrixXd stressIN, int idIN){
    Eigen::Matrix tmpVec = job->jobVector<double>(stressIN.data());
    body->points.T.row(idIN) = tmpVec.transpose();
    return;
}

/*----------------------------------------------------------------------------*/

void materialAssignPressure(Job* job, Body* body, double pressureIN, int idIN){
    Eigen::Matrix tmpMat;
    Eigen::Matrix tmpVec = job->jobVector<double>();
    double trT;

    tmpVec << body->points.T.row(idIN).transpose();
    tmpMat = job->jobTensor<double>(tmpVec.data());

    trT = tmpMat.trace();
    tmpMat = tmpMat - (trT/3.0 + pressureIN) * job->jobTensor<double>(Job::IDENTITY);
    tmpVec = job->jobVector<double>(tmpMat.data());

    body->points.T.row(idIN) = tmpVec.transpose();
    return;
}

/*----------------------------------------------------------------------------*/
//
// Created by aaron on 6/20/17.
// rigid_linear_motion.cpp
//

//
// Created by aaron on 5/26/17.
// cartesian_box.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <grid.hpp>

#include "job.hpp"

#include "serializer.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "boundary.hpp"

Eigen::VectorXd velocity_constraint(0,1);

extern "C" void boundaryWriteFrame(Job* job, Body* body, Serializer* serializer); //write frame to serializer
extern "C" std::string boundarySaveState(Job* job, Body* body, Serializer* serializer, std::string filepath); //save state to serializer folder with returned filename
extern "C" int boundaryLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath); //read state from given full path

extern "C" void boundaryInit(Job* job, Body* body); //initialize boundary object
extern "C" void boundaryGenerateRules(Job* job, Body* body); //generate the rules given job and body state
extern "C" void boundaryApplyRules(Job* job, Body* body); //apply the rules given job and body state

/*----------------------------------------------------------------------------*/

void boundaryInit(Job* job, Body* body){
    velocity_constraint = job->jobVector<double>();

    if (body->boundary.fp64_props.size() < job->DIM){
        std::cout << body->boundary.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i properties defined (<linear_velocity>).\n",
                __FILE__, __func__,job->DIM);
        exit(0);
    } else {
        for (size_t i=0;i<job->DIM;i++){
            velocity_constraint(i) = body->boundary.fp64_props[i];
            body->points.x_t.col(i).setOnes();
            body->points.x_t.col(i) *= velocity_constraint(i);
            body->points.mx_t.col(i) = body->points.m * velocity_constraint(i);
        }
        printf("Boundary properties (linear_speed = %g).\n",
               velocity_constraint.norm());
    }

    std::cout << "Boundary Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void boundaryWriteFrame(Job* job, Body* body, Serializer* serializer){
    //write nodal mask to frame output
    Eigen::MatrixXd tmpMat = job->jobVectorArray<double>(body->points.x.rows());
    tmpMat.setOnes();
    for (size_t pos=0; pos < tmpMat.cols(); pos++){
        tmpMat.col(pos) *= velocity_constraint(pos);
    }
    serializer->serializerWriteVectorArray(tmpMat,("bc_constained_velocity_"+body->name));
    return;
}

/*----------------------------------------------------------------------------*/

std::string boundarySaveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename
    std::stringstream s;
    s << "mpm_v2."  << body->name << "." << body->id << ".boundary." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    //write data
    if (ffile.is_open()) {
        ffile << "# mpm_v2 boundaries/cartesian_box.so\n"; //header
        for (size_t i=0;i<velocity_constraint.rows();i++){
            ffile << velocity_constraint(i) << "\n";
        }

        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Boundary Saved: [" << body->name << "]." << std::endl;

    return filename;
}

/*----------------------------------------------------------------------------*/

int boundaryLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    std::string line; //read line

    std::ifstream fin(fullpath); //file to load from

    if (fin.is_open()) {
        std::getline(fin,line); //first line (header)
        velocity_constraint = job->jobVector<double>();
        for (size_t i=0;i<velocity_constraint.rows();i++){
            std::getline(fin,line);
            velocity_constraint(i) = std::stod(line);
        }
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Boundary Loaded: [" << body->name << "]." << std::endl;

    return 1;
}

/*----------------------------------------------------------------------------*/

void boundaryGenerateRules(Job* job, Body* body){
    //nothing to do here
    //bc mask has already been set
    return;
}

/*----------------------------------------------------------------------------*/

void boundaryApplyRules(Job* job, Body* body){
    for (size_t i=0;i<body->nodes.x_t.rows();i++){
        for (size_t pos=0;pos<body->nodes.x_t.cols();pos++){
            //zero out velocity on boundary
            body->nodes.x_t(i,pos) = velocity_constraint(pos);
            body->nodes.mx_t(i,pos) = body->nodes.m(i) * velocity_constraint(pos);
            body->nodes.f(i,pos) = 0;
        }
    }
    for (size_t pos=0;pos<job->DIM;pos++){
        body->points.x_t.col(pos).setOnes();
        body->points.x_t.col(pos) *= velocity_constraint(pos);
        body->points.mx_t.col(pos) = body->points.m * velocity_constraint(pos);
    }
    return;
}

/*----------------------------------------------------------------------------*/


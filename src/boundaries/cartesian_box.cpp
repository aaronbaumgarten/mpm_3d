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

Eigen::MatrixXi bcNodalMask; //store which dof to control

extern "C" void boundaryWriteFrame(Job* job, Body* body, Serializer* serializer); //write frame to serializer
extern "C" std::string boundarySaveState(Job* job, Body* body, Serializer* serializer, std::string filepath); //save state to serializer folder with returned filename
extern "C" int boundaryLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath); //read state from given full path

extern "C" void boundaryInit(Job* job, Body* body); //initialize boundary object
extern "C" void boundaryGenerateRules(Job* job, Body* body); //generate the rules given job and body state
extern "C" void boundaryApplyRules(Job* job, Body* body); //apply the rules given job and body state

/*----------------------------------------------------------------------------*/

void boundaryWriteFrame(Job* job, Body* body, Serializer* serializer){
    //write nodal mask to frame output
    serializer->serializerWriteVectorArray(bcNodalMask,"bc_nodal_mask");
    return;
}

/*----------------------------------------------------------------------------*/

std::string boundarySaveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);

    //create filename
    std::ostringstream s;
    s << "mpm_v2."  << body->name << "." << body->id << ".boundary." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    std::string filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    //write data
    if (ffile.is_open()) {
        ffile << "# mpm_v2 boundaries/cartesian_box.so\n"; //header
        ffile << bcNodalMask.rows(); //lines to read later
        job->jobVectorArrayToFile(bcNodalMask, ffile);

        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filename << "\" !\n";
        return "ERR";
    }
    return filename;
}

/*----------------------------------------------------------------------------*/

int boundaryLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    std::string line; //read line

    std::ifstream fin(fullpath); //file to load from

    if (fin.is_open()) {
        std::getline(fin,line); //first line (header)
        std::getline(fin,line); //length of file to be read
        bcNodalMask = job->jobVectorArray(std::stoi(line)); //initialize vector array
        job->jobVectorArrayFromFile(bcNodalMask, fin); //read in from file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Boundary Loaded: [" << body->id << "]." << std::endl;

    return 1;
}

/*----------------------------------------------------------------------------*/

void boundaryInit(Job* job, Body* body){
    if (job->grid.filename.compare("cartesian.so") != 0){
        std::cout << "\nBOUNDARY CONDITION WARNING!" << std::endl;
        std::cout << "\"cartesian_box.so\" boundary expects \"cartesian.so\" grid NOT \"" << job->grid.filename << "\"!\n" << std::endl;
    }

    //find bounds of box
    Eigen::VectorXd Lx;
    Lx << body->nodes.x.colwise().maxCoeff();

    //set bounding mask
    double len = body->nodes.x.rows();
    bcNodalMask = job->jobVectorArray(len);
    bcNodalMask.setZero();
    for (size_t i=0;i<len;i++){
        for (size_t pos=0;pos<bcNodalMask.cols();pos++){
            if (body->nodes.x(i,pos) == 0 || body->nodes.x(i,pos) == Lx(pos)) {
                bcNodalMask.row(i).setOnes();
                break;
            }
        }
    }

    return;
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
            if (bcNodalMask(i,pos) == 1){
                //zero out velocity on boundary
                body->nodes.x_t(i,pos) == 0;
            }
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/


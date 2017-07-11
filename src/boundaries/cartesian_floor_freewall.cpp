//
// Created by aaron on 7/7/17.
// cartesian_floor_freewall.cpp
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

void boundaryInit(Job* job, Body* body){
    if (job->grid.filename.compare("cartesian.so") != 0){
        std::cout << "\nBOUNDARY CONDITION WARNING!" << std::endl;
        std::cout << "\"cartesian_floor_freewall.so\" boundary expects \"cartesian.so\" grid NOT \"" << job->grid.filename << "\"!\n" << std::endl;
    }

    //find bounds of box
    Eigen::VectorXd Lx = body->nodes.x.colwise().maxCoeff();

    //set bounding mask
    double len = body->nodes.x.rows();
    bcNodalMask = job->jobVectorArray<int>(len);
    bcNodalMask.setZero();

    for (size_t i=0;i<len;i++){
        for (size_t pos=0;pos<body->nodes.x.cols();pos++){
            if (body->nodes.x(i,pos) == 0 && pos == (body->nodes.x.cols()-1)){
                bcNodalMask.row(i).setOnes(); //lock all directions for floor
            }

            if (body->nodes.x(i,pos) == 0 || body->nodes.x(i,pos) == Lx(pos)) {
                bcNodalMask(i,pos) = 1;
            }
        }

    }

    std::cout << "Boundary Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void boundaryWriteFrame(Job* job, Body* body, Serializer* serializer){
    //write nodal mask to frame output
    Eigen::MatrixXd tmpMat = bcNodalMask.cast<double>();
    serializer->serializerWriteVectorArray(tmpMat,"bc_nodal_mask");
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
        ffile << bcNodalMask.rows() << "\n"; //lines to read later
        job->jobVectorArrayToFile(bcNodalMask, ffile);

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
        std::getline(fin,line); //length of file to be read
        bcNodalMask = job->jobVectorArray<int>(std::stoi(line)); //initialize vector array
        job->jobVectorArrayFromFile(bcNodalMask, fin); //read in from file
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
            if (bcNodalMask(i,pos) == 1){
                //zero out velocity on boundary
                body->nodes.x_t(i,pos) = 0;
                body->nodes.mx_t(i,pos) = 0;
                body->nodes.f(i,pos) = 0;
            }
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/


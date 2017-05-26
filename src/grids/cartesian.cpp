//
// Created by aaron on 5/24/17.
// cartesian.cpp
//


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "job.hpp"
#include "serializer.hpp"

#include "grid.hpp"

#include "body.hpp"
#include "nodes.hpp"

Eigen::VectorXd Lx;
Eigen::VectorXi Nx;
Eigen::VectorXd hx;
Eigen::MatrixXi nodeIDs; //element to node map
Eigen::MatrixXi A; //0,1 for directions in
size_t npe; //nodes per element
Eigen::MatrixXd x_n;

extern "C" void gridInit(Job* job);

extern "C" void gridWriteFrame(Job* job, Serializer* serializer);
extern "C" std::string gridSaveState(Job* job, Serializer* serializer,std::string filepath);
extern "C" int gridLoadState(Job* job, Serializer* serializer, std::string fullpath);

extern "C" int gridWhichElement(Job* job, Eigen::VectorXd xIN);
extern "C" Eigen::VectorXd gridNodeIDToPosition(Job* job, int idIN);
extern "C" void gridEvaluateShapeFnValue(Job* job, Eigen::VectorXd xIN, std::vector<int>& nID, std::vector<double>& nVAL);
extern "C" void gridEvaluateShapeFnGradient(Job* job, Eigen::VectorXd xIN, std::vector<int>& nID, std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& nGRAD);

/*----------------------------------------------------------------------------*/

void hiddenInit(Job* job){
    //called from loading or initializing functions
    //needs to have Lx,Nx,hx defined

    //count nodes
    job->grid.node_count = 1;
    job->grid.element_count = 1;
    npe = 1;
    for (size_t i=0;i<Nx.rows();i++){
        job->grid.node_count *= Nx(i);
        job->grid.element_count *= (Nx(i)-1);
        npe *= 2;
    }

    x_n = job->jobVectorArray<double>(job->grid.node_count);
    for (size_t i=0;i<x_n.rows();i++){
        x_n.row(i) = job->grid.gridNodeIDToPosition(job,i);
    }

    //initialize A matrix
    //for mapping position in cube to id
    //0 -> +0,+0,+0
    //1 -> +1,+0,+0
    //...
    //8 -> +1,+1,+1
    Eigen::VectorXi onoff = job->jobVector<int>();
    A = job->jobVectorArray<int>(npe);
    onoff.setZero();
    for (size_t n=0; n<npe;n++){
        for (size_t i=0;i<onoff.rows();i++){
            A(n,i) = onoff(i);
        }
        for (size_t i=0;i<onoff.rows();i++) {
            if (onoff(i) == 0){
                onoff(i) == 1;
                break;
            } else {
                onoff(i) == 0;
            }
        }
    }

    //setup element to node map
    Eigen::VectorXi ijk = job->jobVector<int>();
    int tmp;
    nodeIDs.resize(job->grid.element_count,npe);
    for (size_t e=0;e<nodeIDs.rows();e++){
        tmp = e;
        //find i,j,k count for element position
        for (size_t i=0;i<ijk.rows();i++){
            ijk(i) = tmp % (Nx(i)-1);
            tmp = tmp / (Nx(i)-1);
        }

        //find node ids for element
        for (size_t n=0;n<nodeIDs.cols();e++){
            for (size_t i=0;i<ijk.rows();i++) {
                nodeIDs(e, n) += (ijk(i)+A(n,i)) * Nx(i);
            }

        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void gridInit(Job* job){
    if (job->grid.fp64_props.size() < job->DIM || job->grid.fp64_props.size() < job->DIM){
        std::cout << job->grid.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i dimensions defined.\n",
                __FILE__, __func__, job->DIM);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        Lx = job->jobVector<double>(job->grid.fp64_props.data());
        Nx = job->jobVector<int>(job->grid.int_props.data());
        hx = Lx.array() / (Nx - job->jobVector(Job::ONES)).array();

        //print grid properties
        std::cout << "Grid properties (Lx = { ";
        for (size_t i=0;i<Lx.rows();i++){
            std::cout << Lx(i) << " ";
        }
        std::cout << "}, Nx = { ";
        for (size_t i=0;i<Nx.rows();i++){
            std::cout << Nx(i) << " ";
        }
        std::cout << "})." << std::endl;
    }

    //call initializing function
    hiddenInit(job);

    std::cout << "Grid Initialized." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void gridWriteFrame(Job* job, Serializer* serializer){
    //nothing to report
    return;
}

/*----------------------------------------------------------------------------*/

std::string gridSaveState(Job* job, Serializer* serializer,std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);

    //create filename
    std::ostringstream s;
    s << "mpm_v2.grid." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    std::string filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    if (ffile.is_open()){
        ffile << "# mpm_v2 grids/cartesian.so\n";
        ffile << Lx.transpose() << "\n" << Nx.transpose() << "\n"; //vectors
        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filepath+filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Grid Saved." << std::endl;

    return filename;
}

/*----------------------------------------------------------------------------*/

int gridLoadState(Job* job, Serializer* serializer, std::string fullpath){
    //job object should be loaded first
    Lx = job->jobVector<double>();
    Nx = job->jobVector<int>();
    hx = job->jobVector<double>();

    std::string line;
    std::stringstream ss;
    std::ifstream fin(fullpath);

    if(fin.is_open()){
        std::getline(fin,line); //first line

        std::getline(fin,line); //Lx
        ss = std::stringstream(line);
        for (size_t i=0;i<Lx.rows();i++){
            ss >> Lx(i);
        }

        std::getline(fin,line); //Nx
        ss = std::stringstream(line);
        for (size_t i=0;i<Nx.rows();i++){
            ss >> Nx(i);
        }

        hx = Lx.array() / (Nx - job->jobVector(Job::ONES)).array();

        //print grid properties
        std::cout << "Grid properties (Lx = { ";
        for (size_t i=0;i<Lx.rows();i++){
            std::cout << Lx(i) << " ";
        }
        std::cout << "}, Nx = { ";
        for (size_t i=0;i<Nx.rows();i++){
            std::cout << Nx(i) << " ";
        }
        std::cout << "})." << std::endl;

        //call hidden initializer
        hiddenInit(job);

        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Grid Loaded." << std::endl;
    return 1;
}

/*----------------------------------------------------------------------------*/

int gridWhichElement(Job* job, Eigen::VectorXd xIN){
    bool inDomain;
    int elementID = floor(xIN[0]/hx[0]); //id in x-dimension
    for (size_t i=0;i<xIN.size();i++){
        inDomain = (xIN[i] <= Lx[i] && xIN[i] >= 0);
        if (!inDomain) { //if xIn is outside domain, return -1
            return -1;
        }
        if (i > 0){
            //add number of elements in next layer of id in higher dimensions
            elementID += floor(xIN[i]/hx[i])*(Nx[i-1]-1);
        }
    }
    return elementID;
}

/*----------------------------------------------------------------------------*/

Eigen::VectorXd gridNodeIDToPosition(Job* job, int idIN){
    Eigen::VectorXi ijk = job->jobVector<int>();
    Eigen::VectorXd tmpVec = job->jobVector<double>();
    int tmp = idIN;
    //find i,j,k representation of node id
    for (size_t i=0;i<ijk.rows();i++){
        ijk(i) = tmp % Nx(i);
        tmp = tmp/Nx(i);
    }
    tmpVec = hx.array() * ijk.array();
    return tmpVec;
}

/*----------------------------------------------------------------------------*/

void gridEvaluateShapeFnValue(Job* job, Eigen::VectorXd xIN, std::vector<int>& nID, std::vector<double>& nVAL){
    Eigen::VectorXd rst = job->jobVector<double>();
    double tmp = 1.0;
    int elementID = gridWhichElement(job,xIN);
    for (size_t n=0;n<nodeIDs.cols();n++){
        //find local coordinates relative to nodal position
        //r = (x_p - x_n)/hx
        rst = (xIN - x_n.row(nodeIDs(elementID,n))).array() / hx.array();
        for (size_t i=0;i<xIN.rows();i++){
            //standard linear hat function
            tmp *= (1 - std::abs(rst(i)));
        }
        nID.push_back(nodeIDs(elementID,n));
        nVAL.push_back(tmp);
        tmp = 1.0;
    }
    return;
}

/*----------------------------------------------------------------------------*/

void gridEvaluateShapeFnGradient(Job* job, Eigen::VectorXd xIN, std::vector<int>& nID, std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& nGRAD){
    Eigen::VectorXd rst = job->jobVector<double>();
    Eigen::VectorXd tmpVec = job->jobVector<double>(Job::ONES);
    double tmp = 1.0;
    int elementID = gridWhichElement(job,xIN);
    for (size_t n=0;n<nodeIDs.cols();n++){
        //find local coordinates relative to nodal position
        //r = (x_p - x_n)/hx
        rst = (xIN - x_n.row(nodeIDs(elementID,n))).array() / hx.array();
        for (size_t i=0;i<xIN.rows();i++){
            //standard linear hat function
            //evaluate at point
            tmp *= (1 - std::abs(rst(i)));
        }
        for (size_t i=0;i<xIN.rows();i++){
            //replace i-direction contribution with sign function
            tmpVec(i) = tmp / (1 - std::abs(rst(i))) * rst(i)/std::abs(rst(i)) / hx(i);
        }
        nID.push_back(nodeIDs(elementID,n));
        nGRAD.push_back(tmpVec);
        tmp = 1.0;
        tmpVec.setZero();
    }
    return;
}


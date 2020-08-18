//
// Created by aaron on 5/15/18.
// cartesian.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "grids.hpp"

/*----------------------------------------------------------------------------*/
//initialize grid assuming that job is set up
void CartesianLinear::init(Job* job){
    if (fp64_props.size() < GRID_DIM || int_props.size() < GRID_DIM){
        std::cout << fp64_props.size() << ", " << int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i dimensions defined.\n",
                __FILE__, __func__, GRID_DIM);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        Lx = KinematicVector(job->JOB_TYPE);
        Nx = Eigen::VectorXi(job->DIM);
        hx = KinematicVector(job->JOB_TYPE);

        for (int pos=0;pos<GRID_DIM;pos++){
            Lx(pos) = fp64_props[pos];
            Nx(pos) = int_props[pos];
            hx(pos) = Lx(pos) / Nx(pos);
        }

        for (int pos=GRID_DIM;pos<hx.size();pos++){
            Lx(pos) = 0;
            Nx(pos) = 0;
            hx(pos) = 0;
        }

        //print grid properties
        std::cout << "Grid properties (Lx = { ";
        for (int i=0;i<Lx.size();i++){
            std::cout << Lx(i) << " ";
        }
        std::cout << "}, Nx = { ";
        for (int i=0;i<Nx.size();i++){
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
//hidden initializer
void CartesianLinear::hiddenInit(Job* job){
    //called from loading or initializing functions
    //needs to have Lx,Nx,hx defined

    //count nodes
    node_count = 1;
    element_count = 1;
    npe = 1;
    for (int i=0;i<GRID_DIM;i++){
        node_count *= (Nx(i)+1);
        element_count *= Nx(i);
        npe *= 2;
    }

    x_n = KinematicVectorArray(node_count,job->JOB_TYPE);
    for (int i=0;i<x_n.size();i++){
        x_n(i) = nodeIDToPosition(job,i);
    }

    //initialize A matrix
    //for mapping position in cube to id
    //0 -> +0,+0,+0
    //1 -> +1,+0,+0
    //...
    //8 -> +1,+1,+1
    Eigen::VectorXi onoff(GRID_DIM);
    A = Eigen::MatrixXi(npe,GRID_DIM);
    onoff.setZero();
    for (int n=0; n<npe;n++){
        for (int i=0;i<onoff.rows();i++){
            A(n,i) = onoff(i);
        }
        for (int i=0;i<onoff.rows();i++) {
            if (onoff(i) == 0){
                onoff(i) = 1;
                break;
            } else {
                onoff(i) = 0;
            }
        }
    }

    //setup element to node map
    Eigen::VectorXi ijk(GRID_DIM);
    int tmp;
    nodeIDs.resize(element_count,npe);
    nodeIDs.setZero();
    for (int e=0;e<nodeIDs.rows();e++){
        tmp = e;
        //find i,j,k count for element position
        for (int i=0;i<ijk.rows();i++){
            ijk(i) = tmp % Nx(i);
            tmp = tmp / Nx(i);
        }

        //find node ids for element
        for (int n=0;n<nodeIDs.cols();n++){
            for (int i=0;i<GRID_DIM;i++) {
                //n = i + j*imax + k*imax*jmax
                //hardcode
                if (i==0) {
                    nodeIDs(e, n) += (ijk(i) + A(n,i));
                } else if (i==1) {
                    nodeIDs(e, n) += (ijk(i) + A(n,i)) * (Nx(i-1)+1);
                } else if (i==2) {
                    nodeIDs(e, n) += (ijk(i) + A(n,i)) * (Nx(i-1)+1) * (Nx(i-2)+1);
                }
            }
        }
    }

    //element volume
    v_e = 1;
    for (int pos=0;pos<GRID_DIM;pos++){
        v_e *= hx(pos);
    }

    //nodal volume
    v_n.resize(x_n.size());
    v_n.setZero();
    for (int e=0;e<nodeIDs.rows();e++){
        for (int c=0;c<nodeIDs.cols();c++) {
            v_n(nodeIDs(e, c)) += v_e / nodeIDs.cols();
        }
    }

    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    //nodal surface integral
    s_n.resize(x_n.size());
    s_n.setZero();
    for (int n=0;n<x_n.size();n++){
        tmp = n;
        for (int i=0;i<ijk.rows();i++){
            ijk(i) = tmp % (Nx(i)+1);
            tmp = tmp/(Nx(i)+1);

            if (ijk(i) == 0 || ijk(i) == Nx(i)){
                s_n(n) = 1;
                for (int pos=0;pos<(GRID_DIM-1);pos++){
                    s_n(n) *= hx(pos);
                }

                //approximate adjustment for axisymetric case
                if (job->JOB_TYPE == job->JOB_AXISYM){
                    tmpVec = nodeIDToPosition(job, n);
                    s_n(n) *= tmpVec(0);
                }
                break;
            }
        }
    }

    return;
}



/*----------------------------------------------------------------------------*/
//
void CartesianLinear::writeFrame(Job* job, Serializer* serializer){
    //nothing to report
    return;
}


/*----------------------------------------------------------------------------*/
//
std::string CartesianLinear::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}


/*----------------------------------------------------------------------------*/
//
int CartesianLinear::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}



/*----------------------------------------------------------------------------*/
//
int CartesianLinear::whichElement(Job* job, KinematicVector& xIN){
    return cartesianWhichElement(job, xIN, Lx, hx, Nx, GRID_DIM);
}

//standard function that other cartesian grids may use
int CartesianLinear::cartesianWhichElement(Job* job, KinematicVector& xIN, KinematicVector& LxIN, KinematicVector& hxIN, Eigen::VectorXi& NxIN, int GRID_DIM_IN){
    assert(xIN.VECTOR_TYPE == LxIN.VECTOR_TYPE && xIN.VECTOR_TYPE == hxIN.VECTOR_TYPE && "cartesianWhichElement failed");
    if (GRID_DIM_IN == -1){
        GRID_DIM_IN = xIN.DIM;
    }

    bool inDomain;
    int elementID = floor(xIN[0] / hxIN[0]); //id in x-dimension
    for (int i = 0; i < GRID_DIM_IN; i++) {
        inDomain = (xIN[i] < LxIN[i] && xIN[i] >= 0);
        if (!inDomain) { //if xIn is outside domain, return -1
            return -1;
        }
        if (i == 1) {
            //add number of elements in next layer of id in higher dimensions
            elementID += floor(xIN[i] / hxIN[i]) * (NxIN[i - 1]);
        } else if (i == 2) {
            elementID += floor(xIN[i] / hxIN[i]) * (NxIN[i - 1]) * (NxIN[i - 2]);
        }
    }
    return elementID;
}


/*----------------------------------------------------------------------------*/
//
bool CartesianLinear::inDomain(Job* job, KinematicVector& xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "inDomain failed");
    for (int i=0;i<GRID_DIM;i++){
        if (!(xIN[i] <= Lx[i] && xIN[i] >= 0)) { //if xIn is outside domain, return -1
            return false;
        }
    }
    return true;
}


/*----------------------------------------------------------------------------*/
//
KinematicVector CartesianLinear::nodeIDToPosition(Job* job, int idIN){
    Eigen::VectorXi ijk(GRID_DIM);
    KinematicVector tmpVec(job->JOB_TYPE);
    int tmp = idIN;
    //find i,j,k representation of node id
    for (int i=0;i<ijk.rows();i++){
        ijk(i) = tmp % (Nx(i)+1);
        tmp = tmp/(Nx(i)+1);
    }
    for (int i=0;i<GRID_DIM;i++){
        tmpVec[i] = hx(i)*ijk(i);
    }
    for (int i=GRID_DIM; i<hx.size(); i++){
        tmpVec[i] = 0;
    }
    return tmpVec;
}



/*----------------------------------------------------------------------------*/
//
void CartesianLinear::evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "evaluateShapeFnValue failed");
    KinematicVector rst(job->JOB_TYPE);
    double tmp = 1.0;
    int elementID = whichElement(job,xIN);
    if (elementID < 0){
        return;
    }
    for (int n=0;n<nodeIDs.cols();n++){
        //find local coordinates relative to nodal position
        for (int i=0;i<GRID_DIM;i++){
            //r = (x_p - x_n)/hx
            rst[i] = (xIN[i] - x_n(nodeIDs(elementID,n),i)) / hx(i);

            //standard linear hat function
            tmp *= (1 - std::abs(rst[i]));
        }
        nID.push_back(nodeIDs(elementID,n));
        nVAL.push_back(tmp);
        tmp = 1.0;
    }
    return;
}


/*----------------------------------------------------------------------------*/
//
void CartesianLinear::evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && nGRAD.VECTOR_TYPE == Lx.VECTOR_TYPE && "evaluateShapeFnGradient failed");
    KinematicVector rst(job->JOB_TYPE);
    KinematicVector tmpVec(job->JOB_TYPE);
    double tmp = 1.0;
    int elementID = whichElement(job,xIN);
    if (elementID < 0){
        return;
    }
    for (int n=0;n<nodeIDs.cols();n++){
        //find local coordinates relative to nodal position
        for (int i=0;i<GRID_DIM;i++){
            //r = (x_p - x_n)/hx
            rst[i] = (xIN[i] - x_n(nodeIDs(elementID,n),i)) / hx(i);

            //standard linear hat function
            //evaluate at point
            tmp *= (1 - std::abs(rst(i)));
        }

        for (int i=0;i<GRID_DIM;i++){
            //replace i-direction contribution with sign function
            //tmpVec(i) = -tmp / (1 - std::abs(rst(i))) * rst(i)/std::abs(rst(i)) / hx(i);
            if (std::abs(rst(i)) > 0){
                tmpVec(i) = -1.0*rst(i)/std::abs(rst(i)) / hx(i);
            } else {
                tmpVec(i) = 0.0;
            }
            for (int j=0; j<GRID_DIM; j++){
                if (j != i){
                    tmpVec(i) *= (1 - std::abs(rst(j)));
                }
            }
        }

        for (int i=GRID_DIM;i<xIN.DIM;i++){
            tmpVec(i) = 0;
        }

        nID.push_back(nodeIDs(elementID,n));
        nGRAD.push_back(tmpVec);
        tmp = 1.0;
        tmpVec.setZero();
    }
    return;
}


/*----------------------------------------------------------------------------*/
//
double CartesianLinear::nodeVolume(Job* job, int idIN){
    return v_n(idIN);
}


/*----------------------------------------------------------------------------*/
//
double CartesianLinear::elementVolume(Job* job, int idIN){
    return v_e;
}


/*----------------------------------------------------------------------------*/
//
int CartesianLinear::nodeTag(Job* job, int idIN){
    return -1;
}

/*----------------------------------------------------------------------------*/
//
double CartesianLinear::nodeSurfaceArea(Job *job, int idIN) {
    return s_n(idIN);
}


/*----------------------------------------------------------------------------*/
//
void CartesianLinear::writeHeader(Job *job, Body *body, Serializer *serializer, std::ofstream &nfile, int SPEC) {
    if (SPEC != Serializer::VTK){
        std::cerr << "ERROR: Unknown file SPEC in writeHeader: " << SPEC  << "! Exiting." << std::endl;
        exit(0);
    }

    int nlen = body->nodes->x.size();

    nfile << "ASCII\n";
    nfile << "DATASET UNSTRUCTURED_GRID\n";

    nfile << "POINTS " << nlen << " double\n";
    for (int i=0;i<nlen;i++){
        //vtk files require x,y,z
        for (int pos = 0; pos < 3; pos++){
            if (pos < body->nodes->x.DIM && (body->nodes->active(i) != 0) && std::isfinite(body->nodes->x(i,pos))){
                nfile << body->nodes->x(i,pos) << " ";
            } else {
                nfile << "0 ";
            }
        }
        nfile << "\n";
    }

    if (GRID_DIM == 1) {
        //use lines
        nfile << "CELLS " << element_count << " " << 3 * element_count << "\n";
        for (int e = 0; e < element_count; e++) {
            nfile << "2 " << nodeIDs(e,0) << " " << nodeIDs(e,1) << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e = 0; e < element_count; e++) {
            nfile << "3\n";
        }
    } else if (GRID_DIM == 2){
        nfile << "CELLS " << element_count << " " << 5 * element_count << "\n";
        for (int e = 0; e < element_count; e++){
            nfile << "4 " << nodeIDs(e,0) << " " << nodeIDs(e,1) << " " << nodeIDs(e,2) << " " << nodeIDs(e,3) << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e=0; e<element_count; e++){
            nfile << "8\n";
        }
    } else if (GRID_DIM == 3){
        nfile << "CELLS " << element_count << " " << 9 * element_count << "\n";
        for (int e = 0; e < element_count; e++){
            nfile << "8 ";
            nfile << nodeIDs(e, 0) << " " << nodeIDs(e, 1) << " ";
            nfile << nodeIDs(e, 2) << " " << nodeIDs(e, 3) << " ";
            nfile << nodeIDs(e, 4) << " " << nodeIDs(e, 5) << " ";
            nfile << nodeIDs(e, 6) << " " << nodeIDs(e, 7) << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e=0; e<element_count; e++){
            nfile << "11\n";
        }
    }

    nfile << "POINT_DATA " << nlen << "\n";
}


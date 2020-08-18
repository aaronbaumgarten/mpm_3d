//
// Created by aaron on 5/25/18.
// cartesian_custom.cpp
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
//
void CartesianCustom::init(Job* job){
    if (fp64_props.size() < job->DIM || int_props.size() < 2*job->DIM){
        std::cout << fp64_props.size() << ", " << int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i dimensions and %i properties defined.\n",
                __FILE__, __func__, job->DIM, job->DIM);
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

        periodic_props = Eigen::VectorXi(GRID_DIM);
        for (int pos=0;pos<job->DIM;pos++){
            periodic_props(pos) = int_props[GRID_DIM + pos];
        }

        //print grid properties
        std::cout << "Grid properties (Lx = { ";
        for (int i=0;i<GRID_DIM;i++){
            std::cout << Lx(i) << " ";
        }
        std::cout << "}, Nx = { ";
        for (int i=0;i<GRID_DIM;i++){
            std::cout << Nx(i) << " ";
        }
        std::cout << "})." << std::endl;
    }

    //call initializing function
    hiddenInit(job);

    std::cout << "Grid Initialized." << std::endl;

    return;
}

void CartesianCustom::hiddenInit(Job* job){
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

    x_n = KinematicVectorArray(node_count, job->JOB_TYPE);
    for (int i=0;i<x_n.size();i++){
        x_n(i) = nodeIDToPosition(job,i);
    }

    //initialize A matrix
    //for mapping position in cube to id
    //0 -> +0,+0,+0
    //1 -> +1,+0,+0
    //...
    //8 -> +1,+1,+1
    Eigen::VectorXi onoff = Eigen::VectorXi(GRID_DIM);
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

    //setup node number to node id map
    Eigen::VectorXi ijk = Eigen::VectorXi(GRID_DIM);
    int tmp;
    nntoni.resize(node_count);
    for (int i=0; i<node_count; i++){
        tmp = i;
        //find i,j,k for node position
        for (int i=0;i<ijk.rows();i++){
            ijk(i) = tmp % (Nx(i)+1);
            tmp = tmp / (Nx(i)+1);
        }

        //wrap x for 2D sim, x and y for 3D sim
        if (GRID_DIM >= 1 && ijk(0) == Nx(0) && periodic_props(0) == PERIODIC){
            ijk(0) = 0;
        }
        if (GRID_DIM >= 2 && ijk(1) == Nx(1) && periodic_props(1) == PERIODIC){
            ijk(1) = 0;
        }
        if (GRID_DIM == 3 && ijk(2) == Nx(2) && periodic_props(2) == PERIODIC){
            ijk(2) = 0;
        }

        tmp = 0;
        //create map
        for (int i=0;i<ijk.rows();i++){
            if (i==0) {
                tmp += ijk(i);
            } else if (i==1) {
                tmp += ijk(i) * (Nx(i-1) + 1);
            } else if (i==2) {
                tmp += ijk(i) * (Nx(i-2) + 1) * (Nx(i-1) + 1);
            }
        }
        nntoni(i) = tmp;
    }

    //setup element to node map
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
            for (int i=0;i<ijk.rows();i++) {
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
            v_n(nntoni(nodeIDs(e, c))) += v_e / nodeIDs.cols();
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
std::string CartesianCustom::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

int CartesianCustom::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//
void CartesianCustom::fixPosition(Job* job, KinematicVector& xIN){
    //wrap position vector around relevant axes

    //int tmp = floor(xIN[1]/hx[1]); //element i,j,k
    //double rem = xIN[1] - tmp*hx[1]; //position relative to element
    //int true_elem = tmp%Nx[1]; //adjusted element
    //double true_pos = true_elem*hx[1] + rem; //adjusted position
    int tmp;
    for (int pos=0; pos<periodic_props.size(); pos++) {
        if (periodic_props(pos) == PERIODIC) {
            tmp = std::floor(xIN[pos] / hx[pos]);
            xIN[pos] += (tmp % Nx[pos] * hx[pos]) - tmp * hx[pos];
            if (xIN[pos] < 0) {
                xIN[pos] += Lx[pos];
            }
        }
    }
    return;
}

int CartesianCustom::whichElement(Job* job, KinematicVector& xIN){
    KinematicVector tmpX = xIN;
    fixPosition(job, tmpX);
    //return standard cartesian which element
    return CartesianLinear::cartesianWhichElement(job, tmpX, Lx, hx, Nx, GRID_DIM);
}

bool CartesianCustom::inDomain(Job* job, KinematicVector& xIN){
    //adjust xIN
    KinematicVector tmpX = xIN;
    fixPosition(job,tmpX);

    for (int i=0;i<GRID_DIM;i++){
        if (!(tmpX[i] <= Lx[i] && tmpX[i] >= 0)) { //if xIn is outside domain, return -1
            return false;
        }
    }
    return true;
}


/*----------------------------------------------------------------------------*/
//
void CartesianCustom::evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL){
    //adjust xIN
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "evaluateShapeFnValue failed");
    KinematicVector tmpVec = xIN;
    fixPosition(job,tmpVec);

    KinematicVector rst = KinematicVector(job->JOB_TYPE);
    double tmp = 1.0;
    int elementID = whichElement(job,tmpVec);
    if (elementID < 0){
        return;
    }
    for (int n=0;n<nodeIDs.cols();n++){
        //find local coordinates relative to nodal position
        //r = (x_p - x_n)/hx
        rst = tmpVec - x_n(nodeIDs(elementID,n));
        for (int i=0;i<GRID_DIM;i++){
            //adjust rst length measures
            rst[i] /= hx[i];
            //standard linear hat function
            tmp *= (1 - std::abs(rst(i)));
        }
        nID.push_back(nntoni(nodeIDs(elementID,n)));
        nVAL.push_back(tmp);
        tmp = 1.0;
    }
    return;
}

void CartesianCustom::evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && nGRAD.VECTOR_TYPE == Lx.VECTOR_TYPE && "evaluateShapeFnGradient failed");
    KinematicVector tmpX = xIN;
    fixPosition(job, tmpX);

    KinematicVector rst(job->JOB_TYPE);
    KinematicVector tmpVec(job->JOB_TYPE);
    double tmp = 1.0;
    int elementID = whichElement(job,tmpX);
    if (elementID < 0){
        return;
    }
    for (int n=0;n<nodeIDs.cols();n++){
        //find local coordinates relative to nodal position
        for (int i=0;i<GRID_DIM;i++){
            //r = (x_p - x_n)/hx
            rst[i] = (tmpX[i] - x_n(nodeIDs(elementID,n),i)) / hx(i);

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

        for (int i=GRID_DIM;i<tmpVec.size();i++){
            tmpVec(i) = 0;
        }
        nID.push_back(nodeIDs(elementID,n));
        nGRAD.push_back(tmpVec);
        tmp = 1.0;
        tmpVec.setZero();
    }
    return;
}
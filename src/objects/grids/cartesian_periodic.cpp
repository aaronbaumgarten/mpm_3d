//
// Created by aaron on 5/25/18.
// cartesian_periodic.cpp
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
void CartesianPeriodic::hiddenInit(Job* job){
    //called from loading or initializing functions
    //needs to have Lx,Nx,hx defined

    //count nodes
    node_count = 1;
    element_count = 1;
    npe = 1;
    for (int i=0;i<Nx.rows();i++){
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
    Eigen::VectorXi onoff = Eigen::VectorXi(job->DIM);
    A = Eigen::MatrixXi(npe,job->DIM);
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
    Eigen::VectorXi ijk = Eigen::VectorXi(job->DIM);
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
        if (job->DIM >= 2 && ijk(0) == Nx(0)){
            ijk(0) = 0;
        }
        if (job->DIM == 3 && ijk(1) == Nx(1)){
            ijk(1) = 0;
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
    for (int pos=0;pos<hx.rows();pos++){
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

    return;
}

/*----------------------------------------------------------------------------*/
//
std::string CartesianPeriodic::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

int CartesianPeriodic::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

/*----------------------------------------------------------------------------*/
//
void CartesianPeriodic::fixPosition(Job* job, KinematicVector& xIN){
    //wrap position vector around relevant axes

    //int tmp = floor(xIN[1]/hx[1]); //element i,j,k
    //double rem = xIN[1] - tmp*hx[1]; //position relative to element
    //int true_elem = tmp%Nx[1]; //adjusted element
    //double true_pos = true_elem*hx[1] + rem; //adjusted position


    int tmp;
    if (job->DIM >= 2){
        tmp = std::floor(xIN[0]/hx[0]);
        xIN[0] += (tmp%Nx[0] * hx[0]) - tmp*hx[0];
        if (xIN[0] < 0){
            xIN[0] += Lx[0];
        }
    }
    if (job->DIM == 3){
        tmp = std::floor(xIN[1]/hx[1]);
        xIN[1] += (tmp%Nx[1] * hx[1]) - tmp*hx[1];
        if (xIN[1] < 0){
            xIN[1] += Lx[1];
        }
    }
    return;
}

int CartesianPeriodic::whichElement(Job* job, KinematicVector& xIN){
    //fix position first, then adjust
    KinematicVector tmp = xIN;
    fixPosition(job, tmp);
    return CartesianLinear::cartesianWhichElement(job, tmp, Lx, hx, Nx);
}

bool CartesianPeriodic::inDomain(Job* job, KinematicVector& xIN){
    //adjust xIN
    KinematicVector tmp = xIN;
    fixPosition(job,tmp);

    for (int i=0;i<tmp.DIM;i++){
        if (!(tmp[i] <= Lx[i] && tmp[i] >= 0)) { //if xIn is outside domain, return -1
            return false;
        }
    }
    return true;
}

/*----------------------------------------------------------------------------*/
//
void CartesianPeriodic::evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL){
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
        for (int i=0;i<tmpVec.rows();i++){
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

void CartesianPeriodic::evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD){
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
        for (int i=0;i<tmpX.DIM;i++){
            //r = (x_p - x_n)/hx
            rst[i] = (tmpX[i] - x_n(nodeIDs(elementID,n),i)) / hx(i);

            //standard linear hat function
            //evaluate at point
            tmp *= (1 - std::abs(rst(i)));
        }
        for (int i=0;i<tmpX.rows();i++){
            //replace i-direction contribution with sign function
            tmpVec(i) = -tmp / (1 - std::abs(rst(i))) * rst(i)/std::abs(rst(i)) / hx(i);
        }
        nID.push_back(nodeIDs(elementID,n));
        nGRAD.push_back(tmpVec);
        tmp = 1.0;
        tmpVec.setZero();
    }
    return;
}
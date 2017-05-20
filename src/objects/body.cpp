//
// Created by aaron on 5/13/17.
// body.cpp
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
#include "points.hpp"

#include "material.hpp"
#include "boundary.hpp"

Body::Body():
        pval(0),
        nval(0),
        phi(0),
        pgrad(0),
        ngrad(0),
        gradphi(0)
{
    id = 0; //when config sets body (set to id)
    name = "default";
    activeMaterial = 0; //when config sets material (set to 1)
    activeBoundary = 0; //when config sets boundary (set to 1)

    points = Points();
    nodes = Nodes();
    material = Material();
    boundary = Boundary();
}

int Body::bodyInit(Job* job){
    points.pointsInit(job,this);
    nodes.nodesInit(job,this);
    material.materialInit(job,this);
    boundary.boundaryInit(job,this);
    return 1;
}

void Body::bodyGenerateMap(Job *job, int use_cpdi = 1) {
    pval.clear();
    nval.clear();
    phi.clear();

    pgrad.clear();
    ngrad.clear();
    gradphi.clear();

    //calculate phi and grad phi
    std::vector<int> nvec(0);
    std::vector<double> valvec(0);
    std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>> gradvec(0);

    for (size_t i=0; i<points.x.rows(); i++){
        job->grid.gridEvaluateShapeFnValue(job,points.x.row(i),nvec,valvec);
        for (size_t j=0; j<nvec.size(); j++){
            pval.push_back((int)i);
            nval.push_back(nvec[j]);
            phi.push_back(valvec[j]);
        }

        job->grid.gridEvaluateShapeFnGradient(job,points.x.row(i),nvec,gradvec);
        for (size_t j=0; j<nvec.size(); j++){
            pgrad.push_back((int)i);
            ngrad.push_back(nvec[j]);
            gradphi.push_back(gradvec[j]);
        }
    }

    return;
}

void Body::bodyCalcNodalValues(Job *job, Eigen::Matrix &nodeVal, Eigen::Matrix &pointVal) {
    //set nodal values based on point values and shape functions
    nodeVal.setZero();
    int nodeID;
    int pointID;

    for (size_t k=0; k<pval.size(); k++){
        pointID = pval[k];
        nodeID = nval[k];
        for (size_t pos=0; pos<pointVal.cols(); pos++){
            nodeVal(nodeID,pos) += phi[k] * pointVal(pointID,pos);
        }
    }

    return;
}

void Body::bodyCalcNodalDivergence(Job *job, Eigen::Matrix &nodeVal, Eigen::Matrix &pointVal) {
    //calculate divergence of a field at nodal positions
    //pass n by 1 nodeVal for p by DIM pointval
    //pass n by DIM nodeVal for p by DIM*DIM pointval
    nodeVal.setZero();
    int nodeID;
    int pointID;

    for (size_t k=0; k<pval.size(); k++){
        pointID = pgrad[k];
        nodeID = ngrad[k];
        for (size_t rpos=0; rpos<nodeVal.cols(); rpos++) {
            for (size_t pos = 0; pos < pointVal.cols(); pos++) {
                //div(u) = dot(grad, u)
                nodeVal(nodeID, rpos) -= gradphi[k](pos) * pointVal(pointID, pos + rpos*(pointVal.cols()));
            }
        }
    }

    return;
}

void Body::bodyCalcPointValues(Job *job, Eigen::Matrix &pointVal, Eigen::Matrix &nodeVal) {
    //set nodal values based on point values and shape functions
    pointVal.setZero();
    int nodeID;
    int pointID;

    for (size_t k=0; k<nval.size(); k++){
        pointID = pval[k];
        nodeID = nval[k];
        for (size_t pos=0; pos<nodeVal.cols(); pos++){
            pointVal(pointID,pos) += phi[k] * nodeVal(nodeID,pos);
        }
    }

    return;
}


void Body::bodyCalcPointGradient(Job *, Eigen::Matrix &pointVal, Eigen::Matrix &nodeVal){
    //calculate gradient of a field at point positions
    //pass n by 1 nodeVal for p by DIM pointval
    //pass n by DIM nodeVal for p by DIM*DIM pointval
    pointVal.setZero();
    int nodeID;
    int pointID;

    for (size_t k=0; k<nval.size(); k++){
        pointID = pgrad[k];
        nodeID = ngrad[k];

        for (size_t rpos=0; rpos<nodeVal.cols(); rpos++) {
            for (size_t pos = 0; pos < pointVal.cols(); pos++) {
                //div(u) = dot(grad, u)
                pointVal(pointID, pos + rpos*(pointVal.cols())) -= gradphi[k](pos) * nodeVal(nodeID, rpos) ;
            }
        }
    }

    return;
}
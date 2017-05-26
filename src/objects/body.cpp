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
        gradphi(0),
        A(0,0)
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
    size_t cp = 1;
    for (size_t i=0;i<job->DIM;i++){
        cp *= 2; //square or cube
    }
    //initialize A matrix
    //for mapping position in cube to id
    //0 -> -1,-1,-1
    //1 -> +1,-1,-1
    //...
    //8 -> +1,+1,+1
    Eigen::VectorXi onoff = -1*job->jobVector<int>(Job::ONES);
    A = job->jobVectorArray(cp);
    for (size_t c=0; c<cp;c++){
        for (size_t i=0;i<onoff.rows();i++){
            A(c,i) = onoff(i);
        }
        for (size_t i=0;i<onoff.rows();i++) {
            if (onoff(i) == -1){
                onoff(i) == 1;
                break;
            } else {
                onoff(i) == -1;
            }
        }
    }

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
    Eigen::VectorXd tmpGrad = job->jobVector<double>();

    for (size_t i=0; i<points.x.rows(); i++){
        if (use_cpdi == 0) {
            nvec.clear();
            valvec.clear();
            job->grid.gridEvaluateShapeFnValue(job,points.x.row(i),nvec,valvec);
            for (size_t j=0; j<nvec.size(); j++){
                pval.push_back((int)i);
                nval.push_back(nvec[j]);
                phi.push_back(valvec[j]);
            }

            nvec.clear();
            gradvec.clear();
            job->grid.gridEvaluateShapeFnGradient(job, points.x.row(i), nvec, gradvec);
            for (size_t j = 0; j < nvec.size(); j++) {
                pgrad.push_back((int) i);
                ngrad.push_back(nvec[j]);
                gradphi.push_back(gradvec[j]);
            }
        } else if (use_cpdi == 1){
            nvec.clear();
            valvec.clear();
            for (size_t c=0;c<A.rows();c++) {
                //spread influence
                job->grid.gridEvaluateShapeFnValue(job, (points.x.row(i) + 0.5*points.extent(i)*A.row(c)), nvec, valvec);
            }
            for (size_t j=0; j<nvec.size(); j++){
                pval.push_back((int)i);
                nval.push_back(nvec[j]);
                phi.push_back(valvec[j]/A.rows()); //reduce influence
            }

            //average gradients along sides of extent
            nvec.clear();
            valvec.clear();
            gradvec.clear();
            for (size_t c=0; c<A.rows(); c++){
                //find shape function value at offset point location
                //add node ids to nodevec
                //add values to valvec
                valvec.clear();
                job->grid.gridEvaluateShapeFnValue(job, (points.x.row(i) + 0.5*points.extent(i)*A.row(c)), nvec, valvec);
                for (size_t v=0; v<valvec.size(); v++){
                    //gradient contribution from corner
                    //G(x) = (S(x+a) - S(x-a))/(2a)
                    tmpGrad = valvec[v]/(A.rows()*0.5*points.extent(i))*A.row(c).transpose();
                    gradvec.push_back(tmpGrad);
                }
            }
            for (size_t j = 0; j < nvec.size(); j++) {
                pgrad.push_back((int) i);
                ngrad.push_back(nvec[j]);
                gradphi.push_back(gradvec[j]);
            }
        } else {
            std::cerr << "Unrecognized Argument for use_cpdi in Body::bodyGenerateMap(): " << use_cpdi << std::endl;
        }
    }

    return;
}

void Body::bodyCalcNodalValues(Job *job, Eigen::Matrix &nodeVal, Eigen::Matrix &pointVal, int SPEC = Body::SET) {
    //set nodal values based on point values and shape functions
    if (SPEC == Body::SET) {
        nodeVal.setZero();
    } else if (SPEC == Body::ADD) {
        // do nothing
    } else {
        std::cerr << "bodyCalcNodalValues(SPEC): Unknown SPEC [" << SPEC << "]" << std::endl;
    }
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

void Body::bodyCalcNodalDivergence(Job *job, Eigen::Matrix &nodeVal, Eigen::Matrix &pointVal, int SPEC = Body::SET) {
    //calculate divergence of a field at nodal positions
    //pass n by 1 nodeVal for p by DIM pointval
    //pass n by DIM nodeVal for p by DIM*DIM pointval
    if (SPEC == Body::SET) {
        nodeVal.setZero();
    } else if (SPEC == Body::ADD) {
        // do nothing
    } else {
        std::cerr << "bodyCalcNodalDivergence(SPEC): Unknown SPEC [" << SPEC << "]" << std::endl;
    }
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

void Body::bodyCalcPointValues(Job *job, Eigen::Matrix &pointVal, Eigen::Matrix &nodeVal, int SPEC = Body::SET) {
    //set nodal values based on point values and shape functions
    if (SPEC == Body::SET) {
        pointVal.setZero();
    } else if (SPEC == Body::ADD) {
        // do nothing
    } else {
        std::cerr << "bodyCalcPointValues(SPEC): Unknown SPEC [" << SPEC << "]" << std::endl;
    }
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


void Body::bodyCalcPointGradient(Job *, Eigen::Matrix &pointVal, Eigen::Matrix &nodeVal, int SPEC = Body::SET){
    //calculate gradient of a field at point positions
    //pass n by 1 nodeVal for p by DIM pointval
    //pass n by DIM nodeVal for p by DIM*DIM pointval
    if (SPEC == Body::SET) {
        pointVal.setZero();
    } else if (SPEC == Body::ADD) {
        // do nothing
    } else {
        std::cerr << "bodyCalcPointValues(SPEC): Unknown SPEC [" << SPEC << "]" << std::endl;
    }
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
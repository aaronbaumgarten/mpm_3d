//
// Created by aaron on 5/5/17.
// body.hpp
//

#ifndef MPM_3D_BODY_HPP
#define MPM_3D_BODY_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

#include "points.hpp"
#include "nodes.hpp"
#include "material.hpp"
#include "boundary.hpp"

class Job;
class Serializer;

class Body{
public:
    //static fields
    static const int SET = 0;
    static const int ADD = 1;

    static const int CPDI_OFF = 0;
    static const int CPDI_ON = 1;

    //body properties here
    int id;
    std::string name;
    int activeMaterial;
    int activeBoundary;

    std::vector<int> pval;
    std::vector<int> nval;
    std::vector<double> phi;

    std::vector<int> pgrad;
    std::vector<int> ngrad;
    std::vector<Eigen::Matrix<double,1,Eigen::Dynamic>, Eigen::aligned_allocator<Eigen::Matrix<double,1,Eigen::Dynamic>>> gradphi;

    Eigen::MatrixXi A; //for mapping cpdi corners to position

    //objects here
    Points points;
    Nodes nodes;
    Material material;
    Boundary boundary;

    //body object specific functions
    Body();
    void bodyInit(Job*);
    std::string bodySaveState(Job*, Serializer*,std::string); //save data to file in serializer directory and return name
    int bodyLoadState(Job*, Serializer*,std::string); //load data from full path

    void bodyGenerateMap(Job*, int SPEC=CPDI_ON); //create shape function mapping vectors (int=1 use cpdi)

    template <typename DerivedA, typename DerivedB> void bodyCalcNodalValues(Job*, Eigen::MatrixBase<DerivedA>& nodeVal, Eigen::MatrixBase<DerivedB>&  pointVal, int SPEC=SET); //calculate nodal value from points
    template <typename DerivedA, typename DerivedB> void bodyCalcNodalGradient(Job*, Eigen::MatrixBase<DerivedA>& nodeVal, Eigen::MatrixBase<DerivedB>&  pointVal, int SPEC=SET); //calculate nodal gradient
    template <typename DerivedA, typename DerivedB> void bodyCalcNodalDivergence(Job*, Eigen::MatrixBase<DerivedA>& nodeVal, Eigen::MatrixBase<DerivedB>&  pointVal, int SPEC=SET); //integrate nodal divergence

    template <typename DerivedA, typename DerivedB> void bodyCalcPointValues(Job*, Eigen::MatrixBase<DerivedA>& pointVal, Eigen::MatrixBase<DerivedB>&  nodeVal, int SPEC=SET); //calculate point values from nodes
    template <typename DerivedA, typename DerivedB> void bodyCalcPointGradient(Job*, Eigen::MatrixBase<DerivedA>& pointVal, Eigen::MatrixBase<DerivedB>&  nodeVal, int SPEC=SET); //calculate point gradients from nodes
    template <typename DerivedA, typename DerivedB> void bodyCalcPointDivergence(Job*, Eigen::MatrixBase<DerivedA>& pointVal, Eigen::MatrixBase<DerivedB>&  nodeVal, int SPEC=SET); //calculate point divergences from nodes

    //other functions
};

template <typename DerivedA, typename DerivedB> void Body::bodyCalcNodalValues(Job *job, Eigen::MatrixBase<DerivedA>& nodeVal, Eigen::MatrixBase<DerivedB>& pointVal, int SPEC /*= Body::SET*/) {
    //set nodal values based on point values and shape functions
    if (nodeVal.cols() > pointVal.cols()){
        std::cerr << "ERROR! Unknown behavior for [" << nodeVal.cols() << "] = interp([" << pointVal.cols() << "]) in bodyCalcNodalValues!" << std::endl;
    }

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

template <typename DerivedA, typename DerivedB> void Body::bodyCalcNodalGradient(Job * job, Eigen::MatrixBase<DerivedA>& nodeVal, Eigen::MatrixBase<DerivedB>& pointVal, int SPEC /*= Body::SET*/){
    //calculate gradient of a field at point positions
    //pass n by 1 nodeVal for p by DIM pointval
    //pass n by DIM nodeVal for p by DIM*DIM pointval
    if (pointVal.cols()*job->DIM < nodeVal.cols()){
        std::cerr << "ERROR! Unknown behavior for [" << nodeVal.cols() << "] = grad([" << pointVal.cols() << "]) with DIM = " << job->DIM << " in bodyCalcNodalGradient!" << std::endl;
    }

    if (SPEC == Body::SET) {
        nodeVal.setZero();
    } else if (SPEC == Body::ADD) {
        // do nothing
    } else {
        std::cerr << "bodyCalcNodalGradient(SPEC): Unknown SPEC [" << SPEC << "]" << std::endl;
    }
    int nodeID;
    int pointID;

    for (size_t k=0; k<nval.size(); k++){
        pointID = pgrad[k];
        nodeID = ngrad[k];

        for (size_t rpos=0; rpos<pointVal.cols(); rpos++) {
            for (size_t pos = 0; pos < nodeVal.cols(); pos++) {
                //div(u) = dot(grad, u)
                nodeVal(nodeID, pos + rpos*(pointVal.cols())) -= gradphi[k](pos) * pointVal(pointID, rpos) ;
            }
        }
    }

    return;
}

template <typename DerivedA, typename DerivedB> void Body::bodyCalcNodalDivergence(Job *job, Eigen::MatrixBase<DerivedA>& nodeVal, Eigen::MatrixBase<DerivedB>& pointVal, int SPEC /*= Body::SET*/) {
    //calculate divergence of a field at nodal positions
    //pass n by 1 nodeVal for p by DIM pointval
    //pass n by DIM nodeVal for p by DIM*DIM pointval
    if (pointVal.cols() < job->DIM*nodeVal.cols()){
        std::cerr << "ERROR! Unknown behavior for [" << nodeVal.cols() << "] = div([" << pointVal.cols() << "]) with DIM = " << job->DIM << " in bodyCalcNodalDivergence!" << std::endl;
    }
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
            for (size_t pos = 0; pos < job->DIM; pos++){    //pointVal.cols(); pos++) {
                //div(u) = dot(grad, u)
                nodeVal(nodeID, rpos) -= gradphi[k](pos) * points.v(pointID) * pointVal(pointID, pos + rpos*(job->DIM));    //pointVal.cols()));
            }
        }
    }

    return;
}

template <typename DerivedA, typename DerivedB> void Body::bodyCalcPointValues(Job *job, Eigen::MatrixBase<DerivedA>& pointVal, Eigen::MatrixBase<DerivedB>& nodeVal, int SPEC /*= Body::SET*/) {
    //set nodal values based on point values and shape functions
    if (nodeVal.cols() < pointVal.cols()){
        std::cerr << "ERROR! Unknown behavior for [" << pointVal.cols() << "] = interp([" << nodeVal.cols() << "]) in bodyCalcPointValues!" << std::endl;
    }

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


template <typename DerivedA, typename DerivedB> void Body::bodyCalcPointGradient(Job * job, Eigen::MatrixBase<DerivedA>& pointVal, Eigen::MatrixBase<DerivedB>& nodeVal, int SPEC /*= Body::SET*/){
    //calculate gradient of a field at point positions
    //pass n by 1 nodeVal for p by DIM pointval
    //pass n by DIM nodeVal for p by DIM*DIM pointval
    if (nodeVal.cols()*job->DIM < pointVal.cols()){
        std::cerr << "ERROR! Unknown behavior for [" << pointVal.cols() << "] = grad([" << nodeVal.cols() << "]) with DIM = " << job->DIM << " in bodyCalcPointGradient!" << std::endl;
    }

    if (SPEC == Body::SET) {
        pointVal.setZero();
    } else if (SPEC == Body::ADD) {
        // do nothing
    } else {
        std::cerr << "bodyCalcPointGradient(SPEC): Unknown SPEC [" << SPEC << "]" << std::endl;
    }
    int nodeID;
    int pointID;

    for (size_t k=0; k<nval.size(); k++){
        pointID = pgrad[k];
        nodeID = ngrad[k];

        for (size_t rpos=0; rpos<nodeVal.cols(); rpos++) {
            for (size_t pos = 0; pos < job->DIM; pos++){    //pointVal.cols(); pos++) {
                pointVal(pointID, pos + rpos*(job->DIM)) += gradphi[k](pos) * nodeVal(nodeID, rpos) ;
            }
        }
    }

    return;
}

template <typename DerivedA, typename DerivedB> void Body::bodyCalcPointDivergence(Job* job, Eigen::MatrixBase<DerivedA>& pointVal, Eigen::MatrixBase<DerivedB>&  nodeVal, int SPEC /*= Body::SET*/){
    //calculate divergence of a field at nodal positions
    //pass n by 1 pointVal for p by DIM nodeVal
    //pass n by DIM pointVal for p by DIM*DIM nodeVal
    if (nodeVal.cols() < job->DIM*pointVal.cols()){
        std::cerr << "ERROR! Unknown behavior for [" << pointVal.cols() << "] = div([" << nodeVal.cols() << "]) with DIM = " << job->DIM << " in bodyCalcPointDivergence!" << std::endl;
    }

    if (SPEC == Body::SET) {
        pointVal.setZero();
    } else if (SPEC == Body::ADD) {
        // do nothing
    } else {
        std::cerr << "bodyCalcPointDivergence(SPEC): Unknown SPEC [" << SPEC << "]" << std::endl;
    }
    int nodeID;
    int pointID;

    for (size_t k=0; k<nval.size(); k++){
        pointID = pgrad[k];
        nodeID = ngrad[k];
        for (size_t rpos=0; rpos<pointVal.cols(); rpos++) {
            for (size_t pos = 0; pos < job->DIM; pos++){    //pointVal.cols(); pos++) {
                //div(u) = dot(grad, u)
                pointVal(pointID, rpos) += gradphi[k](pos) * nodeVal(nodeID, pos + rpos*(job->DIM));
            }
        }
    }
};

#endif //MPM_3D_BODY_HPP

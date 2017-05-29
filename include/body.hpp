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

    //other functions
};

#endif //MPM_3D_BODY_HPP

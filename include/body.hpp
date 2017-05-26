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

class Job;
class Serializer;
class Points;
class Nodes;
class Material;
class Boundary;

class Body{
public:
    //static fields
    static int SET = 0;
    static int ADD = 1;

    static int CPDI_OFF = 0;
    static int CPDI_ON = 1;

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
    int bodyInit(Job*);

    void bodyGenerateMap(Job*, int=CPDI_ON); //create shape function mapping vectors (int=1 use cpdi)

    void bodyCalcNodalValues(Job*, Eigen::Matrix& nodeVal, Eigen::Matrix& pointVal, int=SET); //calculate nodal value from points
    void bodyCalcNodalGradient(Job*, Eigen::Matrix& nodeVal, Eigen::Matrix& pointVal, int=SET); //calculate nodal gradient
    void bodyCalcNodalDivergence(Job*, Eigen::Matrix& nodeVal, Eigen::Matrix& pointVal, int=SET); //integrate nodal divergence
    void bodyCalcPointValues(Job*, Eigen::Matrix& pointVal, Eigen::Matrix& nodeVal, int=SET); //calculate point values from nodes
    void bodyCalcPointGradient(Job*, Eigen::Matrix& pointVal, Eigen::Matrix& nodeVal, int=SET); //calculate point gradients from nodes

    //other functions
};

#endif //MPM_3D_BODY_HPP

//
// Created by aaron on 8/26/16.
// element.hpp
//

#ifndef MPM_3D_ELEMENT_HPP
#define MPM_3D_ELEMENT_HPP

#include <stdlib.h>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <cmath>
#include "body.hpp"
#include "node.hpp"
#include "particle.hpp"

#define S_L(x) (1.0-(x))//right hand shape function
#define S_R(x) (x)//left hand shape function

class Elements{
public:
    //unique id
    //size_t id;

    //node info
    size_t numElements;
    size_t numNodesPerElement;
    Eigen::MatrixXi nodeID;

    //particle info
    //size_t numParticles;
    //std::vector<size_t> particleID;

    //element color (for threading)
    //int color;

    Elements() {
        numElements = 0;
        numNodesPerElement = 0;
        nodeID.resize(0,0);
    }
    Elements(size_t,size_t);

    void resizeElements(size_t,size_t);
    void addElement(Eigen::VectorXi nodeIDs, size_t idIn);
    void calculatePhic(Body *body, size_t ide, size_t idp, size_t idc, int use_cpdi);
};

#endif //MPM_3D_ELEMENT_HPP

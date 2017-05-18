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
    //body properties here
    int activeMaterial;
    int activeBoundary;

    //objects here
    Points points;
    Nodes nodes;
    Material material;
    Boundary boundary;

    //job object specific functions
    Body();
    int bodyInit(Job*);

    //other functions
};

#endif //MPM_3D_BODY_HPP

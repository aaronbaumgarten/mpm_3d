//
// Created by aaron on 5/15/18.
// boundaries.hpp
//

#ifndef MPM_V3_BOUNDARIES_HPP
#define MPM_V3_BOUNDARIES_HPP

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

/*
 * IN THIS FILE, DEFINE NODES OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
 * class Boundary : public MPMObject{
public:
    virtual void init(Job*, Body*) = 0;             //initialize from Job and Body
    virtual void generateRules(Job*, Body*) = 0;    //generate boundary rules
    virtual void applyRules(Job*, Body*) = 0;       //apply boundary rules

    virtual void writeFrame(Job*, Body*, Serializer*) = 0;                      //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;   //save state to file
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;           //load from fullpath
};
 */

class CartesianBox : public Boundary{
public:
    CartesianBox(){
        object_name = "CartesianBox";
    }

    KinematicVectorArray bcNodalMask; //store dof to control

    void init(Job* job, Body* body);
    void generateRules(Job* job, Body* body);
    void applyRules(Job* job, Body* body);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

#endif //MPM_V3_BOUNDARIES_HPP

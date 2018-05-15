//
// Created by aaron on 5/15/18.
// nodes.hpp
//

#ifndef MPM_V3_NODES_HPP
#define MPM_V3_NODES_HPP

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>
#include <sstream>
#include <Eigen/Core>
#include <ctime>

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
 * class Nodes : public MPMObject{
public:
    KinematicVectorArray x, u, x_t, a, mx_t, f; //state vectors
    Eigen::VectorXd m, V;                       //weight measures
    Eigen::VectorXi active;                     //active?

    virtual void init(Job*, Body*) = 0;                                         //initialize from Job and Body
    virtual void writeFrame(Job*, Body*, Serializer*) = 0;                      //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;   //save to file (in given directory)
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;           //load from full path
};
 */

class DefaultNodes : public Nodes{
public:
    DefaultNodes(){
        object_name = "DefaultNodes";
    }

    KinematicVectorArray diff_x_t; //dummy vector array for standard mpm

    void init(Job* job, Body* body);
    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

#endif //MPM_V3_NODES_HPP

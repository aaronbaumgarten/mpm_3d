//
// Created by aaron on 5/16/18.
// bodies.hpp
//

#ifndef MPM_V3_BODIES_HPP
#define MPM_V3_BODIES_HPP

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
 * IN THIS FILE, DEFINE BODY OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
 * class Body : public MPMObject{
public:
    int id;             //numerical id of body in simulation
    std::string name;   //name of body in simulation
    int activeMaterial; //is the continuum body material defined?
    int activeBoundary; //is the continuum body boundary defined?

    MPMScalarSparseMatrix S; //S_ip maps ith node to pth point
    KinematicVectorSparseMatrix gradS; //gradS_ip maps ith node gradient to pth point

    std::unique_ptr<Points> points;
    std::unique_ptr<Nodes> nodes;
    std::unique_ptr<Material> material;
    std::unique_ptr<Boundary> boundary;

    virtual void init(Job*) = 0;                                        //initialiaze from Job
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save to file (in given directory)
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load data from full path
    virtual void generateMap(Job*, int) = 0;                            //generate S and gradS
};
 */


class DefaultBody : public Body{
public:
    DefaultBody(){
        object_name = "DefaultBody";
    }

    //not true CPDI, see Dunatunga and Kamrin 2017
    static const int CPDI_ON = 1;
    static const int CPDI_OFF = 0;

    Eigen::MatrixXi A; //for mapping corners to position

    void init(Job* job);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);
    void generateMap(Job* job, int SPEC);
};

#endif //MPM_V3_BODIES_HPP

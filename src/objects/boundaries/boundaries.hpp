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
 * IN THIS FILE, DEFINE BOUNDARY OBJECTS.
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

/*----------------------------------------------------------------------------*/

class CartesianSmoothBox : public Boundary{
public:
    CartesianSmoothBox(){
        object_name = "CartesianSmoothBox";
    }

    KinematicVectorArray bcNodalMask; //store dof control

    void init(Job* job, Body* body);
    void generateRules(Job* job, Body* body);
    void applyRules(Job* job, Body* body);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

/*----------------------------------------------------------------------------*/

class CartesianFrictionalBox : public Boundary{
public:
    CartesianFrictionalBox(){
        object_name = "CartesianFrictionalBox";
        mu_f = 0;
    }

    double mu_f;
    KinematicVectorArray bcNodalMask; //store dof control
    KinematicVectorArray bcNodalForce; //store forces
    KinematicVectorArray pvec;
    MaterialTensorArray tmp;
    MaterialVectorArray nvec;

    void init(Job* job, Body* body);
    void generateRules(Job* job, Body* body);
    void applyRules(Job* job, Body* body);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

/*----------------------------------------------------------------------------*/

class CartesianBoxCustom : public Boundary{
public:
    CartesianBoxCustom(){
        object_name = "CartesianBoxCustom";
        mu_f = 0;
    }

    static const int NO_SLIP_WALL       = 0;
    static const int FRICTION_LESS_WALL = 1;
    static const int FRICTIONAL_WALL    = 2;
    static const int PERIODIC           = 3;
    static const int DRIVEN_VELOCITY    = 4;

    double mu_f;
    double v_set;
    KinematicVector Lx;
    Eigen::VectorXi limit_props;
    KinematicVectorArray bcNodalMask, bcNodalForce;
    MaterialTensorArray tmp;
    MaterialVectorArray nvec;
    KinematicVectorArray pvec;

    void init(Job* job, Body* body);
    void generateRules(Job* job, Body* body);
    void applyRules(Job* job, Body* body);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

/*----------------------------------------------------------------------------*/

class Regular2DTaylorCouetteCustom : public Boundary{
public:
    Regular2DTaylorCouetteCustom(){
        object_name = "Regular2DTaylorCouetteCustom";
        inner_fp64_prop = 0;
        outer_fp64_prop = 0;
    }

    static const int NO_SLIP_WALL       = 0;
    static const int FRICTION_LESS_WALL = 1;
    static const int FRICTIONAL_WALL    = 2;
    static const int DRIVEN_VELOCITY    = 3;

    int inner_int_prop, outer_int_prop;
    double inner_fp64_prop, outer_fp64_prop;
    Eigen::VectorXi bcNodalMask;
    KinematicVectorArray bcNodalForce, bcReactionForce;
    MaterialTensorArray tmp;
    MaterialVectorArray nvec;
    KinematicVectorArray pvec;

    void init(Job* job, Body* body);
    void generateRules(Job* job, Body* body);
    void applyRules(Job* job, Body* body);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

#endif //MPM_V3_BOUNDARIES_HPP

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

    bool generate_friction = false;
    bool generate_tractions = false;
    bool generate_periodic = false;

    static const int FREE_BOUNDARY      = -1;
    static const int NO_SLIP_WALL       = 0;
    static const int FRICTION_LESS_WALL = 1;
    static const int FRICTIONAL_WALL    = 2;
    static const int PERIODIC           = 3;
    static const int DRIVEN_VELOCITY    = 4;
    static const int DRIVEN_TRACTION    = 5;
    static const int DRIVEN_VELOCITY_BOUNDED_TRACTION = 6;

    double mu_f, tau_max;
    KinematicVectorArray v_set;
    KinematicVector Lx_min, Lx_max; //just for axisym...
    Eigen::VectorXi limit_props;
    KinematicVectorArray bcNodalMask, bcNodalForce;
    MaterialTensorArray tmp;
    MaterialVectorArray nvec;
    KinematicVectorArray pvec;
    Eigen::VectorXd v_n;

    void init(Job* job, Body* body);
    void generateRules(Job* job, Body* body);
    void applyRules(Job* job, Body* body);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

/*----------------------------------------------------------------------------*/

class AssignedVelocity : public Boundary{
public:
    AssignedVelocity(){
        object_name = "AssignedVelocity";
    }

    KinematicVector v_set;

    void init(Job* job, Body* body);
    void generateRules(Job* job, Body* body);
    void applyRules(Job* job, Body* body);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

/*----------------------------------------------------------------------------*/

class GeneralCustomBoundary : public Boundary{
public:
    GeneralCustomBoundary(){
        object_name = "GeneralCustomBoundary";
    }

    bool USE_Y_AXIS_SYMMETRY = false;
    static const int Y_AXIS_SLIP        = 255;

    static const int FREE_BOUNDARY      = -1;
    static const int DIRICHLET          = 0;
    static const int SLIP_BOUNDARY      = 1;
    static const int MAX_PROP_VALUE     = 1;

    Eigen::VectorXi limit_props, bcTag;
    std::vector<KinematicVector> limit_vals;
    KinematicVectorArray bcNodalForce;
    KinematicVectorArray bcValues;
    MaterialTensorArray tmp;
    MaterialVectorArray nvec;
    KinematicVectorArray pvec;
    Eigen::VectorXd v_n;

    void init(Job* job, Body* body);
    void generateRules(Job* job, Body* body);
    void applyRules(Job* job, Body* body);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

#endif //MPM_V3_BOUNDARIES_HPP

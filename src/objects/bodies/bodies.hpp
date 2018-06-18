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

    virtual void generateLoads(Job*) = 0;
    virtual void applyLoads(Job*) = 0;
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

    virtual void init(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);
    virtual void generateMap(Job* job, int SPEC);

    virtual void generateLoads(Job* job);
    virtual void applyLoads(Job* job);
};

/*----------------------------------------------------------------------------*/

class WheelBody : public DefaultBody{
public:
    WheelBody(){
        object_name = "WheelBody";
    }

    double omega, t_start, m_cp;
    MaterialVector a, b, c;
    KinematicVector r, mv_cp, v_cp, mx_cp, x_cp, dv_cp, v_n, v_t;

    void init(Job* job);

    void generateLoads(Job* job);
    void applyLoads(Job* job);
};

/*----------------------------------------------------------------------------*/

class HydrostaticBody : public DefaultBody{
public:
    HydrostaticBody(){
        object_name = "HydrostaticBody";
    }

    double g, effective_density;
    void init(Job* job);
};

#endif //MPM_V3_BODIES_HPP

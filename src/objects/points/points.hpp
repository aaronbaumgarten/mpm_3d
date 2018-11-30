//
// Created by aaron on 5/14/18.
// points.hpp
//

#ifndef MPM_V3_POINTS_HPP
#define MPM_V3_POINTS_HPP

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
 * IN THIS FILE, DEFINE POINTS OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
 * class Points : public MPMObject{
public:
    std::string file;

    KinematicVectorArray x, u, x_t, mx_t, b;    //state vectors
    KinematicTensorArray L;                     //velocity gradient
    MaterialTensorArray T;                      //cauchy stress
    Eigen::VectorXd m, v, v0, extent;           //weight measures
    Eigen::VectorXi active;                     //active?

    virtual void init(Job*, Body*) = 0;                         //initialize from Job and Body
    virtual void readFromFile(Job*, Body*, std::string) = 0;    //construct points from given file
    virtual void generateMap(Job*, Body*, int) = 0;             //generate S and gradS
    virtual void updateIntegrators(Job*, Body*) = 0;            //update integrators (extent, etc.)

    virtual void writeHeader(Job*, Body*, Serializer*, std::ofstream&, int) = 0;
    virtual void writeFrame(Job*, Body*, Serializer*) = 0;                      //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;   //save to file (in given directory)
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;           //load data from full path

    virtual void generateLoads(Job*, Body*) = 0;   //arbitrary loading during simulation
    virtual void applyLoads(Job*, Body*) = 0;
};
 */

class DefaultPoints : public Points{
public:
    DefaultPoints(){
        object_name = "DefaultPoints";
    }

    Eigen::VectorXd extent;
    Eigen::MatrixXi A; //for mapping corners to position

    virtual void init(Job* job, Body* body);
    virtual void readFromFile(Job* job, Body* body, std::string fileIN);
    virtual void generateMap(Job* job, Body* body, int SPEC);             //generate S and gradS
    virtual void updateIntegrators(Job* job, Body* body);                 //update integrators (extent, etc.)

    virtual void writeHeader(Job* job, Body* body, Serializer* serializer, std::ofstream& pfile, int SPEC);
    virtual void writeFrame(Job* job, Body* body, Serializer* serializer);
    virtual std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);

    virtual void generateLoads(Job* job, Body* body);   //arbitrary loading during simulation
    virtual void applyLoads(Job* job, Body* body);
};

/*----------------------------------------------------------------------------*/

//general object for creating sets of points
class Part{
public:
    virtual bool encompasses(KinematicVector&) = 0;
};

class Ball : public Part{
public:
    Ball(KinematicVector oIN, double rIN){
        r = rIN;
        o = oIN;
    }

    double r;           //radius
    KinematicVector o;  //origin

    bool encompasses(KinematicVector& xIN);
};

class Box : public Part{
public:
    Box(KinematicVector x_minIN, double x_maxIN){
        x_min = x_minIN;
        x_max = x_maxIN;
    }

    KinematicVector x_min, x_max;  //bounds

    bool encompasses(KinematicVector& xIN);
};

/*----------------------------------------------------------------------------*/

class GmshPoints : public DefaultPoints{
public:
    GmshPoints(){
        object_name = "GmshPoints";
    }

    //linear material points per cell and density of points
    double lmpp, rho;

    void init(Job* job, Body* body);
    void readFromFile(Job* job, Body* body, std::string fileIN);
};

#endif //MPM_V3_POINTS_HPP

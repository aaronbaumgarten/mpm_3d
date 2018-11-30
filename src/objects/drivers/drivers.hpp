//
// Created by aaron on 5/10/18.
// drivers.hpp
//

#ifndef MPM_V3_DRIVERS_HPP
#define MPM_V3_DRIVERS_HPP

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "job.hpp"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <time.h>

/*
 * IN THIS FILE, DEFINE DRIVER OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
class Driver : public MPMObject{
public:
    //functions which must be implemented by every driver
    virtual void init(Job*) = 0;                                        //initialize from Job
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save to file (in given directory) and return filename
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load from file
    virtual void run(Job*) = 0;                                         //run mpm according to problem
    virtual void generateGravity(Job*) = 0;                             //generate gravity
    virtual void applyGravity(Job*) = 0;                                //apply gravity
};
 */

class DefaultDriver : public Driver{
public:
    DefaultDriver(){
        object_name = "DefaultDriver"; //set object name here
    }

    double stop_time;
    KinematicVector gravity;

    virtual void init(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    virtual void run(Job* job);
    virtual void generateGravity(Job* job);
    virtual void applyGravity(Job* job);
};

/*----------------------------------------------------------------------------*/

class ColumnCollapseDriver : public Driver{
public:
    ColumnCollapseDriver(){
        object_name = "ColumnCollapseDriver";
    }

    double stop_time;
    KinematicVector gravity;

    void setPressure(Job* job);

    void init(Job* job);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void run(Job* job);
    void generateGravity(Job* job);
    void applyGravity(Job* job);
};

/*----------------------------------------------------------------------------*/

class UserDefinedGravityDriver : public Driver{
public:
    UserDefinedGravityDriver(){
        object_name = "UserDefinedGravityDriver";
    }

    double stop_time;
    KinematicVector gravity;

    void init(Job* job);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void run(Job* job);
    void generateGravity(Job* job);
    void applyGravity(Job* job);
};

/*----------------------------------------------------------------------------*/

class CavityFlowDriver : public Driver{
public:
    CavityFlowDriver(){
        object_name = "CavityFlowDriver";
    }

    double stop_time;
    double Ly, hy;
    double v_set;

    void init(Job* job);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void run(Job* job);
    void generateGravity(Job* job);
    void applyGravity(Job* job);
};


/*----------------------------------------------------------------------------*/
//driver for ballistic impact problems
//will launch ballistic object at specified velocity at t=0

class BallisticDriver : public DefaultDriver{
public:
    BallisticDriver(){
        object_name = "BallisticDriver"; //set object name here
    }

    int ballistic_id;
    double stop_time;
    KinematicVector gravity, velocity;

    void init(Job* job);

    void run(Job* job);
};

#endif //MPM_V3_DRIVERS_HPP

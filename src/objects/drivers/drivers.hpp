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

    void init(Job* job);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void run(Job* job);
    void generateGravity(Job* job);
    void applyGravity(Job* job);
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

    void setPressure(Job* job);

    void init(Job* job);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void run(Job* job);
    void generateGravity(Job* job);
    void applyGravity(Job* job);
};

#endif //MPM_V3_DRIVERS_HPP

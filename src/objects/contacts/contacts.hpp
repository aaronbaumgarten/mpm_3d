//
// Created by aaron on 5/15/18.
// contacts.hpp
//

#ifndef MPM_V3_CONTACTS_HPP
#define MPM_V3_CONTACTS_HPP

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
 * IN THIS FILE, DEFINE CONTACT OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
 * class Contact : public MPMObject{
public:
    static const int IMPLICIT = 1; //explicit contact rules
    static const int EXPLICIT = 0; //implicit contact rules

    int id;             //numerical id of contact in simulation
    std::string name;   //name of contact in simulation

    virtual void init(Job*) = 0;            //initialize from Job
    virtual void generateRules(Job*) = 0;   //generate contact rules
    virtual void applyRules(Job*, int) = 0; //apply contact rules

    virtual void writeFrame(Job*, Serializer*) = 0;                     //send frame data to Serializer
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save state to file
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load from file
};
 */

class ContactHuang : public Contact{
public:
    ContactHuang(){
        object_name = "ContactHuang";
    }

    double mu_f = 0.4;
    std::vector<int> bodyIDs = {-1,-1};
    KinematicVectorArray contact_normal, contact_force;

    void init(Job* job);
    void generateRules(Job* job);
    void applyRules(Job* job, int SPEC);

    void writeFrame(Job* job, Serializer* serializer);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);
};

/*----------------------------------------------------------------------------*/

class SlurryMixture : public Contact{
public:
    SlurryMixture(){
        object_name = "SlurryMixture";
    }

    double grains_rho = 2500; //kg/m^3
    double grains_d = 0.001; //m
    double eta_0 = 8.9e-4; //Pa*s
    double fluid_rho = 1000; //kg/m^3

    std::vector<int> bodyIDs = {-1,-1};
    int solid_body_id = -1;
    int fluid_body_id = -1;
    int spec_override = -1;

    Eigen::VectorXd n, V;
    KinematicVectorArray divT;

    void init(Job* job);
    void generateRules(Job* job);
    void applyRules(Job* job, int SPEC);

    void writeFrame(Job* job, Serializer* serializer);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);
};

/*----------------------------------------------------------------------------*/

class SlurryContact : public Contact{
public:
    SlurryContact(){
        object_name = "SlurryContact";
    }

    double mu_f = 0.4;
    std::vector<int> bodyIDs = {-1,-1,-1};
    KinematicVectorArray contact_normal, contact_force;

    virtual void init(Job* job);
    virtual void generateRules(Job* job);
    virtual void applyRules(Job* job, int SPEC);

    virtual void writeFrame(Job* job, Serializer* serializer);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);
};

/*----------------------------------------------------------------------------*/

class SlurryContact_ReflectedBoundary : public SlurryContact{
public:
    SlurryContact_ReflectedBoundary(){
        object_name = "SlurryContact_ReflectedNormal";
    }

    KinematicVector Lx;

    void init(Job* job);
    void generateRules(Job* job);
};

/*----------------------------------------------------------------------------*/

class ContactHuang_ReflectedBoundary : public ContactHuang{
public:
    ContactHuang_ReflectedBoundary(){
        object_name = "ContactHuang_ReflectedNormal";
    }

    KinematicVector Lx;

    void init(Job* job);
    void generateRules(Job* job);
};


#endif //MPM_V3_CONTACTS_HPP

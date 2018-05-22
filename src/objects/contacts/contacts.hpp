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
    KinematicVectorArray contact_normal, mv1_k, mv2_k;

    void init(Job* job);
    void generateRules(Job* job);
    void applyRules(Job* job, int SPEC);

    void writeFrame(Job* job, Serializer* serializer);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);
};

#endif //MPM_V3_CONTACTS_HPP

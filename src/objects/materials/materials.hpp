//
// Created by aaron on 5/15/18.
// materials.hpp
//

#ifndef MPM_V3_MATERIALS_HPP
#define MPM_V3_MATERIALS_HPP

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
 * IN THIS FILE, DEFINE MATERIAL OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
 * class Material : public MPMObject{
public:
    static const int UPDATE = 1;    //update internal variables
    static const int TRIAL = 0;     //don't update internal variables

    virtual void init(Job*, Body*) = 0;                                     //initialize from Job and Body
    virtual void calculateStress(Job*, Body*, int) = 0;                     //calculate stress update
    virtual void assignStress(Job*, Body*, MaterialTensor&, int, int) = 0;  //assign stress to given id
    virtual void assignPressure(Job*, Body*, double, int, int) = 0;         //assign pressure to given id

    virtual void writeFrame(Job*, Body*, Serializer*) = 0;          //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std:;string) = 0;    //save state to file
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;            //load from full path
};
 */

class IsotropicLinearElasticity : public Material {
public:
    IsotropicLinearElasticity(){
        object_name = "IsotropicLinearElasticity";
    }

    //material properties
    double E, nu, G, K, lambda;

    void init(Job* job, Body* body);
    void calculateStress(Job* job, Body* body, int SPEC);
    void assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC);
    void assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};


/*----------------------------------------------------------------------------*/
//
class Sand_SachithLocal : public Material {
public:
    Sand_SachithLocal(){
        object_name = "Sand_SachithLocal";
    }

    //material properties
    double E, nu, G, K, lambda;
    double grains_d;
    Eigen::VectorXd gammap, gammap_dot;

    double MU_S = 0.280;
    double GRAINS_RHO = 2500;
    double RHO_CRITICAL = (GRAINS_RHO*0.6);
    double MU_2 = MU_S;
    double DELTA_MU = (MU_2 - MU_S);
    double I_0 = 0.278;

    void quadratic_roots(double *x1, double *x2, double a, double b, double c);
    double negative_root(double a, double b, double c);

    void init(Job* job, Body* body);
    void calculateStress(Job* job, Body* body, int SPEC);
    void assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC);
    void assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

#endif //MPM_V3_MATERIALS_HPP

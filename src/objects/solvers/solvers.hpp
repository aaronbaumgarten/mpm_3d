//
// Created by aaron on 5/12/18.
// solvers.hpp
//

#ifndef MPM_V3_SOLVERS_HPP
#define MPM_V3_SOLVERS_HPP

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

/*
 * IN THIS FILE, DEFINE SOLVER OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
class Solver : public MPMObject{
public:
    //functions that must be implemented by every solver (not many)
    virtual void init(Job*) = 0;                                        //initialize from Job
    virtual void step(Job*) = 0;                                        //perform single mpm step
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save to file (in given directory) and return filename
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load from file
};
 */

class ExplicitUSL : public Solver{
public:
    ExplicitUSL(){
        object_name = "ExplicitUSL"; //set name of object from registry
    }

    int cpdi_spec = 1;
    int contact_spec = Contact::IMPLICIT;

    void init(Job* job);
    void step(Job* job);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void createMappings(Job *job);
    void mapPointsToNodes(Job* job);
    void generateContacts(Job* job);
    void addContacts(Job* job);
    void generateBoundaryConditions(Job* job);
    void addBoundaryConditions(Job* job);
    void moveGrid(Job* job);
    void movePoints(Job* job);
    void calculateStrainRate(Job* job);
    void updateDensity(Job* job);
    void updateStress(Job* job);
    void generateLoads(Job* job);
    void applyLoads(Job* job);
};

#endif //MPM_V3_SOLVERS_HPP

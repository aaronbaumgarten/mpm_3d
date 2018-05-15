//
// Created by aaron on 5/12/18.
// solvers.hpp
//

#ifndef MPM_V3_SOLVERS_HPP
#define MPM_V3_SOLVERS_HPP

#ifndef MPM_V3_DRIVERS_HPP
#define MPM_V3_DRIVERS_HPP

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
 * IN THIS FILE, DEFINE DRIVER OBJECTS.
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

    void init(Job* job);
    void step(Job* job);

    void createMappings(Job *job);
    void mapPointsToNodes(Job* job);
    void generateContacts(Job* job);
    void addContacts(Job* job);
    void generateBoundaryConditions(Job* job);
    void addBoundaryCondition(Job* job);
    //??
};

#endif //MPM_V3_SOLVERS_HPP

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
#include <thread>

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

    virtual void init(Job* job);
    virtual void step(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    virtual void createMappings(Job *job);
    virtual void mapPointsToNodes(Job* job);
    virtual void generateContacts(Job* job);
    virtual void addContacts(Job* job);
    virtual void generateBoundaryConditions(Job* job);
    virtual void addBoundaryConditions(Job* job);
    virtual void moveGrid(Job* job);
    virtual void movePoints(Job* job);
    virtual void calculateStrainRate(Job* job);
    virtual void updateDensity(Job* job);
    virtual void updateStress(Job* job);
    virtual void generateLoads(Job* job);
    virtual void applyLoads(Job* job);
};

/*----------------------------------------------------------------------------*/

class ParallelExplicitUSL : public ExplicitUSL{
public:
    ParallelExplicitUSL(){
        object_name = "ParallelExplicitUSL"; //set name of object from registry
    }

    int cpdi_spec = 1;
    int contact_spec = Contact::IMPLICIT;
    int num_threads = 1;

    //thread container
    std::vector<std::thread> threads;

    //scalar sparse matrix operations
    static void ssmOperateStoS(const MPMScalarSparseMatrix &S, const Eigen::VectorXd &x, Eigen::VectorXd &lhs, int SPEC, int k_begin, int k_end);
    static void ssmOperateVtoV(const MPMScalarSparseMatrix &S, const KinematicVectorArray &x, KinematicVectorArray &lhs, int SPEC, int k_begin, int k_end);

    //kinematic vector sparse matrix operations
    static void kvsmLeftMultiply(const KinematicVectorSparseMatrix &gradS, const MaterialTensorArray &T, MaterialVectorArray &lhs, int SPEC, int k_begin, int k_end);
    static void kvsmTensorProductT(const KinematicVectorSparseMatrix &gradS, const KinematicVectorArray &x, KinematicTensorArray L, int SPEC, int k_begin, int k_end);

    //parallel scalar add
    static void scalarAdd(const std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& list, Eigen::VectorXd &sum, int i_begin, int i_end);
    static void vectorAddK(const std::vector<KinematicVectorArray>& list, KinematicVectorArray &sum, int i_begin, int i_end);
    static void vectorAddM(const std::vector<MaterialVectorArray>& list, MaterialVectorArray &sum, int i_begin, int i_end);
    static void tensorAdd(const std::vector<KinematicTensorArray>& list, KinematicTensorArray &sum, int i_begin, int i_end);


    virtual void init(Job* job);
    //virtual void step(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    //parallel functions
    void parallelMultiply(const MPMScalarSparseMatrix &S, const Eigen::VectorXd &x, Eigen::VectorXd &lhs, int SPEC);
    void parallelMultiply(const MPMScalarSparseMatrix &S, const KinematicVectorArray &x, KinematicVectorArray &lhs, int SPEC);
    void parallelMultiply(const KinematicVectorSparseMatrix &gradS, const MaterialTensorArray &T, MaterialVectorArray &lhs, int SPEC);
    void parallelMultiply(const KinematicVectorSparseMatrix &gradS, const KinematicVectorArray &x, KinematicTensorArray L, int SPEC);

    //virtual void createMappings(Job *job);
    virtual void mapPointsToNodes(Job* job);
    //virtual void generateContacts(Job* job);
    //virtual void addContacts(Job* job);
    //virtual void generateBoundaryConditions(Job* job);
    //virtual void addBoundaryConditions(Job* job);
    virtual void moveGrid(Job* job);
    virtual void movePoints(Job* job);
    virtual void calculateStrainRate(Job* job);
    virtual void updateDensity(Job* job);
    //virtual void updateStress(Job* job);
    //virtual void generateLoads(Job* job);
    //virtual void applyLoads(Job* job);
};

#endif //MPM_V3_SOLVERS_HPP

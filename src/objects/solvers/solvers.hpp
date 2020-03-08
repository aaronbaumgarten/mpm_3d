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

    int num_threads = 1;
    bool debug = false;

    //thread container
    std::vector<std::thread> threads;

    //scalar sparse matrix operations
    static void ssmOperateStoS(const MPMScalarSparseMatrix &S, const Eigen::VectorXd &x, Eigen::VectorXd &lhs, int SPEC, int k_begin, int k_end);
    static void ssmOperateVtoV(const MPMScalarSparseMatrix &S, const KinematicVectorArray &x, KinematicVectorArray &lhs, int SPEC, int k_begin, int k_end);

    //kinematic vector sparse matrix operations
    static void kvsmLeftMultiply(const KinematicVectorSparseMatrix &gradS, const MaterialTensorArray &T, MaterialVectorArray &lhs, int SPEC, int k_begin, int k_end);
    static void kvsmTensorProductT(const KinematicVectorSparseMatrix &gradS, const KinematicVectorArray &x, KinematicTensorArray &L, int SPEC, int k_begin, int k_end);
    static void kvsmOperateStoV(const KinematicVectorSparseMatrix &gradS, const Eigen::VectorXd &x, KinematicVectorArray &lhs, int SPEC, int k_begin, int k_end);

    //parallel scalar add
    static void scalarAdd(const std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& list, Eigen::VectorXd &sum, int i_begin, int i_end, bool clear);
    static void vectorAddK(const std::vector<KinematicVectorArray>& list, KinematicVectorArray &sum, int i_begin, int i_end, bool clear);
    static void vectorAddM(const std::vector<MaterialVectorArray>& list, MaterialVectorArray &sum, int i_begin, int i_end, bool clear);
    static void tensorAdd(const std::vector<KinematicTensorArray>& list, KinematicTensorArray &sum, int i_begin, int i_end, bool clear);

    //vectorized operations
    static void multiplyKVbyS(KinematicVectorArray &kv, Eigen::VectorXd &s, double scale, KinematicVectorArray &out, int i_begin, int i_end);
    static void divideKVbyS(KinematicVectorArray &kv, Eigen::VectorXd &s, double scale, KinematicVectorArray &out, int i_begin, int i_end);
    static void multiplySbyS(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, int i_begin, int i_end);
    static void multiplyMTbyS(MaterialTensorArray &mt, Eigen::VectorXd &s, double scale, MaterialTensorArray &out, int i_begin, int i_end);
    static void addStoS(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, int i_begin, int i_end);
    static void subtractSfromS(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, int i_begin, int i_end);

    //zero operations
    static void zeroKV(KinematicVectorArray &kv, int i_begin, int i_end);
    static void zeroS(Eigen::VectorXd &s, int i_begin, int i_end);
    static void zeroMT(MaterialTensorArray &mt, int i_begin, int i_end);

    virtual void init(Job* job);
    virtual void step(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    //parallel functions
    virtual void parallelMultiply(const MPMScalarSparseMatrix &S, const Eigen::VectorXd &x, Eigen::VectorXd &lhs, int SPEC, bool clear = true, int memUnitID = 0);
    virtual void parallelMultiply(const MPMScalarSparseMatrix &S, const KinematicVectorArray &x, KinematicVectorArray &lhs, int SPEC, bool clear = true, int memUnitID = 0);
    virtual void parallelMultiply(const KinematicVectorSparseMatrix &gradS, const MaterialTensorArray &T, MaterialVectorArray &lhs, int SPEC, bool clear = true, int memUnitID = 0);
    virtual void parallelMultiply(const KinematicVectorSparseMatrix &gradS, const KinematicVectorArray &x, KinematicTensorArray &L, int SPEC, bool clear = true, int memUnitID = 0);

    virtual void parallelMultiply(KinematicVectorArray &kv, Eigen::VectorXd &s, double scale, KinematicVectorArray &out, bool clear = true);
    virtual void parallelDivide(KinematicVectorArray &kv, Eigen::VectorXd &s, double scale, KinematicVectorArray &out, bool clear = true);
    virtual void parallelMultiply(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, bool clear = true);
    virtual void parallelMultiply(MaterialTensorArray &mt, Eigen::VectorXd &s, double scale, MaterialTensorArray &out, bool clear = true);
    virtual void parallelAdd(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, bool clear = true);
    virtual void parallelSubtract(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, bool clear = true);

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

/*----------------------------------------------------------------------------*/

class ThreadPoolExplicitUSL : public ParallelExplicitUSL{
public:
    ThreadPoolExplicitUSL(){
        object_name = "ThreadPoolExplicitUSL"; //set name of object from registry
    }

    //pointer to job->threadPool
    ThreadPool* jobThreadPool;

    //memory container for parallel operations
    struct parallelMemoryUnit {
        std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>> s;
        std::vector<KinematicVectorArray> kv;
        std::vector<MaterialVectorArray> mv;
        std::vector<KinematicTensorArray> kt;
    };

    //instantiate
    std::vector<parallelMemoryUnit> memoryUnits;

    //scalar sparse matrix operations
    static void ssmOperateStoSwithFlag(const MPMScalarSparseMatrix &S, const Eigen::VectorXd &x, Eigen::VectorXd &lhs, int SPEC, int k_begin, int k_end, volatile bool& done);
    static void ssmOperateVtoVwithFlag(const MPMScalarSparseMatrix &S, const KinematicVectorArray &x, KinematicVectorArray &lhs, int SPEC, int k_begin, int k_end, volatile bool& done);
    static void kvsmOperateStoVwithFlag(const KinematicVectorSparseMatrix &gradS, const Eigen::VectorXd &x, KinematicVectorArray &lhs, int SPEC, int k_begin, int k_end, volatile bool& done);

    //kinematic vector sparse matrix operations
    static void kvsmLeftMultiplywithFlag(const KinematicVectorSparseMatrix &gradS, const MaterialTensorArray &T, MaterialVectorArray &lhs, int SPEC, int k_begin, int k_end, volatile bool& done);
    static void kvsmTensorProductTwithFlag(const KinematicVectorSparseMatrix &gradS, const KinematicVectorArray &x, KinematicTensorArray &L, int SPEC, int k_begin, int k_end, volatile bool& done);

    //parallel scalar add
    static void scalarAddwithFlag(const std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& list, Eigen::VectorXd &sum, int i_begin, int i_end, bool clear, volatile bool& done);
    static void vectorAddKwithFlag(const std::vector<KinematicVectorArray>& list, KinematicVectorArray &sum, int i_begin, int i_end, bool clear, volatile bool& done);
    static void vectorAddMwithFlag(const std::vector<MaterialVectorArray>& list, MaterialVectorArray &sum, int i_begin, int i_end, bool clear, volatile bool& done);
    static void tensorAddwithFlag(const std::vector<KinematicTensorArray>& list, KinematicTensorArray &sum, int i_begin, int i_end, bool clear, volatile bool& done);

    //vectorized operations
    static void multiplyKVbySwithFlag(KinematicVectorArray &kv, Eigen::VectorXd &s, double scale, KinematicVectorArray &out, int i_begin, int i_end, volatile bool& done);
    static void divideKVbySwithFlag(KinematicVectorArray &kv, Eigen::VectorXd &s, double scale, KinematicVectorArray &out, int i_begin, int i_end, volatile bool& done);
    static void multiplySbySwithFlag(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, int i_begin, int i_end, volatile bool& done);
    static void multiplyMTbySwithFlag(MaterialTensorArray &mt, Eigen::VectorXd &s, double scale, MaterialTensorArray &out, int i_begin, int i_end, volatile bool& done);
    static void addStoSwithFlag(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, int i_begin, int i_end, volatile bool& done);
    static void subtractSfromSwithFlag(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, int i_begin, int i_end, volatile bool& done);

    //zero operations
    static void zeroKVwithFlag(KinematicVectorArray &kv, int i_begin, int i_end, volatile bool& done);
    static void zeroSwithFlag(Eigen::VectorXd &s, int i_begin, int i_end, volatile bool& done);
    static void zeroMTwithFlag(MaterialTensorArray &mt, int i_begin, int i_end, volatile bool& done);

    virtual void init(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    //parallel functions
    virtual void parallelMultiply(const MPMScalarSparseMatrix &S, const Eigen::VectorXd &x, Eigen::VectorXd &lhs, int SPEC, bool clear = true, int memUnitID = 0);
    virtual void parallelMultiply(const MPMScalarSparseMatrix &S, const KinematicVectorArray &x, KinematicVectorArray &lhs, int SPEC, bool clear = true, int memUnitID = 0);
    virtual void parallelMultiply(const KinematicVectorSparseMatrix &gradS, const MaterialTensorArray &T, MaterialVectorArray &lhs, int SPEC, bool clear = true, int memUnitID = 0);
    virtual void parallelMultiply(const KinematicVectorSparseMatrix &gradS, const KinematicVectorArray &x, KinematicTensorArray &L, int SPEC, bool clear = true, int memUnitID = 0);

    virtual void parallelMultiply(KinematicVectorArray &kv, Eigen::VectorXd &s, double scale, KinematicVectorArray &out, bool clear = true);
    virtual void parallelDivide(KinematicVectorArray &kv, Eigen::VectorXd &s, double scale, KinematicVectorArray &out, bool clear = true);
    virtual void parallelMultiply(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, bool clear = true);
    virtual void parallelMultiply(MaterialTensorArray &mt, Eigen::VectorXd &s, double scale, MaterialTensorArray &out, bool clear = true);
    virtual void parallelAdd(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, bool clear = true);
    virtual void parallelSubtract(Eigen::VectorXd &a, Eigen::VectorXd &b, double scale, Eigen::VectorXd &out, bool clear = true);

    //virtual void createMappings(Job *job);
    virtual void mapPointsToNodes(Job* job);
    //virtual void generateContacts(Job* job);
    //virtual void addContacts(Job* job);
    //virtual void generateBoundaryConditions(Job* job);
    //virtual void addBoundaryConditions(Job* job);
    virtual void moveGrid(Job* job);
    virtual void movePoints(Job* job);
    virtual void calculateStrainRate(Job* job);
    //virtual void updateDensity(Job* job);
    //virtual void updateStress(Job* job);
    //virtual void generateLoads(Job* job);
    //virtual void applyLoads(Job* job);
};

#endif //MPM_V3_SOLVERS_HPP

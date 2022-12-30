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

/*----------------------------------------------------------------------------*/
//taylor-green vortex error calculator

class TGVErrorSolver : public Solver{
public:
    TGVErrorSolver(){
        object_name = "TGVErrorSolver"; //set name of object from registry
    }

    //domain should be square and periodic
    double u_max = 0;
    double a = 2.0*M_PI;
    double density = 1000;

    //need to create and write an output file
    std::string output_filename;

    //for some math, need domain size
    KinematicVector Lx;

    //for other math I need exact solution coefficients
    KinematicVectorArray a_M;

    virtual void init(Job* job);
    virtual void step(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    virtual KinematicVector getVelocity(Job* job, KinematicVector const &x);
    virtual KinematicVector getAcceleration(Job* job, KinematicVector const &x);
    virtual double getPressure(Job* job, KinematicVector const &x);
    virtual void setGridLengths(Job* job);
    virtual void setSolutionCoeffs(Job* job);

    virtual void createMappings(Job* job);
    virtual void assignPressure(Job* job);
    virtual void calculateAcceleration(Job* job);
    virtual void assignVelocity(Job* job);
    virtual void movePoints(Job* job);
    virtual void calculateStrainRate(Job* job);
    virtual void generateBoundaryConditions(Job* job); //wrap points
    virtual void updateDensity(Job* job);
    virtual void writeErrorInfo(Job* job);
};

class ShiftedTGVErrorSolver : public TGVErrorSolver{
public:
    ShiftedTGVErrorSolver(){
        object_name = "ShiftedTGVErrorSolver"; //set name of object from registry
    }

    virtual KinematicVector getVelocity(Job* job, KinematicVector const &x);
    virtual KinematicVector getAcceleration(Job* job, KinematicVector const &x);
    virtual double getPressure(Job* job, KinematicVector const &x);

    virtual void calculateAcceleration(Job* job);
};

class GeneralizedVortexErrorSolver : public Solver{
public:
    GeneralizedVortexErrorSolver(){
        object_name = "GeneralizedVortexErrorSolver"; //set name of object from registry
    }

    //domain should be 3m x 3m w/ 3m x 3m solid square
    double a = 0.75;    //inner radii
    double b = 1.25;    //outer radii
    double density = 1000;

    //material parameters
    double E = 1e3;
    double nu = 0.3;
    double G = E / (2.0 * (1.0 + nu));
    double K = E / (3.0 * (1.0 - 2 * nu));

    //simulation parameters
    double A = 1.0;
    double B = 1.0;

    //need to create and write an output file
    std::string output_filename;

    //for some math, need domain size
    double x0 = 1.5;
    double y0 = 1.5;
    KinematicVector Lx;

    virtual void init(Job* job);
    virtual void step(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    virtual double alpha(double r, double t);
    virtual double dalpha_dr(double r, double t);
    virtual double d2alpha_dr2(double r, double t);
    virtual double dalpha_dt(double r, double t);
    virtual double d2alpha_dt2(double r, double t);

    virtual KinematicVector getDisplacement(Job* job, KinematicVector const &x);
    virtual KinematicVector getVelocity(Job* job, KinematicVector const &x);
    virtual KinematicVector getAcceleration(Job* job, KinematicVector const &x);
    virtual KinematicVector getBodyForce(Job* job, KinematicVector const &x);
    virtual MaterialTensor getStress(Job* job, KinematicVector const &x);
    virtual void setGridLengths(Job* job);
    virtual void setInitialVelocity(Job* job);

    virtual void createMappings(Job *job);
    virtual void mapPointsToNodes(Job* job);
    virtual void addBoundaryConditions(Job* job);
    virtual void moveGrid(Job* job);
    virtual void movePoints(Job* job);
    virtual void calculateStrainRate(Job* job);
    virtual void updateDensity(Job* job);
    virtual void updateStress(Job* job);
    virtual void addLoads(Job* job);
    virtual void writeErrorInfo(Job* job);
};

class ExplicitUSLwithVolumetricStrainSmoothing : public ExplicitUSL{
public:
    ExplicitUSLwithVolumetricStrainSmoothing(){
        object_name = "ExplicitUSLwithVolumetricStrainSmoothing"; //set name of object from registry
    }

    //booleans for different methods
    bool use_zhang_cells = true;
    bool anti_locking_only = false;
    bool use_mast_cells = false;
    bool use_zhang_flux = false;

    //string for writing velocity magnitude and error metrics
    bool write_error_metrics = false;
    std::string output_file;
    int skip_counter, stride;

    //zhang flux calculations:
    //only for 2D problems, need to create artificial cells
    KinematicVector x_min, x_max;
    KinematicVector Lx, hx;
    Eigen::VectorXi Nx;
    Eigen::VectorXd vc, v0c; //cell volumes
    Eigen::VectorXd dc, ddc; //cell strains
    Eigen::VectorXd Sc; //cell-wise strain derivative (one direction at a time)

    virtual void init(Job* job);
    virtual void updateDensity(Job* job);
    virtual void updateStress(Job* job);
};

/*--------------------------------------------*/
// Solver for Impact Into Mixtures

class BarotropicViscousFluid;
class TillotsonEOSFluid;
class CompressibleBreakageMechanicsSand;

class ExplicitMixtureSolver : public ExplicitUSL{
public:
    ExplicitMixtureSolver(){
        object_name = "ExplicitMixtureSolver"; //set name of object from registry
    }

    // input parameters
    double grains_d, rhof_0;

    // bodies
    Body *fluid_body;
    Body *granular_body;
    Body *solid_body;

    // fluid material model pointers
    BarotropicViscousFluid *barotropic_viscous_fluid_model;
    TillotsonEOSFluid *tillotson_eos_fluid_model;

    // sand material model pointers
    CompressibleBreakageMechanicsSand *compressible_breakage_mechanics_sand_model;

    // porosity field
    Eigen::VectorXd n;

    // fluid density
    Eigen::VectorXd rho_f;

    // solid "true" density rate of change in solid material frame
    Eigen::VectorXd drhos_dt, drhof_dt;

    // fluid simulation flag
    // 0 -- BarotropicViscousFluid
    // 1 -- TillotsonEOSFluid
    int fluid_model = 0;

    // impact simulation flags
    bool use_reflected_boundary = false;
    bool is_adiabatic = true;
    bool is_compressible = true;

    // reflected boundary information
    KinematicVector Lx;
    bool reflected_boundary_initialized = false;

    // standard functions
    virtual void init(Job* job);
    virtual void step(Job* job);

    // adjusted solver functions
    void addInteractionForces(Job* job);
    virtual void updateDensity(Job* job);           //+ porosity and true fluid density
    virtual void updateStress(Job* job);            //+ fluid density update

};

#endif //MPM_V3_SOLVERS_HPP

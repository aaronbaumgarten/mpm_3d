//
// Created by aaron on 12/12/19.
// fvm_grids.hpp
//

#ifndef MPM_V3_FVM_GRIDS_HPP
#define MPM_V3_FVM_GRIDS_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <fstream>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"
#include "mpm_objects.hpp"
#include "fvm_objects.hpp"

#include "threadpool.hpp"

/*
//finite volume grid class
//responsible for tracking element, face definitions and cross integrals with mpm grid
class FiniteVolumeGrid: public MPMObject{
public:

    //boundary condition types
    /*
    static const int DIRICHLET = 0;
    static const int NEUMANN = 1;
    static const int PERIODIC = 2;
    static const int NEUMANN_DAMPING = 3;
    static const int SUPERSONIC_INLET = 4;


    //boundary conditions
    static const int VELOCITY_INLET         = 0;    //dirichlet
    static const int VELOCITY_TEMP_INLET    = 1;    //dirichlet
    static const int VELOCITY_DENSITY_INLET = 2;    //dirichlet
    static const int PRESSURE_INLET         = 3;
    static const int PRESSURE_OUTLET        = 4;
    static const int DAMPED_OUTLET          = 5;
    static const int ADIABATIC_WALL         = 6;    //dirichlet
    static const int THERMAL_WALL           = 7;    //dirichlet
    static const int SYMMETRIC_WALL         = 8;    //dirichlet
    static const int SUPERSONIC_INLET       = 9;    //dirichlet
    static const int SUPERSONIC_OUTLET      = 10;
    static const int PERIODIC               = 11;
    static const int DAMPED_WALL            = 12;

    //Harten entropy correction scale
    static constexpr double delta = 0.1;

    //grid definions
    int face_count;      //number of faces which define grid
    int node_count;      //number of nodes which define grid
    int element_count;   //number of elements in grid
    int int_quad_count;  //number of interior quadrature points
    int ext_quad_count;  //number of boundary/face quadrature points
    int qpe, qpf;        //number of quadrature points per element/face

    int GRID_DIM = -1; //dimension of grid (might be different than simulation dimension)

    MPMScalarSparseMatrix M = MPMScalarSparseMatrix(0,0); //M_ji maps jth element to ith mpm node
    MPMScalarSparseMatrix Q = MPMScalarSparseMatrix(0,0); //Q_ji maps ith mpm node to jth quadrature point
    KinematicVectorSparseMatrix gradQ = KinematicVectorSparseMatrix(0,0); //gradQ_ji maps ith node gradient to jth quadrature point

    virtual void init(Job*, FiniteVolumeDriver*) = 0; //initialize from mpm Job

    virtual void writeHeader(std::ofstream&, int) = 0; //write cell types

    virtual double getElementVolume(int) = 0;      //return element volume of id
    virtual int getElementTag(int) = 0;            //return 'tag' of id
    virtual double getFaceArea(int) = 0;           //return face area of id
    virtual int getFaceTag(int) = 0;               //return 'tag of face id
    virtual KinematicVector getFaceNormal(Job*, int) = 0; //return normal vector associated with given id

    virtual std::vector<int> getElementQuadraturePoints(int e) = 0; //return vector with global quadrature indices
    virtual std::vector<int> getFaceQuadraturePoints(int f) = 0;    //return vector with global quadrature indices

    virtual std::vector<int> getElementFaces(int) = 0; //return faces associated with given element id
    virtual std::array<int,2> getOrientedElementsByFace(int) = 0; //return elements associated with given face id A -|-> B
    virtual std::vector<int> getElementNeighbors(int) = 0; //return list of elements in neighborhood
    virtual KinematicVector getElementCentroid(Job*, int) = 0; //return element centroid
    virtual KinematicVector getFaceCentroid(Job*, int) = 0;
    virtual KinematicVector getQuadraturePosition(Job* job, int q) = 0;
    virtual double getQuadratureWeight(int q) = 0;

    virtual void generateMappings(Job*, FiniteVolumeDriver*) = 0; //generate M, Q, gradQ maps
    virtual void mapMixturePropertiesToQuadraturePoints(Job*, FiniteVolumeDriver*) = 0; //map porosity and solid velocity to quadrature points
    virtual void constructMomentumField(Job*, FiniteVolumeDriver*) = 0; //construct momentum from FV body
    virtual void constructDensityField(Job*, FiniteVolumeDriver*) = 0; //construct momentum from FV body
    virtual void constructEnergyField(Job*, FiniteVolumeDriver*) = 0; //construct energy field from FV body
    virtual KinematicTensorArray getVelocityGradients(Job*, FiniteVolumeDriver*) = 0; //return velocity gradient in each element

    //functions to compute element mass flux
    virtual Eigen::VectorXd calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver) = 0;

    //functions to compute element momentum fluxes
    virtual KinematicVectorArray calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver) = 0;

    //function to compute element energy fluxes
    virtual Eigen::VectorXd calculateElementEnergyFluxes(Job* job, FiniteVolumeDriver* driver) = 0;

    //function to calculate interphase force
    virtual KinematicVectorArray calculateInterphaseForces(Job* job, FiniteVolumeDriver* driver) = 0;

    //interphase force functions when f_d_e != M*f_d_i
    virtual void calculateSplitIntegralInterphaseForces(Job* job,
                                                        FiniteVolumeDriver* driver,
                                                        KinematicVectorArray& f_i,
                                                        KinematicVectorArray& f_e) = 0;

    virtual void calculateSplitIntegralBuoyantForces(Job* job,
                                                     FiniteVolumeDriver* driver,
                                                     KinematicVectorArray& f_i,
                                                     KinematicVectorArray& f_e) = 0;

    virtual void calculateSplitIntegralDragForces(Job* job,
                                                  FiniteVolumeDriver* driver,
                                                  KinematicVectorArray& f_i,
                                                  KinematicVectorArray& f_e) = 0;

    virtual void calculateSplitIntegralCorrectedDragForces(Job* job,
                                                            FiniteVolumeDriver* driver,
                                                            KinematicVectorArray& f_i,
                                                            KinematicVectorArray& f_e,
                                                            const Eigen::VectorXd &K_n) = 0;

    virtual Eigen::VectorXd calculateInterphaseEnergyFlux(Job* job, FiniteVolumeDriver* driver) = 0;
    virtual Eigen::VectorXd calculateInterphaseEnergyFluxUsingNodeBasedDrag(Job* job,
                                                                            FiniteVolumeDriver* driver,
                                                                            const KinematicVectorArray &f_d) = 0;
}
*/

namespace FiniteVolumeMethod {
    struct BCContainer{
        int tag;
        std::array<double,2> values;
        KinematicVector vector;
    };
}

class FVMGridBase : public FiniteVolumeGrid {
public:
    FVMGridBase() {
        object_name = "FVMGridBase";
    }

    //counter for number of threads
    int num_threads = 1;

    //pointer to job->threadPool
    ThreadPool *jobThreadPool;

    //memory container for parallel operations
    struct parallelMemoryUnit {
        std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>> s;
        std::vector<KinematicVectorArray> kv;
        std::vector<MaterialVectorArray> mv;
        std::vector<KinematicTensorArray> kt;
    };

    //booleans for additional simulation options
    bool USE_REDUCED_QUADRATURE = false;
    bool USE_LOCAL_GRADIENT_CORRECTION = false;

    //instantiate
    std::vector<parallelMemoryUnit> memoryUnits;

    //damping coefficient
    double damping_coefficient = 0.999;

    //grid definitions
    Eigen::VectorXd face_areas;
    KinematicVectorArray face_normals; //x,y,z only

    KinematicVectorArray x_e, x_f, x_q;  //element centroids, face_centroid
    Eigen::VectorXd w_q;            //quadrature weights
    Eigen::VectorXd v_e;            //element volumes
    std::vector<bool> q_b;          //flag for bounding quadrature point or interior quadrature point

    //quadrature point porosity and solid velocity
    Eigen::VectorXd n_q;
    KinematicVectorArray gradn_q, v_sq;

    //boundary conditions
    std::vector<FiniteVolumeMethod::BCContainer> bc_info;

    //grid definitions
    std::vector<std::array<int, 2>> face_elements;   //oriented face list (A -|-> B)
    std::vector<std::vector<int>> element_faces;    //unordered list of faces
    std::vector<std::vector<int>> element_neighbors; //list of nearest neighbors

    //gradient reconstruction vectors for solving least squares problem
    std::vector<Eigen::MatrixXd> A_e;
    std::vector<Eigen::MatrixXd> A_inv; //psuedo inverse
    std::vector<Eigen::VectorXd> b_e;

    virtual void init(Job *job, FiniteVolumeDriver *driver);

    virtual void writeHeader(std::ofstream &file, int TYPE) = 0;

    virtual void generateMappings(Job *job, FiniteVolumeDriver *driver);

    virtual void constructMomentumField(Job *job, FiniteVolumeDriver *driver) = 0;

    virtual void constructDensityField(Job *job, FiniteVolumeDriver *driver) = 0;

    virtual void constructEnergyField(Job *job, FiniteVolumeDriver *driver) = 0;

    virtual void constructPorosityField(Job *, FiniteVolumeDriver *) = 0;

    virtual double getElementVolume(int e);

    virtual int getElementTag(int e);

    virtual double getFaceArea(int f);

    virtual int getFaceTag(int f);

    virtual KinematicVector getFaceNormal(Job *job, int f);

    virtual std::vector<int> getElementQuadraturePoints(int e);

    virtual std::vector<int> getFaceQuadraturePoints(int f);

    virtual std::vector<int> getElementFaces(int e);

    virtual std::array<int, 2> getOrientedElementsByFace(int e);

    virtual std::vector<int> getElementNeighbors(int e);

    virtual KinematicVector getElementCentroid(Job *job, int e); //return element centroid
    virtual KinematicVector getFaceCentroid(Job *job, int f);

    virtual KinematicVector getQuadraturePosition(Job *job, int q);

    virtual double getQuadratureWeight(int q);

    virtual void mapMixturePropertiesToQuadraturePoints(Job *job, FiniteVolumeDriver *driver);

    virtual KinematicTensorArray getVelocityGradients(Job *job, FiniteVolumeDriver *driver);

    //functions to compute element mass flux
    virtual Eigen::VectorXd calculateElementMassFluxes(Job *job, FiniteVolumeDriver *driver);
    virtual Eigen::VectorXd calculateElementMassFluxes(Job *job, FiniteVolumeDriver *driver, int f_start, int f_end);

    //functions to compute element momentum fluxes
    virtual KinematicVectorArray calculateElementMomentumFluxes(Job *job, FiniteVolumeDriver *driver);
    virtual KinematicVectorArray calculateElementMomentumFluxes(Job *job, FiniteVolumeDriver *driver, int f_start, int f_end);

    //function to comput element energy fluxes
    virtual Eigen::VectorXd calculateElementEnergyFluxes(Job *job, FiniteVolumeDriver *driver);
    virtual Eigen::VectorXd calculateElementEnergyFluxes(Job *job, FiniteVolumeDriver *driver, int f_start, int f_end);

    //function to calculate interphase force
    virtual KinematicVectorArray calculateInterphaseForces(Job *job, FiniteVolumeDriver *driver);
    virtual KinematicVectorArray calculateBuoyantForces(Job* job, FiniteVolumeDriver* driver);
    virtual KinematicVectorArray calculateDragForces(Job* job, FiniteVolumeDriver* driver);
    virtual void calculateElementIntegrandsForInterphaseForce(Job *job,
                                                         FiniteVolumeDriver *driver,
                                                         KinematicVectorArray& kv,
                                                         Eigen::VectorXd& v,
                                                         int e_start, int e_end);

    virtual void calculateFaceIntegrandsForInterphaseForce(Job *job,
                                                         FiniteVolumeDriver *driver,
                                                         KinematicVectorArray& kv,
                                                         Eigen::VectorXd& v,
                                                         int f_start, int f_end);

    virtual void calculateElementIntegrandsForBuoyantForce(Job *job,
                                                              FiniteVolumeDriver *driver,
                                                              Eigen::VectorXd& v,
                                                              int e_start, int e_end);

    virtual void calculateFaceIntegrandsForBuoyantForce(Job *job,
                                                           FiniteVolumeDriver *driver,
                                                           KinematicVectorArray& kv,
                                                           int f_start, int f_end,
                                                           bool BOUNDARY_ONLY = true);

    virtual void calculateElementIntegrandsForDragForce(Job *job,
                                                           FiniteVolumeDriver *driver,
                                                           KinematicVectorArray& kv,
                                                           int e_start, int e_end);


    virtual KinematicVectorArray calculateCorrectedDragForces(Job *job,
                                                              FiniteVolumeDriver *driver,
                                                              const Eigen::VectorXd &K_n); // <- drag coefficient at each quad point

    virtual Eigen::VectorXd getCorrectedDragCoefficients(Job* job, FiniteVolumeDriver* driver); // <- function to collect those coefficients

    virtual void calculateElementIntegrandsForCorrectedDragForces(Job *job,
                                                                  FiniteVolumeDriver *driver,
                                                                  KinematicVectorArray &kv,
                                                                  const Eigen::VectorXd &K_n,
                                                                  int e_start, int e_end);           // <- integrand calculator


    //interphase force functions when f_d_e != M*f_d_i
    virtual void calculateSplitIntegralInterphaseForces(Job* job,
                                                        FiniteVolumeDriver* driver,
                                                        KinematicVectorArray& f_i,
                                                        KinematicVectorArray& f_e);

    virtual void calculateSplitIntegralBuoyantForces(Job* job,
                                                     FiniteVolumeDriver* driver,
                                                     KinematicVectorArray& f_i,
                                                     KinematicVectorArray& f_e);

    virtual void calculateSplitIntegralDragForces(Job* job,
                                                  FiniteVolumeDriver* driver,
                                                  KinematicVectorArray& f_i,
                                                  KinematicVectorArray& f_e);

    virtual void calculateSplitIntegralCorrectedDragForces(Job* job,
                                                           FiniteVolumeDriver* driver,
                                                           KinematicVectorArray& f_i,
                                                           KinematicVectorArray& f_e,
                                                           const Eigen::VectorXd &K_n);


    //interphase energy flux functions
    virtual Eigen::VectorXd calculateInterphaseEnergyFlux(Job* job, FiniteVolumeDriver* driver);
    virtual Eigen::VectorXd calculateInterphaseEnergyFluxUsingNodeBasedDrag(Job* job,
                                                                            FiniteVolumeDriver* driver,
                                                                            const KinematicVectorArray &f_d);

    virtual Eigen::VectorXd calculateElementIntegrandsForInterphaseEnergyFlux(Job* job,
                                                                              FiniteVolumeDriver* driver,
                                                                              Eigen::VectorXd &result,
                                                                              int e_begin, int e_end);

    virtual Eigen::VectorXd calculateFaceIntegrandsForInterphaseEnergyFlux(Job* job,
                                                                           FiniteVolumeDriver* driver,
                                                                           Eigen::VectorXd &result,
                                                                           int f_begin, int f_end);

    //parallel functions
    virtual void
    parallelMultiply(const MPMScalarSparseMatrix &S, const Eigen::VectorXd &x, Eigen::VectorXd &lhs, int SPEC,
                     bool clear = true, int memUnitID = 0);

    virtual void
    parallelMultiply(const MPMScalarSparseMatrix &S, const KinematicVectorArray &x, KinematicVectorArray &lhs, int SPEC,
                     bool clear = true, int memUnitID = 0);

    virtual void
    parallelMultiply(const KinematicVectorSparseMatrix &gradS, const Eigen::VectorXd &x, KinematicVectorArray &lhs,
                     int SPEC, bool clear = true, int memUnitID = 0);

    //parallel energy flux functions
    static void calcMassFluxes(Job *job,
                               FiniteVolumeDriver *driver,
                               FVMGridBase *grid,
                               Eigen::VectorXd &lhs,
                               int k_begin, int k_end,
                               volatile bool &done);

    static void calcMomentumFluxes(Job *job,
                               FiniteVolumeDriver *driver,
                                   FVMGridBase *grid,
                               KinematicVectorArray &lhs,
                               int k_begin, int k_end,
                               volatile bool &done);

    static void calcEnergyFluxes(Job *job,
                               FiniteVolumeDriver *driver,
                                 FVMGridBase *grid,
                               Eigen::VectorXd &lhs,
                               int k_begin, int k_end,
                               volatile bool &done);

    static void calcElementIntegrandForInterphaseForces(Job *job,
                                                        FiniteVolumeDriver *driver,
                                                        FVMGridBase *grid,
                                                        KinematicVectorArray &kv,
                                                        Eigen::VectorXd& v,
                                                        int k_begin, int k_end,
                                                        volatile bool &done);


    static void calcFaceIntegrandForInterphaseForces(Job *job,
                                                     FiniteVolumeDriver *driver,
                                                     FVMGridBase *grid,
                                                     KinematicVectorArray &kv,
                                                     Eigen::VectorXd& v,
                                                     int k_begin, int k_end,
                                                     volatile bool &done);

    static void calcElementIntegrandsForInterphaseEnergyFlux(Job* job,
                                                             FiniteVolumeDriver* driver,
                                                             FVMGridBase* grid,
                                                             Eigen::VectorXd& v,
                                                             int k_begin, int k_end,
                                                             volatile bool &done);


    static void calcFaceIntegrandsForInterphaseEnergyFlux(Job* job,
                                                             FiniteVolumeDriver* driver,
                                                             FVMGridBase* grid,
                                                             Eigen::VectorXd& v,
                                                             int k_begin, int k_end,
                                                             volatile bool &done);
};

class FVMCartesian : public FVMGridBase{
public:
    FVMCartesian(){
        object_name = "FVMCartesian";
    }

    //grid dimensions
    std::vector<int> Nx;
    std::vector<double> Lx, hx;

    //mapping functions
    std::vector<int> e_to_ijk(int e); //convert element id to ijk element definition
    int ijk_to_e(std::vector<int> ijk); //convert ijk element definition to e
    std::vector<int> n_to_ijk(int n); //convert node id to ijk
    int ijk_to_n(std::vector<int> ijk); //convert ijk to node id

    virtual void init(Job* job, FiniteVolumeDriver* driver);

    virtual void writeHeader(std::ofstream& file, int TYPE);

    virtual void constructMomentumField(Job* job, FiniteVolumeDriver* driver);
    virtual void constructDensityField(Job* job, FiniteVolumeDriver* driver);
    virtual void constructEnergyField(Job* job, FiniteVolumeDriver* driver);
    virtual void constructPorosityField(Job* job, FiniteVolumeDriver* driver);
};

/*----------------------------------------------------------------------------*/

class FVMGmsh2D : public FVMGridBase{
public:
    FVMGmsh2D(){
        object_name = "FVMGmsh2D";
    }

    std::string filename;
    Eigen::MatrixXi nodeIDs;        //element to node map
    int npe = 3;                    //nodes per element
    KinematicVectorArray x_n;       //node positions

    //grid definitions
    std::vector<std::array<int,2>> face_nodes;      //un-oriented list of nodes a---b

    virtual void init(Job* job, FiniteVolumeDriver* driver);

    virtual void writeHeader(std::ofstream& file, int TYPE);
    virtual void constructMomentumField(Job* job, FiniteVolumeDriver* driver);
    virtual void constructDensityField(Job* job, FiniteVolumeDriver* driver);
    virtual void constructEnergyField(Job* job, FiniteVolumeDriver* driver);
    virtual void constructPorosityField(Job* job, FiniteVolumeDriver* driver);
};

#endif //MPM_V3_FVM_GRIDS_HPP

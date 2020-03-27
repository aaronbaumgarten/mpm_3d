//
// Created by aaron on 12/23/19.
// fvm_solvers.hpp
//

#ifndef MPM_V3_FVM_SOLVERS_HPP
#define MPM_V3_FVM_SOLVERS_HPP

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

/*
class FiniteVolumeSolver : public MPMObject{
public:
    //functions that must be implemented by every finite volume solver (not many)
    virtual void init(Job*, FiniteVolumeDriver*) = 0;                                        //initialize from Job
    virtual void step(Job*, FiniteVolumeDriver*) = 0;                                        //perform single mpm step
};
 */

class FVMDefaultSolver : public FiniteVolumeSolver{
public:
    FVMDefaultSolver(){
        object_name = "FVMDefaultSolver";
    }

    Eigen::VectorXd density_fluxes;
    KinematicVectorArray momentum_fluxes;
    Eigen::VectorXd energy_fluxes;

    virtual void init(Job* job, FiniteVolumeDriver* driver);                                        //initialize from Job
    virtual void step(Job* job, FiniteVolumeDriver* driver);                                        //perform single mpm step
};

class FVMRungeKuttaSolver : public FiniteVolumeSolver {
public:
    FVMRungeKuttaSolver(){
        object_name = "FVMRungeKuttaSolver";
    }

    int order = 4;

    Eigen::VectorXd density_fluxes;
    KinematicVectorArray momentum_fluxes;
    Eigen::VectorXd energy_fluxes;

    Eigen::VectorXd rho_k1, rho_k2, rho_k3, rho_k4;
    KinematicVectorArray p_k1, p_k2, p_k3, p_k4;
    Eigen::VectorXd rhoE_k1, rhoE_k2, rhoE_k3, rhoE_k4;

    virtual void init(Job* job, FiniteVolumeDriver* driver);                                        //initialize from Job
    virtual void step(Job* job, FiniteVolumeDriver* driver);                                        //perform single mpm step
};


class FVMSteadyStateSolver : public FiniteVolumeSolver{
public:
    FVMSteadyStateSolver(){
        object_name = "FVMSteadyStateSolver";
    }

    static const int BICGSTAB = 0;
    static const int DIRECT = 1;
    static const int JACOBI = 2;

    Eigen::VectorXd density_fluxes, tmp_df;
    KinematicVectorArray momentum_fluxes, tmp_mf;
    Eigen::VectorXd energy_fluxes, tmp_ef;

    Eigen::VectorXd u_0; //state vector for sth psuedo-transient continuation
    Eigen::VectorXd u_n; //state vector for nth newton step
    Eigen::VectorXd r_n; //residual vector for nth newton step
    double h = 1e-5;     //step length for numerical derivatives
    double nu = 0;     //iterative regularizaton
    double REL_TOL = 1e-10;
    double ABS_TOL = 1e-7;
    int max_i = 100;
    int max_n = 100;
    int max_s = 1;
    int INNER_SOLVER = 0;

    double BiCGSTAB_TOL = 100;
    double LineSearch_TOL = 1e-3;

    int vector_size = 0;

    virtual void init(Job* job, FiniteVolumeDriver* driver);                                        //initialize from Job
    virtual void step(Job* job, FiniteVolumeDriver* driver);                                        //perform single step

    Eigen::VectorXd F(Job* job, FiniteVolumeDriver* driver, const Eigen::VectorXd& u, double nu_i);          //calculate flux function value
    Eigen::VectorXd DF(Job* job, FiniteVolumeDriver* driver, const Eigen::VectorXd& x_i, double h_i);        //calculate flux derivative in x direction

    void convertVectorToStateSpace(Job* job, FiniteVolumeDriver* driver,
                                   const Eigen::VectorXd& v,
                                   Eigen::VectorXd& rho,
                                   KinematicVectorArray& p,
                                   Eigen::VectorXd& rhoE);


    void convertStateSpaceToVector(Job* job, FiniteVolumeDriver* driver,
                                   Eigen::VectorXd& v,
                                   const Eigen::VectorXd& rho,
                                   const KinematicVectorArray& p,
                                   const Eigen::VectorXd& rhoE);
};

class FVMMixtureSolver : public FiniteVolumeSolver{
public:
    FVMMixtureSolver(){
        object_name = "FVMMixtureSolver";
    }

    int cpdi_spec = 1;
    int contact_spec = Contact::IMPLICIT;
    int solid_body_id = -1;

    Eigen::VectorXd density_fluxes;
    KinematicVectorArray momentum_fluxes;
    Eigen::VectorXd energy_fluxes;

    KinematicVectorArray f_i;
    KinematicVectorArray f_e;

    virtual void init(Job* job, FiniteVolumeDriver* driver);                                        //initialize from Job
    virtual void step(Job* job, FiniteVolumeDriver* driver);                                        //perform single mpm step

    virtual void createMappings(Job *job);
    virtual void mapPointsToNodes(Job* job);
    virtual void calculateElementGradients(Job* job, FiniteVolumeDriver* driver);
    virtual void mapMixturePropertiesToElements(Job* job, FiniteVolumeDriver* driver);
    virtual void generateFluxes(Job* job, FiniteVolumeDriver* driver);
    virtual void applyFluxes(Job* job, FiniteVolumeDriver* driver);
    virtual void generateMixtureForces(Job* job, FiniteVolumeDriver* driver);
    virtual void applyMixtureForces(Job* job, FiniteVolumeDriver* driver);
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

class FVMMixtureSolverRK4 : public FVMMixtureSolver{
public:
    FVMMixtureSolverRK4(){
        object_name = "FVMMixtureSolverRK4";
    }

    int vector_size = 0;
    Eigen::VectorXd u_0, u_n, k1, k2, k3, k4;
    KinematicVectorArray f_i1, f_i2, f_i3, f_i4;
    KinematicVectorArray f_d, f_b, f_d_e;
    Eigen::VectorXd K_n; //storage vector for drag coefficients
    KinematicVectorArray first_drag_correction; //drag correction component from t = t_n
    KinematicVectorArray second_drag_correction; //drag correction component form t = t_star;

    void adjustSolidVelocity(Job* job, FiniteVolumeDriver* driver, double h);
    Eigen::VectorXd F(Job* job, FiniteVolumeDriver* driver, const Eigen::VectorXd& u);          //calculate flux function value

    void convertVectorToStateSpace(Job* job, FiniteVolumeDriver* driver,
                                   const Eigen::VectorXd& v,
                                   Eigen::VectorXd& rho,
                                   KinematicVectorArray& p,
                                   Eigen::VectorXd& rhoE);


    void convertStateSpaceToVector(Job* job, FiniteVolumeDriver* driver,
                                   Eigen::VectorXd& v,
                                   const Eigen::VectorXd& rho,
                                   const KinematicVectorArray& p,
                                   const Eigen::VectorXd& rhoE);

    virtual void init(Job* job, FiniteVolumeDriver* driver);                                        //initialize from Job
    virtual void step(Job* job, FiniteVolumeDriver* driver);                                        //perform single mpm step


    virtual void updateDensity(Job* job, int SPEC); //overload to avoid updating internal state
    virtual void updateStress(Job* job, int SPEC); //overload function
};

#endif //MPM_V3_FVM_SOLVERS_HPP

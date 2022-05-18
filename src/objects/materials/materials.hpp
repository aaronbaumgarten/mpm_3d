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

class CompressibleNeohookeanElasticity : public IsotropicLinearElasticity {
public:
    CompressibleNeohookeanElasticity(){
        object_name = "CompressibleNeohookeanElasticity";
    }

    //material properties
    MaterialTensorArray F;

    void init(Job* job, Body* body);
    void calculateStress(Job* job, Body* body, int SPEC);
    void assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC);
    void assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
};

/*----------------------------------------------------------------------------*/

class Sand_SachithLocal : public Material {
public:
    Sand_SachithLocal(){
        object_name = "Sand_SachithLocal";
    }

    //boolean flag for pressure smoothing (i.e. combining with ExplicitUSLwithVolumetricStrainSmoothing)
    //assign pressure by density
    bool use_density_to_determine_pressure = false;

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

/*----------------------------------------------------------------------------*/

class SlurryGranularPhase : public Material{
public:
    SlurryGranularPhase(){
        object_name = "SlurryGranularPhase";
    }

    //residual tolerance
    double ABS_TOL = 1e-6;
    double REL_TOL = 1e-6;
    double h = 1e-5;

    //material properties
    double E, nu, G, K, lambda;

    //mu(Im)
    double mu_1, mu_2, K_3, K_4, K_5, phi_m, grains_rho, eta_0, grains_d, fluid_rho;
    double a, I_0;

    Eigen::VectorXd gammap, gammap_dot, I_v, I, I_m, phi, eta;
    int fluid_body_id = -1;

    virtual double getBeta(double phi, double phi_eq);
    virtual void calcState(double &gdp, double &p, double &eta_in, double &phi_in, double &I_out, double &Iv_out, double &Im_out, double &mu_out, double &phi_eq, double &beta_out);
    virtual double getStep(double val, double ub, double lb);

    virtual void init(Job* job, Body* body);
    virtual void calculateStress(Job* job, Body* body, int SPEC);
    virtual void assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC);
    virtual void assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

    virtual void writeFrame(Job* job, Body* body, Serializer* serializer);
    virtual std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);

    //function for parallel calculation of stresses
    virtual void cSOS(Job* job, Body* body, int SPEC, int k_begin, int k_end);

    static void calculateSubsetOfStresses(Job* job,
                                               Body* body,
                                               SlurryGranularPhase* material,
                                               int SPEC,
                                               int k_begin, int k_end,
                                               volatile bool &done);
};

/*----------------------------------------------------------------------------*/

class BarotropicViscousFluid : public Material {
public:
    BarotropicViscousFluid(){
        object_name = "BarotropicViscousFluid";
    }

    bool use_delta_correction = false;

    double alpha = 0.0;
    double eta, K;

    KinematicVector Lx;
    double h;
    Eigen::VectorXd e;
    KinematicVectorArray grad_e, del_pos;
    Eigen::VectorXd V_i, v_i;

    void init(Job* job, Body* body);
    void calculateStress(Job* job, Body* body, int SPEC);
    void assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC);
    void assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

/*----------------------------------------------------------------------------*/

class SlurryFluidPhase : public Material {
public:
    SlurryFluidPhase(){
        object_name = "SlurryFluidPhase";
    }

    double alpha = 0.0;
    double eta_0, K;
    double solid_rho;
    int solid_body_id = -1;

    Eigen::VectorXd J, n_p;

    KinematicVector Lx;
    double h;
    Eigen::VectorXd e;
    KinematicVectorArray grad_e, del_pos;
    Eigen::VectorXd V_i, v_i;

    void init(Job* job, Body* body);
    void calculateStress(Job* job, Body* body, int SPEC);
    void assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC);
    void assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

/*----------------------------------------------------------------------------*/

class Cornstarch : public SlurryGranularPhase{
public:
    Cornstarch(){
        object_name = "Cornstarch";
    }

    //new variable
    double a_0, a_inf, phi_j, phi_c, K_6, K_7, tau_star;
    double Delta, alpha, phi_star;

    Eigen::VectorXd c, phi_m_vec;

    void calcState(double &gdp, double &p, double &eta_in, double &phi_in, double &I_out, double &Iv_out, double &Im_out, double &mu_out, double &phi_eq, double &beta_out);
    double getBeta(double phi, double phi_eq);
    void init(Job* job, Body* body);
    void calculateStress(Job* job, Body* body, int SPEC);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
};

/*----------------------------------------------------------------------------*/

class SlurryGranularPhase_wUnderCompaction : public SlurryGranularPhase{
public:
    SlurryGranularPhase_wUnderCompaction(){
        object_name = "SlurryGranularPhase_wUnderCompaction";
    }

    //material properties
    double phi_c;

    virtual void init(Job* job, Body* body);
    virtual void calculateStress(Job* job, Body* body, int SPEC);
};

/*----------------------------------------------------------------------------*/
//this is the material that defines a fish's body motion (6.883 spring 2020 project)
class Fish : public Material {
public:
    Fish(){
        object_name = "Fish";
    }

    bool USE_RIGID_TAIL = false;
    bool USE_RIGID_HEAD = false;

    static const int STARBOARD_HEAD = 0;
    static const int STARBOARD_ABDOMEN = 1;
    static const int STARBOARD_TAIL = 2;
    static const int PORT_HEAD = 3;
    static const int PORT_ABDOMEN = 4;
    static const int PORT_TAIL = 5;

    //material properties
    double E, nu, G, K, lambda;
    MaterialTensorArray F, F_e, F_m; //deformation tensors
    std::vector<int> label;           //label for fish parts

    double alpha_h, alpha_a, alpha_t;
    double beta_h, beta_a, beta_t;
    double gamma_h, gamma_a, gamma_t;

    void init(Job* job, Body* body);
    void calculateStress(Job* job, Body* body, int SPEC);
    void assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC);
    void assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};


class BreakageMechanicsSand : public Material {
public:
    BreakageMechanicsSand(){
        object_name = "BreakageMechanicsSand";
    }

    //material properties
    double pr, K, G, Ec, M_0, g, phi_l, phi_u, l, u, theta, rho_0;

    //history variables
    Eigen::VectorXd B, phi;     //breakage measure, internal porosity
    MaterialTensorArray Be;     //elastic left cauchy-green tensor

    //scalar outputs
    Eigen::VectorXd evDot, esDot, BDot;

    //debug variables
    Eigen::VectorXd kVec, yVec;

    //state struct for computation
    struct MaterialState {
        double B;       //Breakage
        double rho;     //Effective Density
        double phi;     //Porosity
        double ev;      //Volumetric Strain (Elastic)
        double es;      //Shear Strain (Elastic)
    };

    //computational crutches apadted from MATLAB
    double CriticalStatePorosityFromMaterialState(MaterialState mat_state);
    double EBFromMaterialState(MaterialState mat_state);
    double FFromMaterialState(MaterialState mat_state);
    std::vector<double> PQFromMaterialState(MaterialState mat_state);
    std::vector<double> RelativePlasticityRatesFromMaterialStateandDeformation(MaterialState mat_state, MaterialTensor Be);
    double YieldFunctionFromMaterialStateandDeformation(MaterialState mat_state, MaterialTensor Be);
    MaterialTensor CauchyStressFromMaterialStateandDeformation(MaterialState mat_state, MaterialTensor Be);

    void init(Job* job, Body* body);
    void calculateStress(Job* job, Body* body, int SPEC);
    void assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC);
    void assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
    std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);
};

#endif //MPM_V3_MATERIALS_HPP

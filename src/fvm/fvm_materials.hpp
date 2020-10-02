//
// Created by aaron on 12/23/19.
// fvm_materials.hpp
//

#ifndef MPM_V3_FVM_MATERIALS_HPP
#define MPM_V3_FVM_MATERIALS_HPP

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
//finite volume material class
//represents a single continuum body filling fluid domain
class FiniteVolumeMaterial : public MPMObject{
public:
    //initialize from job and driver
    virtual void init(Job*, FiniteVolumeDriver*) = 0;

    //functions to calculate fluid fields
    virtual MaterialTensor getStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n) = 0;
    virtual MaterialTensor getShearStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n) = 0;
    virtual double getPressure(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n) = 0;
    virtual double getTemperature(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n) = 0;
    virtual double getSpeedOfSound(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n) = 0;
    virtual void calculateElementPressures(Job*, FiniteVolumeDriver*) = 0;
    virtual void calculateElementShearStresses(Job*, FiniteVolumeDriver*) = 0;
    virtual void calculateElementTemperatures(Job*, FiniteVolumeDriver*) = 0;

    //fluid equations of state
    virtual double getDensityFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n) = 0;
    virtual double getInternalEnergyFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n) = 0;
    virtual double getPressureFromDensityAndTemperature(Job*, FiniteVolumeDriver*, double rho, double theta, double n) = 0;
    virtual KinematicVector getHeatFlux(Job*, FiniteVolumeDriver*, double rho, double theta, const KinematicVector& theta_x, double n) = 0;
    virtual double getSpeedOfSoundFromEnthalpy(Job *, FiniteVolumeDriver *, double rho, const KinematicVector &p, double rhoH, double n) = 0;
    virtual double getPressureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Pt, double M) = 0;
    virtual double getTemperatureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Tt, double M) = 0;

    //mixture model functions
    virtual int calculatePorosity(Job*, FiniteVolumeDriver*) = 0; //return 1 if mixture problem, return 0 if FVM problem only
    virtual int updateSolidPhaseVelocity(Job*, FiniteVolumeDriver*) = 0; //return 1 if mixture problem, return 0 if FVM problem only
    virtual KinematicVector getInterphaseDrag(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n) = 0;
};
 */

class FVMBarotropicViscousFluid : public FiniteVolumeMaterial{
public:
    FVMBarotropicViscousFluid(){
        object_name = "FVMBarotropicViscousFluid";
    }

    //material properties
    double kappa, eta, rho_0;

    //initialize from job and driver
    virtual void init(Job* job, FiniteVolumeDriver* driver);

    //functions to calculate fluid fields
    virtual MaterialTensor getStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n);
    virtual MaterialTensor getShearStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getPressure(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getTemperature(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getSpeedOfSound(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual void calculateElementPressures(Job*, FiniteVolumeDriver*);
    virtual void calculateElementShearStresses(Job*, FiniteVolumeDriver*);
    virtual void calculateElementTemperatures(Job*, FiniteVolumeDriver*);

    //fluid equations of state
    virtual double getDensityFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n);
    virtual double getInternalEnergyFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n);
    virtual double getPressureFromDensityAndTemperature(Job*, FiniteVolumeDriver*, double rho, double theta, double n);
    virtual KinematicVector getHeatFlux(Job*, FiniteVolumeDriver*, double rho, double theta, const KinematicVector& theta_x, double n);
    virtual double getSpeedOfSoundFromEnthalpy(Job *, FiniteVolumeDriver *, double rho, const KinematicVector &p, double rhoH, double n);
    virtual double getPressureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Pt, double M);
    virtual double getTemperatureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Tt, double M);

    //mixture model functions
    virtual int calculatePorosity(Job*, FiniteVolumeDriver*); //return 1 if mixture problem, return 0 if FVM problem only
    virtual int updateSolidPhaseVelocity(Job*, FiniteVolumeDriver*); //return 1 if mixture problem, return 0 if FVM problem only
    virtual KinematicVector getInterphaseDrag(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n);
    virtual double getInterphaseDragCoefficient(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n, int SPEC = REGULAR_DRAG);
};

class FVMSlurryFluidPhase : public FiniteVolumeMaterial{
public:
    FVMSlurryFluidPhase(){
        object_name = "FVMSlurryFluidPhase";
    }

    //material properties
    double kappa, eta, rho_0;

    //solid properties
    double solid_rho, grain_diam;
    int solid_body_id = -1;

    //initialize from job and driver
    virtual void init(Job* job, FiniteVolumeDriver* driver);

    //functions to calculate fluid fields
    virtual MaterialTensor getStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n);
    virtual MaterialTensor getShearStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getPressure(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getTemperature(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getSpeedOfSound(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual void calculateElementPressures(Job*, FiniteVolumeDriver*);
    virtual void calculateElementShearStresses(Job*, FiniteVolumeDriver*);
    virtual void calculateElementTemperatures(Job*, FiniteVolumeDriver*);

    //fluid equations of state
    virtual double getDensityFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n);
    virtual double getInternalEnergyFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n);
    virtual double getPressureFromDensityAndTemperature(Job*, FiniteVolumeDriver*, double rho, double theta, double n);
    virtual KinematicVector getHeatFlux(Job*, FiniteVolumeDriver*, double rho, double theta, const KinematicVector& theta_x, double n);
    virtual double getSpeedOfSoundFromEnthalpy(Job *, FiniteVolumeDriver *, double rho, const KinematicVector &p, double rhoH, double n);
    virtual double getPressureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Pt, double M);
    virtual double getTemperatureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Tt, double M);

    //mixture model functions
    virtual int calculatePorosity(Job*, FiniteVolumeDriver*); //return 1 if mixture problem, return 0 if FVM problem only
    virtual int updateSolidPhaseVelocity(Job*, FiniteVolumeDriver*); //return 1 if mixture problem, return 0 if FVM problem only
    virtual KinematicVector getInterphaseDrag(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n);
    virtual double getInterphaseDragCoefficient(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n, int SPEC = REGULAR_DRAG);
};

class FVMIdealGas : public FiniteVolumeMaterial{
public:
    FVMIdealGas(){
        object_name = "FVMIdealGas";
    }

    //material properties
    double heat_capacity_ratio, R, eta, thermal_conductivity;

    //initialize from job and driver
    virtual void init(Job* job, FiniteVolumeDriver* driver);

    //functions to calculate fluid fields
    virtual MaterialTensor getStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n);
    virtual MaterialTensor getShearStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getPressure(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getTemperature(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getSpeedOfSound(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual void calculateElementPressures(Job*, FiniteVolumeDriver*);
    virtual void calculateElementShearStresses(Job*, FiniteVolumeDriver*);
    virtual void calculateElementTemperatures(Job*, FiniteVolumeDriver*);

    //fluid equations of state
    virtual double getDensityFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n);
    virtual double getInternalEnergyFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n);
    virtual double getPressureFromDensityAndTemperature(Job*, FiniteVolumeDriver*, double rho, double theta, double n);
    virtual KinematicVector getHeatFlux(Job*, FiniteVolumeDriver*, double rho, double theta, const KinematicVector& theta_x, double n);
    virtual double getSpeedOfSoundFromEnthalpy(Job *, FiniteVolumeDriver *, double rho, const KinematicVector &p, double rhoH, double n);
    virtual double getPressureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Pt, double M);
    virtual double getTemperatureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Tt, double M);

    //mixture model functions
    virtual int calculatePorosity(Job*, FiniteVolumeDriver*); //return 1 if mixture problem, return 0 if FVM problem only
    virtual int updateSolidPhaseVelocity(Job*, FiniteVolumeDriver*); //return 1 if mixture problem, return 0 if FVM problem only
    virtual KinematicVector getInterphaseDrag(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n);
    virtual double getInterphaseDragCoefficient(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n, int SPEC = REGULAR_DRAG);
};

class FVMSlurryGasPhase : public FiniteVolumeMaterial{
public:
    FVMSlurryGasPhase(){
        object_name = "FVMSlurryGasPhase";
    }

    //use eddy viscosity correction
    bool USE_EDDY_VISCOSITY = false;
    double a_T, h_T;

    //material properties
    double heat_capacity_ratio, R, eta, thermal_conductivity;

    //solid properties
    double solid_rho, grain_diam;
    int solid_body_id = -1;

    //initialize from job and driver
    virtual void init(Job* job, FiniteVolumeDriver* driver);

    //functions to calculate fluid fields
    virtual MaterialTensor getStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n);
    virtual MaterialTensor getShearStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getPressure(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getTemperature(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual double getSpeedOfSound(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n);
    virtual void calculateElementPressures(Job*, FiniteVolumeDriver*);
    virtual void calculateElementShearStresses(Job*, FiniteVolumeDriver*);
    virtual void calculateElementTemperatures(Job*, FiniteVolumeDriver*);

    //fluid equations of state
    virtual double getDensityFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n);
    virtual double getInternalEnergyFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n);
    virtual double getPressureFromDensityAndTemperature(Job*, FiniteVolumeDriver*, double rho, double theta, double n);
    virtual KinematicVector getHeatFlux(Job*, FiniteVolumeDriver*, double rho, double theta, const KinematicVector& theta_x, double n);
    virtual double getSpeedOfSoundFromEnthalpy(Job *, FiniteVolumeDriver *, double rho, const KinematicVector &p, double rhoH, double n);
    virtual double getPressureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Pt, double M);
    virtual double getTemperatureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Tt, double M);

    //mixture model functions
    virtual int calculatePorosity(Job*, FiniteVolumeDriver*); //return 1 if mixture problem, return 0 if FVM problem only
    virtual int updateSolidPhaseVelocity(Job*, FiniteVolumeDriver*); //return 1 if mixture problem, return 0 if FVM problem only
    virtual KinematicVector getInterphaseDrag(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n);
    virtual double getInterphaseDragCoefficient(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n, int SPEC = REGULAR_DRAG);
};

#endif //MPM_V3_FVM_MATERIALS_HPP

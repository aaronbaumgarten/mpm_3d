//
// Created by aaron on 12/23/19.
// fvm_barotropic_viscous_fluid.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <fstream>
#include <job.hpp>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "mpm_objects.hpp"
#include "fvm_objects.hpp"
#include "fvm_materials.hpp"

/*----------------------------------------------------------------------------*/
//initialize from job and driver
void FVMBarotropicViscousFluid::init(Job* job, FiniteVolumeDriver* driver){
    //check that simulation type is ISOTHERMAL
    if (driver->TYPE != FiniteVolumeDriver::ISOTHERMAL){
        std::cout << "WARNING! FVMBarotropicViscousFluid is an ISOTHERMAL model. Simulation TYPE "
                  << driver->TYPE << " != ISOTHERMAL." << std::endl;
    }

    if (fp64_props.size() < 3) {
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 3 properties defined (kappa, eta, rho_0).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        kappa = fp64_props[0]; //bulk modulus (Pa)
        eta = fp64_props[1];   //shear viscosity (Pa*s)
        rho_0 = fp64_props[2]; //reference density (kg/m^3)
    }
    std::cout << "FiniteVolumeMaterial properties (kappa = " << kappa << " Pa, eta = " << eta << " Pa*s, rho_0 = " << rho_0;
    std::cout << " kg/m^3)." << std::endl;
    std::cout << "FiniteVolumeMaterial initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//functions to calculate fluid fields
//get full stress tensor in fluid
MaterialTensor FVMBarotropicViscousFluid::getStress(Job* job,
                                                    FiniteVolumeDriver* driver,
                                                    const KinematicTensor& L,
                                                    double rho,
                                                    const KinematicVector& p,
                                                    double rhoE,
                                                    double n){
    //sigma = tau - p1
    return eta*(L + L.transpose() - 2.0/3.0*L.trace()*MaterialTensor::Identity())
           - kappa*std::log(rho/rho_0)*MaterialTensor::Identity();
}

//get shear components of stress tensor
MaterialTensor FVMBarotropicViscousFluid::getShearStress(Job* job,
                                                         FiniteVolumeDriver* driver,
                                                         const KinematicTensor& L,
                                                         double rho,
                                                         const KinematicVector& p,
                                                         double rhoE,
                                                         double n){
    //tau = 2*eta*(D - 1/3 trD I)
    return eta*(L + L.transpose() - 2.0/3.0*L.trace()*MaterialTensor::Identity());
}

//get volumetric component of stress tensor
double FVMBarotropicViscousFluid::getPressure(Job* job,
                                              FiniteVolumeDriver* driver,
                                              double rho,
                                              const KinematicVector& p,
                                              double rhoE,
                                              double n){
    //P = kappa*log(rho/rho_0)
    return kappa*std::log(rho/rho_0);
}

//get temperature of fluid from state
double FVMBarotropicViscousFluid::getTemperature(Job* job,
                                                 FiniteVolumeDriver* driver,
                                                 double rho,
                                                 const KinematicVector& p,
                                                 double rhoE,
                                                 double n){
    //isothermal model
    return 0.0;
}

//return speed of sound
double FVMBarotropicViscousFluid::getSpeedOfSound(Job* job,
                                                  FiniteVolumeDriver* driver,
                                                  double rho,
                                                  const KinematicVector& p,
                                                  double rhoE,
                                                  double n){
    //c^2 = kappa/rho
    return std::sqrt(kappa/rho);
}

//loop over elements and fill in pressures
void FVMBarotropicViscousFluid::calculateElementPressures(Job* job, FiniteVolumeDriver* driver){
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        driver->fluid_body->P(e) = kappa*std::log(driver->fluid_body->rho(e)/rho_0);
    }
    return;
}

//loop over elements and fill in shear stresses
void FVMBarotropicViscousFluid::calculateElementShearStresses(Job* job, FiniteVolumeDriver* driver){
    //get strain rates from grid reconstruction (not flux limited)
    KinematicTensorArray L = driver->fluid_body->L; //driver->fluid_grid->getVelocityGradients(job, driver);

    //fill in stresses accordingly
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        driver->fluid_body->tau[e] = eta*(L[e] + L[e].transpose() - 2.0/3.0*L[e].trace()*MaterialTensor::Identity());
    }
    return;
}

//loop over elements and fill in temperatures
void FVMBarotropicViscousFluid::calculateElementTemperatures(Job* job, FiniteVolumeDriver* driver){
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        driver->fluid_body->theta(e) = 0.0;
    }
    return;
}

/*----------------------------------------------------------------------------*/
//fluid equations of state
double FVMBarotropicViscousFluid::getDensityFromPressureAndTemperature(Job* job,
                                                                       FiniteVolumeDriver* driver,
                                                                       double pressure,
                                                                       double theta,
                                                                       double n){
    return rho_0*std::exp(pressure/kappa);
}

double FVMBarotropicViscousFluid::getInternalEnergyFromPressureAndTemperature(Job* job,
                                                                              FiniteVolumeDriver* driver,
                                                                              double pressure,
                                                                              double theta,
                                                                              double n){
    return 0.0;
}

double FVMBarotropicViscousFluid::getPressureFromDensityAndTemperature(Job* job,
                                                                       FiniteVolumeDriver* driver,
                                                                       double rho,
                                                                       double theta,
                                                                       double n){
    //P = kappa*log(rho/rho_0)
    return kappa*std::log(rho/rho_0);
}

KinematicVector FVMBarotropicViscousFluid::getHeatFlux(Job* job,
                                                       FiniteVolumeDriver* driver,
                                                       double rho,
                                                       double theta,
                                                       const KinematicVector& theta_x,
                                                       double n){
    //no heat flux for now...
    return KinematicVector(job->JOB_TYPE);
}

double FVMBarotropicViscousFluid::getSpeedOfSoundFromEnthalpy(Job *job, FiniteVolumeDriver *driver, double rho, const KinematicVector &p, double rhoH,
                                                              double n) {
    std::cerr << "ERROR: FVMBarotropicViscousFluid does not have getSpeedOfSoundFromEnthalpy() implemented. Exiting." << std::endl;
    exit(0);

    return 0.0;
}

double FVMBarotropicViscousFluid::getPressureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Pt, double M){
    std::cerr << "ERROR: FVMBarotropicViscousFluid does not have getPressureFromStagnationProperties() implemented. Exiting." << std::endl;
    exit(0);

    return 0.0;
}

double FVMBarotropicViscousFluid::getTemperatureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Tt, double M){
    std::cerr << "ERROR: FVMBarotropicViscousFluid does not have getTemperatureFromStagnationProperties() implemented. Exiting." << std::endl;
    exit(0);

    return 0.0;
}

//mixture model functions
int FVMBarotropicViscousFluid::calculatePorosity(Job* job, FiniteVolumeDriver* driver){
    //mixture not implemented
    return 0;
} //return 1 if mixture problem, return 0 if FVM problem only

int FVMBarotropicViscousFluid::updateSolidPhaseVelocity(Job* job, FiniteVolumeDriver* driver){
    //mixture not implemented
    return 0;
} //return 1 if mixture problem, return 0 if FVM problem only

KinematicVector FVMBarotropicViscousFluid::getInterphaseDrag(Job* job, FiniteVolumeDriver* driver,
                                                             double rho,
                                                             const KinematicVector& v_f,
                                                             const KinematicVector& v_s,
                                                             double n){
    //mixture not implemented
    return KinematicVector(job->JOB_TYPE);
}

double FVMBarotropicViscousFluid::getInterphaseDragCoefficient(Job* job, FiniteVolumeDriver* driver,
                                                                        double rho,
                                                                        const KinematicVector& v_f,
                                                                        const KinematicVector& v_s,
                                                                        double n,
                                                                        int SPEC){
    return 0;
}

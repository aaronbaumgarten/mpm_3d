//
// Created by aaron on 4/28/20.
// fvm_ideal_gas.cpp
//

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
void FVMIdealGas::init(Job* job, FiniteVolumeDriver* driver){
    //check that simulation type is THERMAL
    if (driver->TYPE != FiniteVolumeDriver::THERMAL){
        std::cout << "WARNING! FVMIdealGas is a THERMAL model. Simulation TYPE "
                  << driver->TYPE << " != THERMAL." << std::endl;
    }

    if (fp64_props.size() < 4) {
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 4 properties defined (heat_capacity_ratio, R, eta, thermal_conductivity).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        heat_capacity_ratio = fp64_props[0];    //gamma
        R = fp64_props[1];                      //specific gas constant (J/(kg*K))
        eta = fp64_props[2];                    //shear viscosity (Pa*s)
        thermal_conductivity = fp64_props[3];   //k (W/(m*K))
    }
    std::cout << "FiniteVolumeMaterial properties:" << std::endl;
    std::cout << "    gamma = " << heat_capacity_ratio << std::endl;
    std::cout << "    R     = " << R << " J/(kg*K)" << std::endl;
    std::cout << "    eta   = " << eta << " Pa*s" << std::endl;
    std::cout << "    k     = " << thermal_conductivity << " W/(m*K)" << std::endl;
    std::cout << "FiniteVolumeMaterial initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//functions to calculate fluid fields
//get full stress tensor in fluid
MaterialTensor FVMIdealGas::getStress(Job* job,
                                                    FiniteVolumeDriver* driver,
                                                    const KinematicTensor& L,
                                                    double rho,
                                                    const KinematicVector& p,
                                                    double rhoE,
                                                    double n){
    //sigma = tau - p1
    return eta*(L + L.transpose() - 2.0/3.0*L.trace()*MaterialTensor::Identity())
           - (rhoE - 0.5*p.dot(p)/rho)*(heat_capacity_ratio - 1.0)*MaterialTensor::Identity();
}

//get shear components of stress tensor
MaterialTensor FVMIdealGas::getShearStress(Job* job,
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
double FVMIdealGas::getPressure(Job* job,
                                              FiniteVolumeDriver* driver,
                                              double rho,
                                              const KinematicVector& p,
                                              double rhoE,
                                              double n){
    //rhoE = p/(gamma - 1) + 0.5*rho*u^2
    return (rhoE - 0.5*p.dot(p)/rho)*(heat_capacity_ratio - 1.0);
}

//get temperature of fluid from state
double FVMIdealGas::getTemperature(Job* job,
                                                 FiniteVolumeDriver* driver,
                                                 double rho,
                                                 const KinematicVector& p,
                                                 double rhoE,
                                                 double n){
    //rhoE = rho*R*T/(gamma - 1) + 0.5*rho*u^2
    return (heat_capacity_ratio - 1.0)*(rhoE - 0.5*p.dot(p)/rho)/(rho*R);
}

//return speed of sound
double FVMIdealGas::getSpeedOfSound(Job* job,
                                                  FiniteVolumeDriver* driver,
                                                  double rho,
                                                  const KinematicVector& p,
                                                  double rhoE,
                                                  double n){
    //c^2 = sqrt((gamma - 1 + n)RT)/n =====> see pg. 154, nb #6
    double theta = getTemperature(job, driver, rho, p, rhoE, n);
    //return std::sqrt((heat_capacity_ratio - 1.0 + n)*R*theta)/n;
    return std::sqrt(heat_capacity_ratio*R*theta);
}

//loop over elements and fill in pressures
void FVMIdealGas::calculateElementPressures(Job* job, FiniteVolumeDriver* driver){
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        driver->fluid_body->P(e) = (driver->fluid_body->rhoE(e)
                                    - 0.5*driver->fluid_body->p[e].dot(driver->fluid_body->p[e])
                                      /driver->fluid_body->rho(e))*(heat_capacity_ratio - 1.0);
    }
    return;
}

//loop over elements and fill in shear stresses
void FVMIdealGas::calculateElementShearStresses(Job* job, FiniteVolumeDriver* driver){
    //get strain rates from grid reconstruction (not flux limited)
    KinematicTensorArray L = driver->fluid_body->L; //driver->fluid_grid->getVelocityGradients(job, driver);

    //fill in stresses accordingly
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        driver->fluid_body->tau[e] = eta*(L[e] + L[e].transpose() - 2.0/3.0*L[e].trace()*MaterialTensor::Identity());
    }
    return;
}

//loop over elements and fill in temperatures
void FVMIdealGas::calculateElementTemperatures(Job* job, FiniteVolumeDriver* driver){
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        driver->fluid_body->theta(e) = (heat_capacity_ratio - 1.0)
                                       *(driver->fluid_body->rhoE(e)
                                         - 0.5*driver->fluid_body->p[e].dot(driver->fluid_body->p[e])
                                           /driver->fluid_body->rho(e))
                                       /(driver->fluid_body->rho(e)*R);
    }
    return;
}

/*----------------------------------------------------------------------------*/
//fluid equations of state
double FVMIdealGas::getDensityFromPressureAndTemperature(Job* job,
                                                                       FiniteVolumeDriver* driver,
                                                                       double pressure,
                                                                       double theta,
                                                                       double n){
    //p = rho*R*T
    return pressure/(R*theta);
}

double FVMIdealGas::getInternalEnergyFromPressureAndTemperature(Job* job,
                                                                              FiniteVolumeDriver* driver,
                                                                              double pressure,
                                                                              double theta,
                                                                              double n){
    return pressure/(heat_capacity_ratio - 1.0);
}

double FVMIdealGas::getPressureFromDensityAndTemperature(Job* job,
                                                                       FiniteVolumeDriver* driver,
                                                                       double rho,
                                                                       double theta,
                                                                       double n){
    //P = rho*R*T
    return rho*R*theta;
}

KinematicVector FVMIdealGas::getHeatFlux(Job* job,
                                                       FiniteVolumeDriver* driver,
                                                       double rho,
                                                       double theta,
                                                       const KinematicVector& theta_x,
                                                       double n){
    //q = -k*grad(theta)
    return -thermal_conductivity*theta_x;
}

double FVMIdealGas::getSpeedOfSoundFromEnthalpy(Job *job,
                                                FiniteVolumeDriver *driver,
                                                double rho,
                                                const KinematicVector &p,
                                                double rhoH,
                                                double n) {
    //a^2 = (gamma - 1)*(H - 0.5*q^2)
    return std::sqrt((heat_capacity_ratio - 1.0) * (rhoH/rho - 0.5*p.dot(p)/(rho*rho)));
}

double FVMIdealGas::getPressureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Pt, double M){
    //(P^t/P) = (1 + M^2 * (gamma - 1)/2)^(gamma/gamma - 1)
    return Pt * std::pow(1.0 + M*M*(heat_capacity_ratio - 1)/2.0, -(heat_capacity_ratio)/(heat_capacity_ratio - 1.0));
}

double FVMIdealGas::getTemperatureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Tt, double M){
    //(T^t/T) = (1 + M^2 * (gamma - 1)/2)
    return Tt / (1.0 + M*M*(heat_capacity_ratio - 1)/2.0);
}


//mixture model functions
int FVMIdealGas::calculatePorosity(Job* job, FiniteVolumeDriver* driver){
    //mixture not implemented
    return 0;
} //return 1 if mixture problem, return 0 if FVM problem only

int FVMIdealGas::updateSolidPhaseVelocity(Job* job, FiniteVolumeDriver* driver){
    //mixture not implemented
    return 0;
} //return 1 if mixture problem, return 0 if FVM problem only

KinematicVector FVMIdealGas::getInterphaseDrag(Job* job, FiniteVolumeDriver* driver,
                                                             double rho,
                                                             const KinematicVector& v_f,
                                                             const KinematicVector& v_s,
                                                             double n){
    //mixture not implemented
    return KinematicVector(job->JOB_TYPE);
}

double FVMIdealGas::getInterphaseDragCoefficient(Job* job, FiniteVolumeDriver* driver,
                                                               double rho,
                                                               const KinematicVector& v_f,
                                                               const KinematicVector& v_s,
                                                               double n,
                                                               int SPEC){
    return 0;
}

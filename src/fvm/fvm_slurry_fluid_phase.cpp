//
// Created by aaron on 12/23/19.
// fvm_barotropic_viscous_fluid.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <job.hpp>
#include <objects/bodies/bodies.hpp>

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
void FVMSlurryFluidPhase::init(Job* job, FiniteVolumeDriver* driver){
    //check that simulation type is ISOTHERMAL
    if (driver->TYPE != FiniteVolumeDriver::ISOTHERMAL){
        std::cout << "WARNING! FVMSlurryFluidPhase is an ISOTHERMAL model. Simulation TYPE "
                  << driver->TYPE << " != ISOTHERMAL." << std::endl;
    }

    if (fp64_props.size() < 5 || (str_props.size() < 1 && int_props.size() < 1)){
        std::cout << fp64_props.size() <<", " << str_props.size() << ", " << int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 5 properties defined ({K, eta, rho_0, solid_rho, grain_diam}, {solid_body}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        kappa = fp64_props[0];
        eta = fp64_props[1];
        rho_0 = fp64_props[2];
        solid_rho = fp64_props[3];
        grain_diam = fp64_props[4];

        //set body ids by name
        if (str_props.size() == 1){
            for (int b = 0; b < job->bodies.size(); b++) {
                if (str_props[0].compare(job->bodies[b]->name) == 0){
                    solid_body_id = b;
                    break;
                }
            }
        }

        // or set body ids by int
        if (solid_body_id < 0){
            if (int_props.size() == 1) {
                solid_body_id = int_props[0];
            } else {
                std::cout << fp64_props.size() <<", " << str_props.size() << ", " << int_props.size() << "\n";
                fprintf(stderr,
                        "%s:%s: Need at least 5 properties defined ({K, eta, rho_0, solid_rho, grain_diam}, {solid_body}).\n",
                        __FILE__, __func__);
                exit(0);
            }
        }
    }

    std::cout << "FiniteVolumeMaterial properties (kappa = " << kappa << " Pa, eta = " << eta << " Pa*s, rho_0 = " << rho_0;
    std::cout << " kg/m^3, solid_rho = " << solid_rho << "kg/m^3, grain_diam = " << grain_diam << " m, solid_body_id = " << solid_body_id << ")." << std::endl;
    std::cout << "FiniteVolumeMaterial initialized." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//functions to calculate fluid fields
//get full stress tensor in fluid
MaterialTensor FVMSlurryFluidPhase::getStress(Job* job,
                                                    FiniteVolumeDriver* driver,
                                                    const KinematicTensor& L,
                                                    double rho,
                                                    const KinematicVector& p,
                                                    double rhoE,
                                                    double n){
    //sigma = tau - p1
    return eta*(1.0 + 5.0/2.0 * (1.0 - n))*(L + L.transpose() - 2.0/3.0*L.trace()*MaterialTensor::Identity())
           - kappa*std::log(rho/(n*rho_0))*MaterialTensor::Identity();
}

//get shear components of stress tensor
MaterialTensor FVMSlurryFluidPhase::getShearStress(Job* job,
                                                         FiniteVolumeDriver* driver,
                                                         const KinematicTensor& L,
                                                         double rho,
                                                         const KinematicVector& p,
                                                         double rhoE,
                                                         double n){
    //tau = 2*eta*(D - 1/3 trD I)
    return eta*(1.0 + 5.0/2.0 * (1.0 - n))*(L + L.transpose() - 2.0/3.0*L.trace()*MaterialTensor::Identity());
}

//get volumetric component of stress tensor
double FVMSlurryFluidPhase::getPressure(Job* job,
                                              FiniteVolumeDriver* driver,
                                              double rho,
                                              const KinematicVector& p,
                                              double rhoE,
                                              double n){
    //P = kappa*log(rho/rho_0)
    return kappa*std::log(rho/(n*rho_0));
}

//get temperature of fluid from state
double FVMSlurryFluidPhase::getTemperature(Job* job,
                                                 FiniteVolumeDriver* driver,
                                                 double rho,
                                                 const KinematicVector& p,
                                                 double rhoE,
                                                 double n){
    //isothermal model
    return 0.0;
}

//return speed of sound
double FVMSlurryFluidPhase::getSpeedOfSound(Job* job,
                                                  FiniteVolumeDriver* driver,
                                                  double rho,
                                                  const KinematicVector& p,
                                                  double rhoE,
                                                  double n){
    //c^2 = kappa/rho
    return std::sqrt(n*kappa/rho);
}

//loop over elements and fill in pressures
void FVMSlurryFluidPhase::calculateElementPressures(Job* job, FiniteVolumeDriver* driver){
    //map porosity field to elements
    //integrated porosity field
    Eigen::VectorXd n_e = driver->fluid_grid->M * driver->fluid_body->n;

    double n = 1.0;
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        //average porosity over element volume
        n = n_e(e) / driver->fluid_grid->getElementVolume(e);

        //fill in centroid pressure
        driver->fluid_body->P(e) = kappa*std::log(driver->fluid_body->rho(e)/(n*rho_0));
    }
    return;
}

//loop over elements and fill in shear stresses
void FVMSlurryFluidPhase::calculateElementShearStresses(Job* job, FiniteVolumeDriver* driver){
    //get strain rates from grid reconstruction (not flux limited)
    KinematicTensorArray L = driver->fluid_body->L; //driver->fluid_grid->getVelocityGradients(job, driver);

    //map porosity field to elements
    //integrated porosity field
    Eigen::VectorXd n_e = driver->fluid_grid->M * driver->fluid_body->n;

    double n = 1.0;
    //fill in stresses accordingly
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        //average porosity over element volume
        n = n_e(e) / driver->fluid_grid->getElementVolume(e);

        //fill in centroid shear stress
        driver->fluid_body->tau[e] = eta*(1.0 + 5.0/2.0 * (1.0 - n))*(L[e] + L[e].transpose() - 2.0/3.0*L[e].trace()*MaterialTensor::Identity());
    }
    return;
}

//loop over elements and fill in temperatures
void FVMSlurryFluidPhase::calculateElementTemperatures(Job* job, FiniteVolumeDriver* driver){
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        driver->fluid_body->theta(e) = 0.0;
    }
    return;
}

/*----------------------------------------------------------------------------*/
//fluid equations of state
double FVMSlurryFluidPhase::getDensityFromPressureAndTemperature(Job* job,
                                                                       FiniteVolumeDriver* driver,
                                                                       double pressure,
                                                                       double theta,
                                                                       double n){
    return n*rho_0*std::exp(pressure/kappa);
}

double FVMSlurryFluidPhase::getInternalEnergyFromPressureAndTemperature(Job* job,
                                                                              FiniteVolumeDriver* driver,
                                                                              double pressure,
                                                                              double theta,
                                                                              double n){
    return 0.0;
}

double FVMSlurryFluidPhase::getPressureFromDensityAndTemperature(Job* job,
                                                                       FiniteVolumeDriver* driver,
                                                                       double rho,
                                                                       double theta,
                                                                       double n){
    //P = kappa*log(rho/rho_0)
    return kappa*std::log(rho/(n*rho_0));
}

KinematicVector FVMSlurryFluidPhase::getHeatFlux(Job* job,
                                                       FiniteVolumeDriver* driver,
                                                       double rho,
                                                       double theta,
                                                       const KinematicVector& theta_x,
                                                       double n){
    //no heat flux for now...
    return KinematicVector(job->JOB_TYPE);
}


double FVMSlurryFluidPhase::getSpeedOfSoundFromEnthalpy(Job *job, FiniteVolumeDriver *driver, double rho, const KinematicVector &p, double rhoH,
                                                        double n) {
    std::cerr << "ERROR: FVMSlurryFluidPhase does not have getSpeedOfSoundFromEnthalpy() implemented. Exiting." << std::endl;
    exit(0);

    return 0.0;
}

double FVMSlurryFluidPhase::getPressureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Pt, double M){
    //std::cerr << "ERROR: FVMSlurryFluidPhase does not have getPressureFromStagnationProperties() implemented. Exiting." << std::endl;
    //exit(0);

    //see pg 123 nb #7
    return Pt + kappa/2.0 * std::log(1.0 - M*M/2.0);
}

double FVMSlurryFluidPhase::getTemperatureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Tt, double M){
    //std::cerr << "ERROR: FVMSlurryFluidPhase does not have getTemperatureFromStagnationProperties() implemented. Exiting." << std::endl;
    //exit(0);

    return 0.0;
}


//mixture model functions
int FVMSlurryFluidPhase::calculatePorosity(Job* job, FiniteVolumeDriver* driver){
    //generate approximate porosity field from MPM body
    double m1;
    for (int i=0;i<driver->fluid_body->n.rows();i++){

        m1 = job->bodies[solid_body_id]->nodes->m[i];
        if (m1 > 0){
            driver->fluid_body->n(i) = 1.0 - (m1 / (job->grid->nodeVolume(job,i)*solid_rho));

            if (driver->fluid_body->n(i) < 1e-2){
                //if porosity drops below 0.2, non-physical
                //std::cout << "WARNING: Porosity over-estimate!" << std::endl;
                driver->fluid_body->n(i) = 1e-2;
            }
        } else {
            driver->fluid_body->n(i) = 1.0;
        }
    }
    return 1;
} //return 1 if mixture problem, return 0 if FVM problem only

int FVMSlurryFluidPhase::updateSolidPhaseVelocity(Job* job, FiniteVolumeDriver* driver){
    //copy solid phase velocity field
    driver->fluid_body->v_s = job->bodies[solid_body_id]->nodes->x_t;

    return 1;
} //return 1 if mixture problem, return 0 if FVM problem only

KinematicVector FVMSlurryFluidPhase::getInterphaseDrag(Job* job, FiniteVolumeDriver* driver,
                                                             double rho,
                                                             const KinematicVector& v_f,
                                                             const KinematicVector& v_s,
                                                             double n){
    //permeability
    double Re = (v_s - v_f).norm() * rho * grain_diam / eta;
    double C = 0;
    //Beetstra
    if (n < 1e-10 || ((1-n) < 1e-10) || eta < 1e-10){
        C = 0;
    } else if (Re < 1e-10){
        C = 18.0 * (1 - n) * eta / (grain_diam * grain_diam) *
            (10.0 * (1 - n) / n + n * n * n * (1.0 + 1.5 * std::sqrt(1 - n)));
    } else {
        C = 18.0 * (1 - n) * eta / (grain_diam * grain_diam) *
            (10.0 * (1 - n) / n + n * n * n * (1.0 + 1.5 * std::sqrt(1 - n)) +
             0.413 * Re / (24.0 * n) *
             (1 / n + 3 * n * (1 - n) + 8.4 * std::pow(Re, -0.343)) /
             (1 + std::pow(10.0, 3 * (1 - n)) * std::pow(Re, -0.5 * (1 + 4 * (1 - n)))));
    }

    /*
    if (C > 1e-10) {
        std::cout << "C: " << C;
    }

    //approximate limit for C in explicit simulations
    if (n >= 1e-10 && ((1-n) >= 1e-10)){
        //std::cout << C << ", " << C / (1 + job->dt*C*(1.0/rho + n/(rho * (1.0 - n)))) << std::endl;
        C /= (1 + job->dt*C*(1.0/rho + n/(rho * (1.0 - n))));
    }

    if (C > 1e-10) {
        std::cout << " ?= " << C << std::endl;
    }
     */

    return C * (v_s - v_f);
}

double FVMSlurryFluidPhase::getInterphaseDragCoefficient(Job* job, FiniteVolumeDriver* driver,
                                                         double rho,
                                                         const KinematicVector& v_f,
                                                         const KinematicVector& v_s,
                                                         double n,
                                                         int SPEC){
    //permeability
    double Re = (v_s - v_f).norm() * rho * grain_diam / eta;
    double C = 0;
    //Beetstra
    if (n < 1e-10 || ((1-n) < 1e-10) || eta < 1e-10){
        return 0;
    } else if (Re < 1e-10){
        C = 18.0 * (1 - n) * eta / (grain_diam * grain_diam) *
            (10.0 * (1 - n) / n + n * n * n * (1.0 + 1.5 * std::sqrt(1 - n)));
    } else {
        C = 18.0 * (1 - n) * eta / (grain_diam * grain_diam) *
            (10.0 * (1 - n) / n + n * n * n * (1.0 + 1.5 * std::sqrt(1 - n)) +
             0.413 * Re / (24.0 * n) *
             (1 / n + 3 * n * (1 - n) + 8.4 * std::pow(Re, -0.343)) /
             (1 + std::pow(10.0, 3 * (1 - n)) * std::pow(Re, -0.5 * (1 + 4 * (1 - n)))));
    }

    if (SPEC == REGULAR_DRAG) {
        return C;
    } else {
        double rho_s = solid_rho * (1.0 - n);
        return C*rho*rho_s / (rho*rho_s + job->dt * C * (rho + rho_s));
    }
}
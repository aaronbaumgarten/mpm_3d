//
// Created by aaron on 4/28/20.
// fvm_ideal_gas.cpp
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
void FVMSlurryGasPhase::init(Job* job, FiniteVolumeDriver* driver){
    //check that simulation type is THERMAL
    if (driver->TYPE != FiniteVolumeDriver::THERMAL){
        std::cout << "WARNING! FVMSlurryGasPhase is a THERMAL model. Simulation TYPE "
                  << driver->TYPE << " != THERMAL." << std::endl;
    }

    if (fp64_props.size() < 6 || (str_props.size() < 1 && int_props.size() < 1)) {
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 7 properties defined ({heat_capacity_ratio, R, eta, thermal_conductivity, solid_rho, grain_diam}, {solid_body_id}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        heat_capacity_ratio = fp64_props[0];    //gamma
        R = fp64_props[1];                      //specific gas constant (J/(kg*K))
        eta = fp64_props[2];                    //shear viscosity (Pa*s)
        thermal_conductivity = fp64_props[3];   //k (W/(m*K))

        solid_rho = fp64_props[4]; //grain density (kg/m^3)
        grain_diam = fp64_props[5]; //grain diameter (m)

        //set body ids by name
        if (str_props.size() >= 1){
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
                        "%s:%s: Need at least 7 properties defined ({heat_capacity_ratio, R, eta, thermal_conductivity, solid_rho, grain_diam}, {solid_body_id}).\n",
                        __FILE__, __func__);
                exit(0);
            }
        }
    }

    //loop over str-props and assign relevant flags
    int len;
    std::vector<std::string> options = {"USE_EDDY_VISCOSITY"};
    for (int i=0; i<str_props.size(); i++){
        switch (Parser::findStringID(options, str_props[i])){
            case 0:
                //USE_EDDY_VISCOSITY
                USE_EDDY_VISCOSITY = true;
                len = fp64_props.size();
                a_T = fp64_props[len - 1];
                h_T = fp64_props[len - 2];
                std::cout << "FVMSlurryGasPhase using eddy viscosity. a = " << a_T << ", h = " << h_T << std::endl;
                break;
            default:
                //do nothing
                break;
        }
    }

    std::cout << "FiniteVolumeMaterial properties:" << std::endl;
    std::cout << "    gamma = " << heat_capacity_ratio << std::endl;
    std::cout << "    R     = " << R << " J/(kg*K)" << std::endl;
    std::cout << "    eta   = " << eta << " Pa*s" << std::endl;
    std::cout << "    k     = " << thermal_conductivity << " W/(m*K)" << std::endl;
    std::cout << "    rho_s = " << solid_rho << " (kg/m^3)" << std::endl;
    std::cout << "    d     = " << grain_diam << " (m)" << std::endl;
    std::cout << "    solid = \"" << job->bodies[solid_body_id]->name << "\"" << std::endl;
    std::cout << "FiniteVolumeMaterial initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//functions to calculate fluid fields
//get full stress tensor in fluid
MaterialTensor FVMSlurryGasPhase::getStress(Job* job,
                                                    FiniteVolumeDriver* driver,
                                                    const KinematicTensor& L,
                                                    double rho,
                                                    const KinematicVector& p,
                                                    double rhoE,
                                                    double n){
    //sigma = tau - p1
    if (!USE_EDDY_VISCOSITY) {
        return eta * (1.0 + 5.0 / 2.0 * (1.0 - n)) *
               (L + L.transpose() - 2.0 / 3.0 * L.trace() * MaterialTensor::Identity())
               - (rhoE - 0.5 * p.dot(p) / rho) * (heat_capacity_ratio - 1.0) / (n) * MaterialTensor::Identity();
    } else {
        double mu_T = a_T * a_T * h_T * h_T * rho/n * (L + L.transpose() - 2.0 / 3.0 * L.trace() * MaterialTensor::Identity()).norm()  / std::sqrt(2);
        return (eta * (1.0 + 5.0 / 2.0 * (1.0 - n)) + mu_T) *
               (L + L.transpose() - 2.0 / 3.0 * L.trace() * MaterialTensor::Identity())
               - (rhoE - 0.5 * p.dot(p) / rho) * (heat_capacity_ratio - 1.0) / (n) * MaterialTensor::Identity();
    }
}

//get shear components of stress tensor
MaterialTensor FVMSlurryGasPhase::getShearStress(Job* job,
                                                         FiniteVolumeDriver* driver,
                                                         const KinematicTensor& L,
                                                         double rho,
                                                         const KinematicVector& p,
                                                         double rhoE,
                                                         double n){
    //tau = 2*eta*(D - 1/3 trD I)
    if (!USE_EDDY_VISCOSITY) {
        return eta * (1.0 + 5.0 / 2.0 * (1.0 - n)) *
               (L + L.transpose() - 2.0 / 3.0 * L.trace() * MaterialTensor::Identity());
    } else {
        double mu_T = a_T * a_T * h_T * h_T * rho/n
                      * (L + L.transpose() - 2.0 / 3.0 * L.trace() * MaterialTensor::Identity()).norm() / std::sqrt(2);
        return (eta * (1.0 + 5.0 / 2.0 * (1.0 - n)) + mu_T) *
               (L + L.transpose() - 2.0 / 3.0 * L.trace() * MaterialTensor::Identity());
    }
}

//get volumetric component of stress tensor
double FVMSlurryGasPhase::getPressure(Job* job,
                                              FiniteVolumeDriver* driver,
                                              double rho,
                                              const KinematicVector& p,
                                              double rhoE,
                                              double n){
    //rhoE = n*p/(gamma - 1) + 0.5*rho*u^2
    return (rhoE - 0.5*p.dot(p)/rho)*(heat_capacity_ratio - 1.0)/n;
}

//get temperature of fluid from state
double FVMSlurryGasPhase::getTemperature(Job* job,
                                                 FiniteVolumeDriver* driver,
                                                 double rho,
                                                 const KinematicVector& p,
                                                 double rhoE,
                                                 double n){
    //rhoE = n*rho*R*T/(gamma - 1) + 0.5*rho*u^2
    return (heat_capacity_ratio - 1.0)*(rhoE - 0.5*p.dot(p)/rho)/(rho*R);
}

//return speed of sound
double FVMSlurryGasPhase::getSpeedOfSound(Job* job,
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
void FVMSlurryGasPhase::calculateElementPressures(Job* job, FiniteVolumeDriver* driver){
    //integrated porosity field
    Eigen::VectorXd n_e = driver->fluid_grid->M * driver->fluid_body->n;

    double n = 1.0;
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        //average porosity over element volume
        n = n_e(e) / driver->fluid_grid->getElementVolume(e);
        driver->fluid_body->P(e) = (driver->fluid_body->rhoE(e)
                                    - 0.5*driver->fluid_body->p[e].dot(driver->fluid_body->p[e])
                                      /driver->fluid_body->rho(e))*(heat_capacity_ratio - 1.0)/n;
    }
    return;
}

//loop over elements and fill in shear stresses
void FVMSlurryGasPhase::calculateElementShearStresses(Job* job, FiniteVolumeDriver* driver){
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
        if (!USE_EDDY_VISCOSITY) {
            driver->fluid_body->tau[e] = eta * (1.0 + 5.0 / 2.0 * (1.0 - n)) * (L[e] + L[e].transpose() -
                                                                                2.0 / 3.0 * L[e].trace() *
                                                                                MaterialTensor::Identity());
        } else {
            double mu_T = a_T * a_T * h_T * h_T * driver->fluid_body->rho(e)/n * (L[e] + L[e].transpose()
                                                                                  - 2.0 / 3.0 * L[e].trace()
                                                                                    * MaterialTensor::Identity()).norm() / std::sqrt(2);
            driver->fluid_body->tau[e] = (mu_T + eta * (1.0 + 5.0 / 2.0 * (1.0 - n))) * (L[e] + L[e].transpose() -
                                                                                         2.0 / 3.0 * L[e].trace() *
                                                                                         MaterialTensor::Identity());
        }
    }
    return;
}

//loop over elements and fill in temperatures
void FVMSlurryGasPhase::calculateElementTemperatures(Job* job, FiniteVolumeDriver* driver){
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
double FVMSlurryGasPhase::getDensityFromPressureAndTemperature(Job* job,
                                                                       FiniteVolumeDriver* driver,
                                                                       double pressure,
                                                                       double theta,
                                                                       double n){
    //p = rho*R*T
    return n*pressure/(R*theta);
}

double FVMSlurryGasPhase::getInternalEnergyFromPressureAndTemperature(Job* job,
                                                                              FiniteVolumeDriver* driver,
                                                                              double pressure,
                                                                              double theta,
                                                                              double n){
    return n*pressure/(heat_capacity_ratio - 1.0);
}

double FVMSlurryGasPhase::getPressureFromDensityAndTemperature(Job* job,
                                                                       FiniteVolumeDriver* driver,
                                                                       double rho,
                                                                       double theta,
                                                                       double n){
    //P = rho*R*T
    return rho*R*theta/n;
}

KinematicVector FVMSlurryGasPhase::getHeatFlux(Job* job,
                                                       FiniteVolumeDriver* driver,
                                                       double rho,
                                                       double theta,
                                                       const KinematicVector& theta_x,
                                                       double n){
    //q = -k*grad(theta)
    return -thermal_conductivity*theta_x*n;
}

double FVMSlurryGasPhase::getSpeedOfSoundFromEnthalpy(Job *job,
                                                FiniteVolumeDriver *driver,
                                                double rho,
                                                const KinematicVector &p,
                                                double rhoH,
                                                double n) {
    //a^2 = (gamma - 1)*(H - 0.5*q^2)
    //return std::sqrt((heat_capacity_ratio - 1.0)/n * (rhoH/rho - 0.5*p.dot(p)/(rho*rho)));
    return std::sqrt((heat_capacity_ratio - 1.0) * (rhoH/rho - 0.5*p.dot(p)/(rho*rho)));
}

double FVMSlurryGasPhase::getPressureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Pt, double M){
    //(P^t/P) = (1 + M^2 * (gamma - 1)/2)^(gamma/gamma - 1)
    return Pt * std::pow(1.0 + M*M*(heat_capacity_ratio - 1)/2.0, -(heat_capacity_ratio)/(heat_capacity_ratio - 1.0));
}

double FVMSlurryGasPhase::getTemperatureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Tt, double M){
    //(T^t/T) = (1 + M^2 * (gamma - 1)/2)
    return Tt / (1.0 + M*M*(heat_capacity_ratio - 1)/2.0);
}

//mixture model functions
int FVMSlurryGasPhase::calculatePorosity(Job* job, FiniteVolumeDriver* driver){
    //generate approximate porosity field from MPM body
    double m1;
    for (int i=0;i<driver->fluid_body->n.rows();i++){

        m1 = job->bodies[solid_body_id]->nodes->m(i);
        if (m1 > 0){
            driver->fluid_body->n(i) = 1.0 - (m1 / (job->grid->nodeVolume(job,i)*solid_rho));

            if (driver->fluid_body->n(i) < 0.2){
                //if porosity drops below 0.2, non-physical
                //std::cout << "WARNING: Porosity over-estimate!" << std::endl;
                driver->fluid_body->n(i) = 0.2;
            }
        } else {
            driver->fluid_body->n(i) = 1.0;
        }
    }
    return 1;
} //return 1 if mixture problem, return 0 if FVM problem only

int FVMSlurryGasPhase::updateSolidPhaseVelocity(Job* job, FiniteVolumeDriver* driver){
    //copy solid phase velocity field
    driver->fluid_body->v_s = job->bodies[solid_body_id]->nodes->x_t;

    return 1;
} //return 1 if mixture problem, return 0 if FVM problem only

KinematicVector FVMSlurryGasPhase::getInterphaseDrag(Job* job, FiniteVolumeDriver* driver,
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

    return C * (v_s - v_f);
}

double FVMSlurryGasPhase::getInterphaseDragCoefficient(Job* job, FiniteVolumeDriver* driver,
                                                               double rho,
                                                               const KinematicVector& v_f,
                                                               const KinematicVector& v_s,
                                                               double n,
                                                               int SPEC){
    //if eta small, C unstable
    if (eta < 1e-10){
        return 0.0;
    }

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

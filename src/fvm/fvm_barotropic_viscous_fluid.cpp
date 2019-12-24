//
// Created by aaron on 12/23/19.
// fvm_barotropic_viscous_fluid.cpp
//

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
#include "fvm_materials.hpp"

/*----------------------------------------------------------------------------*/
//initialize from job and driver
void FVMBarotropicViscousFluid::init(Job* job, FiniteVolumeDriver* driver){
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
//get full stress tensor in fluid at given spatial position
MaterialTensor FVMBarotropicViscousFluid::getStress(Job* job,
                                                    FiniteVolumeDriver* driver,
                                                    KinematicVector x,
                                                    KinematicTensor L,
                                                    double rho,
                                                    double theta){
    //sigma = tau - p1
    return eta*(L + L.transpose() - 2.0/3.0*L.trace()*MaterialTensor::Identity())
           - kappa*std::log(rho/rho_0)*MaterialTensor::Identity();
}

//get shear component of fluid stress
MaterialTensor FVMBarotropicViscousFluid::getShearStress(Job* job,
                                                         FiniteVolumeDriver* driver,
                                                         KinematicVector x,
                                                         KinematicTensor L,
                                                         double rho,
                                                         double theta){
    //tau = 2*eta*(D - 1/3 trD I)
    return eta*(L + L.transpose() - 2.0/3.0*L.trace()*MaterialTensor::Identity());
}

//get volumetric component of fluid stress
double FVMBarotropicViscousFluid::getPressure(Job* job,
                                              FiniteVolumeDriver* driver,
                                              KinematicVector x,
                                              double rho,
                                              double theta){
    //P = kappa*log(rho/rho_0)
    return kappa*std::log(rho/rho_0);
}

//return speed of sound (dP/drho)
double FVMBarotropicViscousFluid::getSpeedOfSound(Job* job,
                                                  FiniteVolumeDriver* driver,
                                                  KinematicVector x,
                                                  double rho,
                                                  double theta){
    //c^2 = kappa/rho
    return std::sqrt(kappa/rho);
}

//loop over elements and calculate pressure at centroids
void FVMBarotropicViscousFluid::calculateElementPressures(Job* job,
                                                          FiniteVolumeDriver* driver){
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        driver->fluid_body->P(e) = kappa*std::log(driver->fluid_body->rho(e)/rho_0);
    }
    return;
}

//loop over elements and calculate shear stresses
void FVMBarotropicViscousFluid::calculateElementShearStresses(Job* job,
                                                              FiniteVolumeDriver* driver){
    //get strain rates from grid reconstruction (not flux limited)
    KinematicTensorArray L = driver->fluid_grid->getVelocityGradients(job, driver);

    //fill in stresses accordingly
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        driver->fluid_body->tau[e] = eta*(L[e] + L[e].transpose() - 2.0/3.0*L[e].trace()*MaterialTensor::Identity());
    }
    return;
}

//advect material properties through domain
void FVMBarotropicViscousFluid::solveMaterialEquations(Job* job, FiniteVolumeDriver* driver){
    //do nothing
    return;
}
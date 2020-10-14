//
// Created by aaron on 10/14/20.
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
#include "fvm_artificial_viscosity.hpp"

double ArtificialViscosityCalculator::getS(Job* job,
            FiniteVolumeDriver* driver,
            int f){
    //function to calculate shock capturing variable s
    return 1.0;
}

ArtificialViscosityCalculator::fluxVector ArtificialViscosityCalculator::getArtificialViscosityFlux(Job* job,
                                  FiniteVolumeDriver* driver,
                                  int f){
    //function to add artificial viscous flux to ALL conserved fields

    //get shock capturing variable
    double s = ArtificialViscosityCalculator::getS(job, driver, f);

    //get elements and face normal
    int e_minus = driver->fluid_grid->getOrientedElementsByFace(f)[0];
    int e_plus = driver->fluid_grid->getOrientedElementsByFace(f)[1];
    KinematicVector normal = driver->fluid_grid->getFaceNormal(job, f);

    //calculate max wave speed
    double rho_plus = driver->fluid_body->rho(e_plus);
    double rho_minus = driver->fluid_body->rho(e_minus);
    KinematicVector u_plus = driver->fluid_body->p(e_plus)/rho_plus;
    KinematicVector u_minus = driver->fluid_body->p(e_minus)/rho_minus;
    double E_plus = driver->fluid_body->rhoE(e_plus)/rho_plus;
    double E_minus = driver->fluid_body->rhoE(e_minus)/rho_minus;

    //scale rho_plus and rho_minus
    rho_plus /= driver->fluid_body->n_e(e_plus);
    rho_minus /= driver->fluid_body->n_e(e_minus);
    double n_bar = (driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e(e_minus)) / 2.0;

    double P_plus = driver->fluid_material->getPressure(job, driver,
                                                 n_bar*rho_plus,
                                                 n_bar*rho_plus*u_plus,
                                                 n_bar*rho_plus*E_plus,
                                                 n_bar); //n_q(q_list[q]));
    double H_plus = E_plus + P_plus/rho_plus;

    double P_minus = driver->fluid_material->getPressure(job, driver,
                                                        n_bar*rho_minus,
                                                        n_bar*rho_minus*u_minus,
                                                        n_bar*rho_minus*E_minus,
                                                        n_bar); //n_q(q_list[q]));
    double H_minus = E_minus + P_minus/rho_minus;

    //approximate Roe advective rate
    double w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus))/2.0;
    double rho_bar = w_1_bar*w_1_bar;
    KinematicVector u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) / (2.0 * w_1_bar);
    double u_bar_squared = u_bar.dot(u_bar);
    double H_bar = (std::sqrt(rho_plus) * H_plus + std::sqrt(rho_minus) * H_minus) / (2.0 * w_1_bar);

    double c = 0;
    if (driver->TYPE == FiniteVolumeDriver::THERMAL) {
        c = driver->fluid_material->getSpeedOfSoundFromEnthalpy(job,
                                                                driver,
                                                                n_bar * rho_bar,
                                                                n_bar * rho_bar * u_bar,
                                                                n_bar * rho_bar * H_bar,
                                                                n_bar);
    } else if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL){
        c = driver->fluid_material->getSpeedOfSound(job, driver,
                                                           n_bar * rho_bar,
                                                           n_bar * rho_bar * u_bar,
                                                           0.0,
                                                           n_bar);
    } else {
        std::cerr << "ERROR: ArtificialViscosityCalculator not implemented for driver types of " << driver->TYPE;
        std::cerr << "! Exiting." << std::endl;
        exit(0);
    }

    //find maximum lambda
    double lambda_max = std::abs(u_bar.dot(normal)) + c;

    //apply artificial viscosity model
    //\epsilon * h * div(\nabla U)
    double V_plus = driver->fluid_grid->getElementVolume(e_plus);
    double V_minus = driver->fluid_grid->getElementVolume(e_minus);
    double A = driver->fluid_grid->getFaceArea(f);

    //return flux along n direction
    fluxVector result;
    result.rho = -lambda_max * (rho_plus - rho_minus) * A * n_bar;
    result.p = -lambda_max * (rho_plus*u_plus - rho_minus*u_minus) * A * n_bar;
    result.rhoE = -lambda_max * (rho_plus*E_plus - rho_minus*E_minus) * A * n_bar;

    return result;
}


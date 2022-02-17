//
// Created by aaron on 10/14/20.
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
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
#include "fvm_serializers.hpp"

double ArtificialViscosityCalculator::getAVCoeff(Job* job,
                                                 FiniteVolumeDriver* driver,
                                                 int f){
    //function to calculate shock capturing variable s
    //return 1.0;

    bool USE_CONSTANT_S = false;
    bool USE_MCCOR_COL = false;
    bool USE_COUCHMAN = false;
    bool USE_MCCOR_COL_3 = false;
    bool USE_CUSTOM = true;
    bool USE_ELEMENT_S = false;

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

    if (USE_MCCOR_COL) { //USE McCorquodale and Collelo shock capturing method
        double c_plus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                n_bar * rho_plus,
                                                                n_bar * rho_plus * u_plus,
                                                                rho_plus * E_plus,
                                                                n_bar);

        double c_minus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                 n_bar * rho_minus,
                                                                 n_bar * rho_minus * u_minus,
                                                                 rho_minus * E_minus,
                                                                 n_bar);

        double h_lambda = (u_plus - u_minus).dot(normal); //?
        double alpha = 3.0;//0.3;
        double beta = 0.3;
        double c_min = std::min(c_plus, c_minus);

        if (h_lambda < 0){
            return alpha * std::abs(h_lambda) * std::min(1.0, h_lambda*h_lambda/(beta * c_min * c_min));
        } else {
            return 0.0;
        }

    } else if (USE_CONSTANT_S) { //USE uniform artificial viscosity

        double P_plus = driver->fluid_material->getPressure(job, driver,
                                                            n_bar * rho_plus,
                                                            n_bar * rho_plus * u_plus,
                                                            n_bar * rho_plus * E_plus,
                                                            n_bar); //n_q(q_list[q]));
        double H_plus = E_plus + P_plus / rho_plus;

        double P_minus = driver->fluid_material->getPressure(job, driver,
                                                             n_bar * rho_minus,
                                                             n_bar * rho_minus * u_minus,
                                                             n_bar * rho_minus * E_minus,
                                                             n_bar); //n_q(q_list[q]));
        double H_minus = E_minus + P_minus / rho_minus;

        //approximate Roe advective rate
        double w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
        double rho_bar = w_1_bar * w_1_bar;
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
        } else if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
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
        double lambda_max = (std::abs(u_bar.dot(normal)) + c);

        return lambda_max;
    } else if (USE_COUCHMAN) {
        double c_plus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                n_bar * rho_plus,
                                                                n_bar * rho_plus * u_plus,
                                                                rho_plus * E_plus,
                                                                n_bar);

        double c_minus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                 n_bar * rho_minus,
                                                                 n_bar * rho_minus * u_minus,
                                                                 rho_minus * E_minus,
                                                                 n_bar);

        KinematicVector h = driver->fluid_grid->getElementCentroid(job,e_plus) -
                            driver->fluid_grid->getElementCentroid(job,e_minus);
        double h_lambda = (u_plus + u_minus).dot(h) / (2.0*h.norm()) * (u_plus.norm()/c_plus - u_minus.norm()/c_minus); //?
        double c_max = std::max(c_plus + u_plus.norm(), c_minus + u_minus.norm());

        if (h_lambda < 0){
            return c_max * std::min(std::abs(h_lambda)/c_max, 1.0);
        } else {
            return 0.0;
        }
    } else if (USE_MCCOR_COL_3) {
        double P_plus = driver->fluid_material->getPressure(job, driver,
                                                            n_bar * rho_plus,
                                                            n_bar * rho_plus * u_plus,
                                                            n_bar * rho_plus * E_plus,
                                                            n_bar); //n_q(q_list[q]));
        double H_plus = E_plus + P_plus / rho_plus;

        double P_minus = driver->fluid_material->getPressure(job, driver,
                                                             n_bar * rho_minus,
                                                             n_bar * rho_minus * u_minus,
                                                             n_bar * rho_minus * E_minus,
                                                             n_bar); //n_q(q_list[q]));
        double H_minus = E_minus + P_minus / rho_minus;

        //approximate Roe advective rate
        double w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
        double rho_bar = w_1_bar * w_1_bar;
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
        } else if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
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
        double lambda_max = (std::abs(u_bar.dot(normal)) + c);
        double beta = 0.01;
        double h_lambda = (u_plus - u_minus).dot(normal);

        if (h_lambda < 0){
            return lambda_max * std::min(1.0, h_lambda*h_lambda/(beta * c * c));
        } else {
            return 0.0;
        }
    } else if (USE_CUSTOM) {
        /*
        double P_plus = driver->fluid_material->getPressure(job, driver,
                                                            n_bar * rho_plus,
                                                            n_bar * rho_plus * u_plus,
                                                            n_bar * rho_plus * E_plus,
                                                            n_bar); //n_q(q_list[q]));
        double H_plus = E_plus + P_plus / rho_plus;

        double P_minus = driver->fluid_material->getPressure(job, driver,
                                                             n_bar * rho_minus,
                                                             n_bar * rho_minus * u_minus,
                                                             n_bar * rho_minus * E_minus,
                                                             n_bar); //n_q(q_list[q]));
        double H_minus = E_minus + P_minus / rho_minus;

        //approximate Roe advective rate
        double w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
        double rho_bar = w_1_bar * w_1_bar;
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
        } else if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
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
         */

        double c_plus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                n_bar * rho_plus,
                                                                n_bar * rho_plus * u_plus,
                                                                rho_plus * E_plus,
                                                                n_bar);

        double c_minus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                 n_bar * rho_minus,
                                                                 n_bar * rho_minus * u_minus,
                                                                 rho_minus * E_minus,
                                                                 n_bar);

        //find maximum lambda
        double lambda_max = std::max(std::abs(u_plus.dot(normal)) + c_plus,
                                     std::abs(u_minus.dot(normal)) + c_minus); //(std::abs(u_bar.dot(normal)) + c);
        double beta = 0.3;
        double h_lambda = (u_plus - u_minus).dot(normal);

        if (h_lambda < 0){
            double c_min = std::min(c_plus, c_minus);
            double epsilon_v = -h_lambda + std::abs(c_plus - c_minus);
            return lambda_max * std::min(1.0, epsilon_v*epsilon_v/(beta * beta * c_min * c_min));
            //return lambda_max * std::min(1.0, epsilon_v / (beta * c_min));
        } else {
            return 0.0;
        }
    } else if (USE_ELEMENT_S) {
        double c_plus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                n_bar * rho_plus,
                                                                n_bar * rho_plus * u_plus,
                                                                rho_plus * E_plus,
                                                                n_bar);

        double c_minus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                 n_bar * rho_minus,
                                                                 n_bar * rho_minus * u_minus,
                                                                 rho_minus * E_minus,
                                                                 n_bar);

        //find maximum lambda
        double lambda_max = std::max(std::abs(u_plus.dot(normal)) + c_plus,
                                     std::abs(u_minus.dot(normal)) + c_minus); //(std::abs(u_bar.dot(normal)) + c);
        double beta = 0.3;
        KinematicVector h = driver->fluid_grid->getElementCentroid(job,e_plus) -
                            driver->fluid_grid->getElementCentroid(job,e_minus);
        double h_lambda = (u_plus - u_minus).dot(normal);
        h_lambda = h.norm() * 0.5 * ((driver->fluid_body->L[e_plus] + driver->fluid_body->L[e_minus])*normal).dot(normal); //

        if (h_lambda < 0){
            double c_min = std::min(c_plus, c_minus);
            double epsilon_v = -h_lambda + std::abs(c_plus - c_minus);
            return lambda_max * std::min(1.0, epsilon_v*epsilon_v/(beta * beta * c_min * c_min));
            //return lambda_max * std::min(1.0, epsilon_v / (beta * c_min));
        } else {
            return 0.0;
        }
    }

    return 0.0;
}

ArtificialViscosityCalculator::fluxVector ArtificialViscosityCalculator::getArtificialViscosityFlux(Job* job,
                                  FiniteVolumeDriver* driver,
                                  int f){
    //function to add artificial viscous flux to ALL conserved fields

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

    //apply artificial viscosity model
    //\epsilon * div(\nabla U)
    double V_plus = driver->fluid_grid->getElementVolume(e_plus);
    double V_minus = driver->fluid_grid->getElementVolume(e_minus);
    double A = driver->fluid_grid->getFaceArea(f);

    //get shock capturing variable
    double epsilon = ArtificialViscosityCalculator::getAVCoeff(job, driver, f);

    //return flux along n direction
    fluxVector result;
    result.rho = -epsilon * (rho_plus - rho_minus) * A * n_bar;
    result.p = -epsilon * (rho_plus*u_plus - rho_minus*u_minus) * A * n_bar;
    result.rhoE = -epsilon * (rho_plus*E_plus - rho_minus*E_minus) * A * n_bar;

    return result;
}


void ArtificialViscosityCalculator::writeFrame(Job* job, FiniteVolumeDriver* driver){
    //function to write average shock sensor value in each cell
    Eigen::VectorXd s(driver->fluid_grid->element_count);

    std::vector<int> e_faces;
    int f;
    int int_faces;
    double tmp_s;
    for (int e=0; e<driver->fluid_grid->element_count; e++) {
        e_faces = driver->fluid_grid->getElementFaces(e);
        tmp_s = 0;
        int_faces = 0;
        for (int i = 0; i < e_faces.size(); i++) {
            f = e_faces[i];

            //calculate equivalent shock sensor value for each method

            bool USE_CONSTANT_S = false;
            bool USE_MCCOR_COL = false;
            bool USE_COUCHMAN = false;
            bool USE_MCCOR_COL_3 = false;
            bool USE_CUSTOM = true;
            bool USE_ELEMENT_S = false;

            //get elements and face normal
            int e_minus = driver->fluid_grid->getOrientedElementsByFace(f)[0];
            int e_plus = driver->fluid_grid->getOrientedElementsByFace(f)[1];
            KinematicVector normal = driver->fluid_grid->getFaceNormal(job, f);

            if (e_plus > -1 && e_minus > -1) {
                int_faces += 1;

                //calculate max wave speed
                double rho_plus = driver->fluid_body->rho(e_plus);
                double rho_minus = driver->fluid_body->rho(e_minus);
                KinematicVector u_plus = driver->fluid_body->p(e_plus) / rho_plus;
                KinematicVector u_minus = driver->fluid_body->p(e_minus) / rho_minus;
                double E_plus = driver->fluid_body->rhoE(e_plus) / rho_plus;
                double E_minus = driver->fluid_body->rhoE(e_minus) / rho_minus;

                //scale rho_plus and rho_minus
                rho_plus /= driver->fluid_body->n_e(e_plus);
                rho_minus /= driver->fluid_body->n_e(e_minus);
                double n_bar = (driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e(e_minus)) / 2.0;

                if (USE_MCCOR_COL) { //USE McCorquodale and Collelo shock capturing method
                    double c_plus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                            n_bar * rho_plus,
                                                                            n_bar * rho_plus * u_plus,
                                                                            rho_plus * E_plus,
                                                                            n_bar);

                    double c_minus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                             n_bar * rho_minus,
                                                                             n_bar * rho_minus * u_minus,
                                                                             rho_minus * E_minus,
                                                                             n_bar);

                    double h_lambda = (u_plus - u_minus).dot(normal); //?
                    double beta = 0.3;
                    double c_min = std::min(c_plus, c_minus);

                    if (h_lambda < 0) {
                        tmp_s += std::min(1.0, h_lambda * h_lambda / (beta * c_min * c_min));
                    } else {
                        //do nothing
                    }

                } else if (USE_CONSTANT_S) { //USE uniform artificial viscosity
                    tmp_s += 1.0;
                } else if (USE_COUCHMAN) {
                    double c_plus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                            n_bar * rho_plus,
                                                                            n_bar * rho_plus * u_plus,
                                                                            rho_plus * E_plus,
                                                                            n_bar);

                    double c_minus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                             n_bar * rho_minus,
                                                                             n_bar * rho_minus * u_minus,
                                                                             rho_minus * E_minus,
                                                                             n_bar);

                    KinematicVector h = driver->fluid_grid->getElementCentroid(job, e_plus) -
                                        driver->fluid_grid->getElementCentroid(job, e_minus);
                    double h_lambda = (u_plus + u_minus).dot(h) / (2.0 * h.norm()) *
                                      (u_plus.norm() / c_plus - u_minus.norm() / c_minus); //?
                    double c_max = std::max(c_plus + u_plus.norm(), c_minus + u_minus.norm());

                    if (h_lambda < 0) {
                        tmp_s += std::min(std::abs(h_lambda) / c_max, 1.0);
                    } else {
                        //do nothing
                    }
                } else if (USE_MCCOR_COL_3) {
                    double P_plus = driver->fluid_material->getPressure(job, driver,
                                                                        n_bar * rho_plus,
                                                                        n_bar * rho_plus * u_plus,
                                                                        n_bar * rho_plus * E_plus,
                                                                        n_bar); //n_q(q_list[q]));
                    double H_plus = E_plus + P_plus / rho_plus;

                    double P_minus = driver->fluid_material->getPressure(job, driver,
                                                                         n_bar * rho_minus,
                                                                         n_bar * rho_minus * u_minus,
                                                                         n_bar * rho_minus * E_minus,
                                                                         n_bar); //n_q(q_list[q]));
                    double H_minus = E_minus + P_minus / rho_minus;

                    //approximate Roe advective rate
                    double w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
                    double rho_bar = w_1_bar * w_1_bar;
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
                    } else if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
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
                    double lambda_max = (std::abs(u_bar.dot(normal)) + c);
                    double beta = 0.01;
                    double h_lambda = (u_plus - u_minus).dot(normal);

                    if (h_lambda < 0){
                        tmp_s += std::min(1.0, h_lambda*h_lambda/(beta * c * c));
                    } else {
                        //do nothing
                    }
                } else if (USE_CUSTOM) {
                    /*
                    double P_plus = driver->fluid_material->getPressure(job, driver,
                                                                        n_bar * rho_plus,
                                                                        n_bar * rho_plus * u_plus,
                                                                        n_bar * rho_plus * E_plus,
                                                                        n_bar); //n_q(q_list[q]));
                    double H_plus = E_plus + P_plus / rho_plus;

                    double P_minus = driver->fluid_material->getPressure(job, driver,
                                                                         n_bar * rho_minus,
                                                                         n_bar * rho_minus * u_minus,
                                                                         n_bar * rho_minus * E_minus,
                                                                         n_bar); //n_q(q_list[q]));
                    double H_minus = E_minus + P_minus / rho_minus;

                    //approximate Roe advective rate
                    double w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
                    double rho_bar = w_1_bar * w_1_bar;
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
                    } else if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
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
                    */

                    double c_plus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                            n_bar * rho_plus,
                                                                            n_bar * rho_plus * u_plus,
                                                                            rho_plus * E_plus,
                                                                            n_bar);

                    double c_minus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                             n_bar * rho_minus,
                                                                             n_bar * rho_minus * u_minus,
                                                                             rho_minus * E_minus,
                                                                             n_bar);

                    //find maximum lambda
                    double lambda_max = std::max(std::abs(u_plus.dot(normal)) + c_plus,
                                                 std::abs(u_minus.dot(normal)) + c_minus); //(std::abs(u_bar.dot(normal)) + c);
                    double beta = 0.3;
                    double h_lambda = (u_plus - u_minus).dot(normal);

                    if (h_lambda < 0){
                        double c_min = std::min(c_plus, c_minus);
                        double epsilon_v = -h_lambda + std::abs(c_plus - c_minus);
                        tmp_s += std::min(1.0, epsilon_v*epsilon_v/(beta * beta * c_min * c_min));
                        //tmp_s += std::min(1.0, epsilon_v / (beta * c_min));
                    } else {
                        //do nothing
                    }
                } else if (USE_ELEMENT_S){
                    double c_plus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                            n_bar * rho_plus,
                                                                            n_bar * rho_plus * u_plus,
                                                                            rho_plus * E_plus,
                                                                            n_bar);

                    double c_minus = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                             n_bar * rho_minus,
                                                                             n_bar * rho_minus * u_minus,
                                                                             rho_minus * E_minus,
                                                                             n_bar);

                    //find maximum lambda
                    double beta = 0.3;
                    KinematicVector h = driver->fluid_grid->getElementCentroid(job,e_plus) -
                                        driver->fluid_grid->getElementCentroid(job,e_minus);
                    double h_lambda = h.norm() * 0.5 * (driver->fluid_body->L[e_plus].trace()
                                                        + driver->fluid_body->L[e_minus].trace()); //(u_plus - u_minus).dot(normal);

                    if (h_lambda < 0){
                        double c_min = std::min(c_plus, c_minus);
                        double epsilon_v = -h_lambda + std::abs(c_plus - c_minus);
                        tmp_s += std::min(1.0, epsilon_v*epsilon_v/(beta * beta * c_min * c_min));
                        //return lambda_max * std::min(1.0, epsilon_v / (beta * c_min));
                    } else {
                        //do nothing
                    }
                }
            }
        }

        //average s
        s(e) = tmp_s / int_faces;
    }

    //write s to file
    driver->serializer->writeScalarArray(s, "shock");

    return;
}


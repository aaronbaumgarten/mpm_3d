//
// Created by aaron on 12/12/19.
// fvm_cartesian.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <Eigen/Dense>
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
#include "fvm_grids.hpp"

/*----------------------------------------------------------------------------*/

double FVMGridBase::getElementVolume(int e){
    return v_e(e);
}

int FVMGridBase::getElementTag(int e){
    //do nothing.
    return 0;
}

double FVMGridBase::getFaceArea(int f){
    //get normal definition
    return face_areas(f);
}

int FVMGridBase::getFaceTag(int f){
    //return boundary tag of faces
    return bc_info[f].tag;
}

KinematicVector FVMGridBase::getFaceNormal(Job* job, int f){
    //return oriented normal associated with face
    return face_normals[f];
}

std::vector<int> FVMGridBase::getElementFaces(int e){
    //get faces associated with each element
    return element_faces[e];
}

std::array<int,2> FVMGridBase::getOrientedElementsByFace(int f){
    //get oriented elements associated with given face
    // A -|-> B
    return face_elements[f];
}

std::vector<int> FVMGridBase::getElementNeighbors(int e){
    //get list of elements that share a node with e
    return element_neighbors[e];
}


KinematicVector FVMGridBase::getElementCentroid(Job* job, int e){
    return x_e(e);
}

KinematicVector FVMGridBase::getFaceCentroid(Job* job, int f){
    return x_f(f);
}

std::vector<int> FVMGridBase::getElementQuadraturePoints(int e){
    std::vector<int> tmp = std::vector<int>(qpe);
    for (int q = 0; q<qpe; q++){
        tmp[q] = e*qpe + q;
    }
    return tmp;
}

std::vector<int> FVMGridBase::getFaceQuadraturePoints(int f){
    std::vector<int> tmp = std::vector<int>(qpf);
    for (int q = 0; q<qpf; q++){
        tmp[q] = int_quad_count + f*qpf + q;
    }
    return tmp;
}

KinematicVector FVMGridBase::getQuadraturePosition(Job* job, int q){
    return x_q(q);
}

double FVMGridBase::getQuadratureWeight(int q){
    return w_q(q);
}

/*----------------------------------------------------------------------------*/
//functions for solving system of equations

void FVMGridBase::mapMixturePropertiesToQuadraturePoints(Job* job, FiniteVolumeDriver* driver){
    //need to map porosity and solid velocity to quadrature points for flux calculations

    //ask material to update body porosity field
    if (driver->fluid_material->calculatePorosity(job, driver) == 1){
        //successful calculation of porosity
        n_q = Q*driver->fluid_body->n;          //n_q = Q_iq * n_i
        gradn_q = gradQ*driver->fluid_body->n;  //gradn_q = gradQ_iq * n_i
    } else {
        //no porosity calculation performed
        //n_q.setConstant(1.0);
        //gradn_q.setZero();
    }

    //ask material to update solid velocity field
    if (driver->fluid_material->updateSolidPhaseVelocity(job, driver) == 1){
        //successful update of solid velocity
        v_sq = Q*driver->fluid_body->v_s;       //v_sq = Q_iq * v_si
    } else {
        //no porosity calculation performed
        //v_sq.setZero();
    }

    return;
}

void FVMGridBase::constructMomentumField(Job* job, FiniteVolumeDriver* driver){
    if (driver->ORDER == 1){
        driver->fluid_body->p_x.setZero(); //let momentum be constant within an element
    } else if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMGridBase does not currently implement higher ORDER momentum reconstruction."
                      << std::endl;
        }

        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        KinematicVector p_0, p, p_max, p_min, min_dif;
        min_dif = KinematicVector(job->JOB_TYPE);
        double tmp_dif, rho_0;
        for (int e = 0; e<element_count; e++){
            //least squares fit of u_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            p_0 = driver->fluid_body->p[e];
            rho_0 = driver->fluid_body->rho(e);
            p_max = p_0; p_min = p_0;
            min_dif.setZero();
            //loop over components of momentum
            for (int mom_index = 0; mom_index<GRID_DIM; mom_index++){
                //create system of equations
                for (int ii=0; ii<element_neighbors[e].size(); ii++){
                    //just fill in b vector
                    p = driver->fluid_body->p[element_neighbors[e][ii]];
                    b_e[e](ii) = p[mom_index] - p_0[mom_index];

                    //update maximum and minimum velocities
                    if (p[mom_index] > p_max[mom_index]){
                        p_max[mom_index] = p[mom_index];
                    } else if (p[mom_index] < p_min[mom_index]){
                        p_min[mom_index] = p[mom_index];
                    }
                }

                //solve for mom_pos component of gradient
                sol = A_inv[e]*b_e[e];
                for (int pos = 0; pos<GRID_DIM; pos++){
                    driver->fluid_body->p_x(e, mom_index, pos) = sol(pos);
                }

                //calculate minimum magnitude difference in velocity components
                min_dif[mom_index] = std::abs(p_0[mom_index] - p_min[mom_index]);
                if (std::abs(p_0[mom_index] - p_max[mom_index]) < min_dif[mom_index]) {
                    min_dif[mom_index] = std::abs(p_0[mom_index] - p_max[mom_index]);
                }
            }

            //limit gradient to ensure monotonicity
            for (int mom_index = 0; mom_index<GRID_DIM; mom_index++){
                //calculate maximum momentum change in cell
                tmp_dif = 0;
                for (int pos = 0; pos < GRID_DIM; pos++){
                    tmp_dif += hx[pos]*std::abs(driver->fluid_body->p_x(e, mom_index, pos)/2.0);
                }

                if (tmp_dif > min_dif[mom_index]){
                    //if maximum velocity change in cell is larger than maximum difference b/w neighbors
                    //need to scale gradient
                    for (int pos=0; pos<GRID_DIM; pos++){
                        driver->fluid_body->p_x(e, mom_index, pos) *= min_dif[mom_index]/tmp_dif;
                    }
                }
            }
        }
    }
    return;
}

void FVMGridBase::constructDensityField(Job* job, FiniteVolumeDriver* driver){
    if (driver->ORDER == 1){
        driver->fluid_body->rho_x.setZero(); //let density be constant within an element
    } if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMGridBase does not currently implement higher ORDER momentum reconstruction."
                      << std::endl;
        }

        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        double rho_0, rho, rho_max, rho_min, min_dif;
        double tmp_dif;
        for (int e = 0; e<element_count; e++){
            //least squares fit of u_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            rho_0 = driver->fluid_body->rho(e);
            rho_max = rho_0;
            rho_min = rho_0;

            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                //fill in b vector
                rho = driver->fluid_body->rho(element_neighbors[e][ii]);
                b_e[e](ii) = rho - rho_0;

                //update maximum and minimum velocities
                if (rho > rho_max){
                    rho_max = rho;
                } else if (rho < rho_min){
                    rho_min = rho;
                }
            }

            //solve for components of gradient
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                driver->fluid_body->rho_x(e, pos) = sol(pos);
            }

            //calculate minimum magnitude difference in velocity components
            min_dif = std::abs(rho_0 - rho_min);
            if (std::abs(rho_0 - rho_max) < min_dif) {
                min_dif = std::abs(rho_0 - rho_max);
            }

            //limit gradient to ensure monotonicity
            //calculate maximum velocity change in cell
            tmp_dif = 0;
            for (int pos = 0; pos < GRID_DIM; pos++){
                tmp_dif += hx[pos]*std::abs(driver->fluid_body->rho_x(e, pos)/2.0);
            }

            if (tmp_dif > min_dif){
                //if maximum velocity change in cell is larger than maximum difference b/w neighbors
                //need to scale gradient
                for (int pos=0; pos<GRID_DIM; pos++){
                    driver->fluid_body->rho_x(e, pos) *= min_dif/tmp_dif;
                }
            }
        }
    }
    return;
}

void FVMGridBase::constructEnergyField(Job* job, FiniteVolumeDriver* driver){
    if (driver->ORDER == 1){
        driver->fluid_body->rhoE_x.setZero(); //let energy be constant within an element
    } if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMGridBase does not currently implement higher ORDER momentum reconstruction."
                      << std::endl;
        }

        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        double rhoE_0, rhoE, rhoE_max, rhoE_min, min_dif;
        double tmp_dif;
        for (int e = 0; e<element_count; e++){
            //least squares fit of u_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            rhoE_0 = driver->fluid_body->rhoE(e);
            rhoE_max = rhoE_0;
            rhoE_min = rhoE_0;

            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                //fill in b vector
                rhoE = driver->fluid_body->rhoE(element_neighbors[e][ii]);
                b_e[e](ii) = rhoE - rhoE_0;

                //update maximum and minimum velocities
                if (rhoE > rhoE_max){
                    rhoE_max = rhoE;
                } else if (rhoE < rhoE_min){
                    rhoE_min = rhoE;
                }
            }

            //solve for components of gradient
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                driver->fluid_body->rhoE_x(e, pos) = sol(pos);
            }

            //calculate minimum magnitude difference in velocity components
            min_dif = std::abs(rhoE_0 - rhoE_min);
            if (std::abs(rhoE_0 - rhoE_max) < min_dif) {
                min_dif = std::abs(rhoE_0 - rhoE_max);
            }

            //limit gradient to ensure monotonicity
            //calculate maximum velocity change in cell
            tmp_dif = 0;
            for (int pos = 0; pos < GRID_DIM; pos++){
                tmp_dif += hx[pos]*std::abs(driver->fluid_body->rho_x(e, pos)/2.0);
            }

            if (tmp_dif > min_dif){
                //if maximum velocity change in cell is larger than maximum difference b/w neighbors
                //need to scale gradient
                for (int pos=0; pos<GRID_DIM; pos++){
                    driver->fluid_body->rhoE_x(e, pos) *= min_dif/tmp_dif;
                }
            }
        }
    }
    return;
}

KinematicTensorArray FVMGridBase::getVelocityGradients(Job* job, FiniteVolumeDriver* driver){
    //reconstruct velocity field
    KinematicTensorArray u_x = KinematicTensorArray(element_count, job->JOB_TYPE);
    KinematicVector u = KinematicVector(job->JOB_TYPE);
    KinematicVector p = KinematicVector(job->JOB_TYPE);
    double rho, rho_0, tmp_dif;

    Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
    KinematicVector x_0, x;
    KinematicVector u_0;

    for (int e=0; e<element_count; e++){
        p = driver->fluid_body->p[e];
        rho_0 = driver->fluid_body->rho(e);
        u_0 = p/rho_0;

        //least squares fit of u_x to neighbors of element e
        x_0 = getElementCentroid(job, e);
        for (int dir = 0; dir<GRID_DIM; dir++){
            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                rho = driver->fluid_body->rho(element_neighbors[e][ii]);
                p = driver->fluid_body->p[element_neighbors[e][ii]];
                b_e[e](ii) = p[dir]/rho - u_0[dir];
            }

            //solve for u_pos component of gradient
            //sol = A_e[e].householderQr().solve(b_e[e]);
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                u_x(e, dir, pos) = sol(pos);
            }
        }

    }
    return u_x;
}

/*----------------------------------------------------------------------------*/
//functions to compute element-wise mass flux
Eigen::VectorXd FVMGridBase::calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver){
    //flux rates
    Eigen::VectorXd result = Eigen::VectorXd(element_count);
    result.setZero();

    int e_plus, e_minus;    //elements
    double flux, area;      //face values

    //quadrature
    KinematicVector x = KinematicVector(job->JOB_TYPE);
    std::vector<int> q_list;

    //isothermal flux function
    KinematicVector p_plus, p_minus;
    KinematicVector u_plus, u_minus, normal, u_bar; //velocity
    double rho_plus, rho_minus, rho_bar;
    double lambda_1, lambda_2, a_1, a_2, c;

    //thermal flux function
    double rhoE_plus, rhoE_minus;
    double E_plus, E_minus, H_plus, H_minus, P_plus, P_minus, H_bar;
    double lambda_3, a_3;
    double w_1_bar, u_bar_squared;


    //loop over faces and use quadrature to reconstruct flux integral
    u_plus = KinematicVector(job->JOB_TYPE);
    u_minus = KinematicVector(job->JOB_TYPE);
    u_bar = KinematicVector(job->JOB_TYPE);

    for (int f = 0; f < face_count; f++) {
        //face dimensions
        area = getFaceArea(f);
        normal = getFaceNormal(job, f);
        q_list = getFaceQuadraturePoints(f);

        //flux calculation depends on whether face is on boundary
        if (bc_info[f].tag == -1 || (bc_info[f].tag == PERIODIC && face_elements[f][0] > -1)) {
            //face is interior to domain; no BCs (or periodic)
            e_minus = face_elements[f][0];
            e_plus = face_elements[f][1];

            //different flux functions for different simulation TYPES
            if (driver->TYPE == FiniteVolumeDriver::INCOMPRESSIBLE){
                std::cerr << "FVMGridBase does not have INCOMPRESSIBLE flux function implemented. Exiting." << std::endl;
                exit(0);
            } else if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
                //loop over quadrature points
                for (int q = 0; q < q_list.size(); q++) {
                    //relative position to centroid of A
                    if (bc_info[f].tag == PERIODIC) {
                        //by convention, centroid of A will be on other side of domain
                        //so use relative distance to B but flip normal direction
                        x = x_q[q_list[q]] - x_e[e_plus];
                        x -= 2.0 * (x.dot(normal)) * normal;
                    } else {
                        x = x_q(q_list[q]) - x_e(e_minus);
                    }

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    u_minus = (driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x) / rho_minus;

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x) / rho_plus;

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) /
                            (std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = std::sqrt(rho_plus * rho_minus);
                    rhoE_minus = driver->fluid_body->rhoE(e_minus);
                    c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                                                                n_q(q_list[q]));

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    if (lambda_1 < delta * c) {
                        lambda_1 = 0.5 * ((lambda_1 * lambda_1) / (delta * c) + (delta * c));
                    }
                    if (lambda_2 < delta * c) {
                        lambda_2 = 0.5 * ((lambda_2 * lambda_2) / (delta * c) + (delta * c));
                    }

                    //calculate Roe eigenvector coefficients
                    a_1 = 0.5 * (rho_plus - rho_minus) - 1.0 / (2.0 * c) * (rho_plus * u_plus - rho_minus * u_minus -
                                                                            (rho_plus - rho_minus) * u_bar).dot(normal);
                    a_2 = (rho_plus - rho_minus) - a_1;

                    //flux in n direction
                    flux = w_q(q_list[q]) * 0.5 * (rho_plus * u_plus.dot(normal) + rho_minus * u_minus.dot(normal)
                                                   - a_1 * lambda_1 - a_2 * lambda_2);

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else if (driver->TYPE == FiniteVolumeDriver::THERMAL) {
                //loop over quadrature points
                for (int q = 0; q < q_list.size(); q++) {
                    //relative position to centroid of A
                    if (bc_info[f].tag == PERIODIC) {
                        //by convention, centroid of A will be on other side of domain
                        //so use relative distance to B but flip normal direction
                        x = x_q[q_list[q]] - x_e[e_plus];
                        x -= 2.0 * (x.dot(normal)) * normal;
                    } else {
                        x = x_q(q_list[q]) - x_e(e_minus);
                    }

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x;
                    u_minus = p_minus / rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    E_minus = rhoE_minus/rho_minus;
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    H_minus = E_minus + P_minus/rho_minus;

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    E_plus = rhoE_plus/rho_plus;
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    H_minus = E_plus + P_plus/rho_plus;

                    //approximate Roe advective rate
                    w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus))/2.0;
                    rho_bar = w_1_bar*w_1_bar;
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) / (2.0 * w_1_bar);
                    u_bar_squared = u_bar.dot(u_bar);
                    H_bar = (std::sqrt(rho_plus) * H_plus + std::sqrt(rho_minus) * H_minus) / (2.0 * w_1_bar);
                    c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                                                                n_q(q_list[q]));

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    lambda_3 = std::abs(u_bar.dot(normal));

                    //Harten's entropy fix
                    if (lambda_1 < delta * c) {
                        lambda_1 = 0.5 * ((lambda_1 * lambda_1) / (delta * c) + (delta * c));
                    }
                    if (lambda_2 < delta * c) {
                        lambda_2 = 0.5 * ((lambda_2 * lambda_2) / (delta * c) + (delta * c));
                    }

                    //calculate Roe eigenvector coefficients
                    a_3 = ((H_bar - u_bar_squared)*(rho_plus - rho_minus)
                            + u_bar.dot(p_plus - p_minus)
                            - (rhoE_plus - rhoE_minus)) / (H_bar - 0.5*u_bar_squared);

                    a_1 = 0.5*((rho_plus - rho_minus)
                                - a_3
                                - (p_plus - p_minus - u_bar*(rho_plus - rho_minus)).dot(normal) / c);

                    a_2 = (rho_plus - rho_minus) - a_3 - a_1;


                    //flux in n direction
                    flux = w_q(q_list[q]) * 0.5 * (rho_plus * u_plus.dot(normal) + rho_minus * u_minus.dot(normal)
                                                   - a_1 * lambda_1 - a_2 * lambda_2 - a_3*lambda_3);

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else {
                std::cerr << "FVMGridBase does not have flux function implemented for TYPE " << driver->TYPE << "! Exiting." << std::endl;
                exit(0);
            }
        } else if (bc_info[f].tag == VELOCITY_INLET) {
            //face has prescribed velocity
            //reconstruct density field at quadrature points
            for (int q=0; q<q_list.size(); q++){
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    flux = rho_minus * w_q(q_list[q]) * bc_info[f].vector.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    flux = rho_plus * w_q(q_list[q]) * bc_info[f].vector.dot(normal);
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == VELOCITY_DENSITY_INLET){
            //density and velocity given
            for (int q=0; q<q_list.size(); q++){
                //flux = rho * u.dot(n)
                flux = bc_info[f].values[0] * w_q(q_list[q]) * bc_info[f].vector.dot(normal);

                //check which side of face has valid element
                if (e_minus > -1){
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    result(e_plus) += flux;
                }
            }
        } else if(bc_info[f].tag == VELOCITY_TEMP_INLET){
            //velocity and temperature given
            //reconstruct density field at quadrature points
            for (int q=0; q<q_list.size(); q++){
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    flux = rho_minus * w_q(q_list[q]) * bc_info[f].vector.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    flux = rho_plus * w_q(q_list[q]) * bc_info[f].vector.dot(normal);
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PRESSURE_INLET) {
            //face has prescribed pressure (and temperature)
            //reconstruct velocity field at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                       bc_info[f].values[0],
                                                                                       bc_info[f].values[1],
                                                                                       n_q(q_list[q]));
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    u_minus = driver->fluid_body->p(e_minus) + driver->fluid_body->p_x[e_minus]*x;
                    u_minus /= rho_minus;
                    flux = rho_bar * w_q(q_list[q]) * u_minus.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = driver->fluid_body->p(e_plus) + driver->fluid_body->p_x[e_plus]*x;
                    u_plus /= rho_plus;
                    flux = rho_bar * w_q(q_list[q]) * u_plus.dot(normal);
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PRESSURE_OUTLET || bc_info[f].tag == DAMPED_OUTLET) {
            //face has prescribed pressure (and temperature)
            //reconstruct density and velocity (but adjust if flow direction changes)
            for (int q = 0; q < q_list.size(); q++) {
                rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                       bc_info[f].values[0],
                                                                                       bc_info[f].values[1],
                                                                                       n_q(q_list[q]));
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    u_minus = driver->fluid_body->p(e_minus) + driver->fluid_body->p_x[e_minus]*x;
                    u_minus /= rho_minus;
                    if (u_minus.dot(normal) < 0){
                        //flow direction is into A
                        flux = rho_bar * w_q(q_list[q]) * u_minus.dot(normal);
                    } else {
                        //flow direction is out of A
                        flux = rho_minus * w_q(q_list[q]) * u_minus.dot(normal);
                    }
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    u_plus = driver->fluid_body->p(e_plus) + driver->fluid_body->p_x[e_plus]*x;
                    u_plus /= rho_bar;
                    if (u_plus.dot(normal) > 0){
                        //flow direction is into B
                        flux = rho_bar * w_q(q_list[q]) * u_plus.dot(normal);
                    } else {
                        //flow direction is out of B
                        flux = rho_plus * w_q(q_list[q]) * u_plus.dot(normal);
                    }
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == ADIABATIC_WALL ||
                    bc_info[f].tag == THERMAL_WALL ||
                    bc_info[f].tag == SYMMETRIC_WALL){
            //do nothing, no flux across this boundary
        } else if (bc_info[f].tag == SUPERSONIC_INLET){
            //density and velocity given
            for (int q=0; q<q_list.size(); q++){
                //flux = rho * u.dot(n)
                flux = bc_info[f].values[0] * w_q(q_list[q]) * bc_info[f].vector.dot(normal);

                //check which side of face has valid element
                if (e_minus > -1){
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == SUPERSONIC_OUTLET){
            //reconstruct flux (CAREFUL! back-flow unstable)
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    u_minus = driver->fluid_body->p(e_minus) + driver->fluid_body->p_x[e_minus]*x;
                    u_minus /= rho_minus;
                    flux = rho_minus * w_q(q_list[q]) * u_minus.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = driver->fluid_body->p(e_plus) + driver->fluid_body->p_x[e_plus]*x;
                    u_plus /= rho_plus;
                    flux = rho_plus * w_q(q_list[q]) * u_plus.dot(normal);
                    result(e_plus) += flux;
                }
            }
        } else {
            std::cerr << "ERROR! FVMGridBase does not have flux defined for bc tag " << bc_info[f].tag << "!" << std::endl;
            //don't exit, without flux defined, essentially a wall
        }
    }
    return result;
}

//functions to compute element momentum fluxes (including tractions)
KinematicVectorArray FVMGridBase::calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver){
    //flux rates
    KinematicVectorArray result = KinematicVectorArray(element_count, job->JOB_TYPE);
    result.setZero();

    int e_plus, e_minus;    //elements
    double area;      //face values
    KinematicVector flux = KinematicVector(job->JOB_TYPE);

    //quadrature
    KinematicVector x = KinematicVector(job->JOB_TYPE);
    std::vector<int> q_list;

    //isothermal flux function
    KinematicVector p_plus, p_minus;
    KinematicVector u_plus, u_minus, normal, u_bar; //velocity
    double rho_plus, rho_minus, rho_bar;
    double lambda_1, lambda_2, a_1, a_2, lambda_4, c;
    KinematicVector a_4;

    //thermal flux function
    double rhoE_plus, rhoE_minus, rhoE_bar;
    double E_plus, E_minus, H_plus, H_minus, P_plus, P_minus, H_bar;
    double lambda_3, a_3;
    double w_1_bar, u_bar_squared;

    //traction calculatons
    KinematicTensorArray L = getVelocityGradients(job, driver);
    KinematicTensor L_tmp = KinematicTensor(job->JOB_TYPE);
    MaterialTensor tau_plus, tau_minus;


    //loop over faces and use quadrature to reconstruct flux integral
    u_plus = KinematicVector(job->JOB_TYPE);
    u_minus = KinematicVector(job->JOB_TYPE);
    u_bar = KinematicVector(job->JOB_TYPE);

    for (int f = 0; f < face_count; f++) {
        //face dimensions
        area = getFaceArea(f);
        normal = getFaceNormal(job, f);
        q_list = getFaceQuadraturePoints(f);

        //flux calculation depends on whether face is on boundary
        if (bc_info[f].tag == -1 || (bc_info[f].tag == PERIODIC && face_elements[f][0] > -1)) {
            //face is interior to domain; no BCs (or periodic)
            e_minus = face_elements[f][0];
            e_plus = face_elements[f][1];

            //different flux functions for different simulation TYPES
            if (driver->TYPE == FiniteVolumeDriver::INCOMPRESSIBLE) {
                std::cerr << "FVMGridBase does not have INCOMPRESSIBLE flux function implemented. Exiting."
                          << std::endl;
                exit(0);
            } else if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
                //loop over quadrature points
                for (int q = 0; q < q_list.size(); q++) {
                    //relative position to centroid of A
                    if (bc_info[f].tag == PERIODIC) {
                        //by convention, centroid of A will be on other side of domain
                        //so use relative distance to B but flip normal direction
                        x = x_q[q_list[q]] - x_e[e_plus];
                        x -= 2.0 * (x.dot(normal)) * normal;
                    } else {
                        x = x_q(q_list[q]) - x_e(e_minus);
                    }

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x;
                    u_minus = p_minus / rho_minus;

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) /
                            (std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = std::sqrt(rho_plus * rho_minus);
                    rhoE_minus = driver->fluid_body->rhoE(e_minus);
                    c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                                                                n_q(q_list[q]));

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    if (lambda_1 < delta*c){
                        lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                    }
                    if (lambda_2 < delta*c){
                        lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                    }
                    lambda_4 = std::abs(u_bar.dot(normal));

                    //calculate Roe eigenvector coefficients
                    a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                    a_2 = (rho_plus - rho_minus) - a_1;
                    a_4 = p_plus - p_minus - (rho_plus-rho_minus)*u_bar;
                    a_4 = a_4 - a_4.dot(normal)*normal; //remove normal component of a_4 vector

                    //flux in n direction
                    flux = w_q(q_list[q]) * 0.5 *(p_plus*u_plus.dot(normal) + p_minus*u_minus.dot(normal)
                                                - a_1*lambda_1*(u_bar - c*normal)
                                                - a_2*lambda_2*(u_bar + c*normal)
                                                - a_4*lambda_4);

                    //add tractions to momentum flux
                    //for now use simple reconstruction of theta
                    tau_minus = driver->fluid_material->getShearStress(job, driver, L[e_minus], rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    tau_plus = driver->fluid_material->getShearStress(job, driver, L[e_plus], rho_plus, p_plus, rhoE_minus, n_q(q_list[q]));
                    P_plus = c*c*(rho_plus - rho_bar) + driver->fluid_material->getPressure(job, driver, rho_bar, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = c*c*(rho_minus - rho_plus) + P_plus;

                    flux += w_q(q_list[q]) * 0.5 * ((P_plus + P_minus)*normal - KinematicVector((tau_plus + tau_minus)*normal, job->JOB_TYPE));

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else if (driver->TYPE == FiniteVolumeDriver::THERMAL) {
                //loop over quadrature points
                for (int q = 0; q < q_list.size(); q++) {
                    //relative position to centroid of A
                    if (bc_info[f].tag == PERIODIC) {
                        //by convention, centroid of A will be on other side of domain
                        //so use relative distance to B but flip normal direction
                        x = x_q[q_list[q]] - x_e[e_plus];
                        x -= 2.0 * (x.dot(normal)) * normal;
                    } else {
                        x = x_q(q_list[q]) - x_e(e_minus);
                    }

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x;
                    u_minus = p_minus / rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    E_minus = rhoE_minus / rho_minus;
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus,
                                                                  n_q(q_list[q]));
                    H_minus = E_minus + P_minus / rho_minus;

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus)  + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    E_plus = rhoE_plus / rho_plus;
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus,
                                                                 n_q(q_list[q]));
                    H_minus = E_plus + P_plus / rho_plus;

                    //approximate Roe advective rate
                    w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
                    rho_bar = w_1_bar * w_1_bar;
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) / (2.0 * w_1_bar);
                    u_bar_squared = u_bar.dot(u_bar);
                    H_bar = (std::sqrt(rho_plus) * H_plus + std::sqrt(rho_minus) * H_minus) / (2.0 * w_1_bar);
                    c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                                                                n_q(q_list[q]));

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    lambda_3 = std::abs(u_bar.dot(normal));
                    lambda_4 = std::abs(u_bar.dot(normal));

                    //Harten's entropy fix
                    if (lambda_1 < delta * c) {
                        lambda_1 = 0.5 * ((lambda_1 * lambda_1) / (delta * c) + (delta * c));
                    }
                    if (lambda_2 < delta * c) {
                        lambda_2 = 0.5 * ((lambda_2 * lambda_2) / (delta * c) + (delta * c));
                    }

                    //calculate Roe eigenvector coefficients
                    a_3 = ((H_bar - u_bar_squared) * (rho_plus - rho_minus)
                           + u_bar.dot(p_plus - p_minus)
                           - (rhoE_plus - rhoE_minus)) / (H_bar - 0.5 * u_bar_squared);

                    a_1 = 0.5 * ((rho_plus - rho_minus)
                                 - a_3
                                 - (p_plus - p_minus - u_bar * (rho_plus - rho_minus)).dot(normal) / c);

                    a_2 = (rho_plus - rho_minus) - a_3 - a_1;

                    a_4 = (p_plus - p_minus) - u_bar*(rho_plus - rho_minus);
                    a_4 -= a_4.dot(normal)*normal;                              //remove normal component of vector

                    //flux in n direction
                    flux = w_q(q_list[q]) * 0.5 *(p_plus*u_plus.dot(normal) + p_minus*u_minus.dot(normal)
                                                - a_1*lambda_1*(u_bar - c*normal)
                                                - a_2*lambda_2*(u_bar + c*normal)
                                                - a_3*lambda_3*(u_bar)
                                                - a_4*lambda_4);

                    //add tractions to momentum flux
                    tau_minus = driver->fluid_material->getShearStress(job, driver, L[e_minus], rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    tau_plus = driver->fluid_material->getShearStress(job, driver, L[e_plus], rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux += w_q(q_list[q]) * 0.5 * ((P_plus + P_minus)*normal - KinematicVector((tau_plus + tau_minus)*normal, job->JOB_TYPE));

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else {
                std::cerr << "FVMGridBase does not have flux function implemented for TYPE " << driver->TYPE
                          << "! Exiting." << std::endl;
                exit(0);
            }
        } else if (bc_info[f].tag == VELOCITY_INLET) {
            //face has prescribed velocity
            //reconstruct density and traction field at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = bc_info[f].vector * rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rho_minus*bc_info[f].vector*bc_info[f].vector.dot(normal)
                                             + P_minus*normal
                                             - KinematicVector(tau_minus*normal, job->JOB_TYPE));

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = bc_info[f].vector*rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rho_plus*bc_info[f].vector*bc_info[f].vector.dot(normal)
                                             + P_plus*normal
                                             - KinematicVector(tau_plus*normal, job->JOB_TYPE));

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == VELOCITY_DENSITY_INLET) {
            //density and velocity given
            //reconstruct traction field at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    p_minus = bc_info[f].vector*bc_info[f].values[0]; //p = u * rho
                    rho_minus = bc_info[f].values[0];
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (p_minus*bc_info[f].vector.dot(normal)
                                             + P_minus*normal
                                             - KinematicVector(tau_minus*normal, job->JOB_TYPE));

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    p_plus = bc_info[f].vector * bc_info[f].values[0]; //p = rho * u
                    rho_plus = bc_info[f].values[0];
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (p_plus*bc_info[f].vector.dot(normal)
                                             + P_plus*normal
                                             - KinematicVector(tau_plus*normal, job->JOB_TYPE));

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == VELOCITY_TEMP_INLET) {
            //velocity and temperature given
            //reconstruct traction field at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = rho_minus * bc_info[f].vector;
                    P_minus = driver->fluid_material->getPressureFromDensityAndTemperature(job, driver, rho_minus, bc_info[f].values[0], n_q(q_list[q]));

                    //rhoE = rho * eps + 0.5*rho*u^2
                    rhoE_minus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver, P_minus, bc_info[f].values[0], n_q(q_list[q]));
                    rhoE_minus += 0.5*rho_minus*bc_info[f].vector.dot(bc_info[f].vector);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (p_minus*bc_info[f].vector.dot(normal)
                                             + P_minus*normal
                                             - KinematicVector(tau_minus*normal, job->JOB_TYPE));

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = rho_plus * bc_info[f].vector;
                    P_plus = driver->fluid_material->getPressureFromDensityAndTemperature(job, driver, rho_plus, bc_info[f].values[0], n_q(q_list[q]));

                    //rhoE = rho * eps + 0.5*rho*u^2
                    rhoE_plus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver, P_plus, bc_info[f].values[0], n_q(q_list[q]));
                    rhoE_plus += 0.5*rho_plus*bc_info[f].vector.dot(bc_info[f].vector);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (p_plus*bc_info[f].vector.dot(normal)
                                             + P_plus*normal
                                             - KinematicVector(tau_plus*normal, job->JOB_TYPE));

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PRESSURE_INLET) {
            //face has prescribed pressure (and temperature)
            //reconstruct velocity field at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                       bc_info[f].values[0],
                                                                                       bc_info[f].values[1],
                                                                                       n_q(q_list[q]));

                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;
                    u_minus = p_minus/rho_minus;

                    //tau_minus = 0
                    P_minus = bc_info[f].values[0];

                    flux = w_q(q_list[q]) * (rho_bar*u_minus*u_minus.dot(normal)
                                             + P_minus*normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;
                    u_plus = p_plus/rho_plus;

                    //tau_minus = 0
                    P_plus = bc_info[f].values[0];

                    flux = w_q(q_list[q]) * (rho_bar*u_plus*u_plus.dot(normal)
                                             + P_plus*normal);

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PRESSURE_OUTLET || bc_info[f].tag == DAMPED_OUTLET) {
            //face has prescribed pressure (and temperature)
            //reconstruct density and velocity (but adjust if flow direction changes)
            for (int q = 0; q < q_list.size(); q++) {
                rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                       bc_info[f].values[0],
                                                                                       bc_info[f].values[1],
                                                                                       n_q(q_list[q]));

                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x;
                    u_minus = p_minus / rho_minus;

                    //tau_minus = 0
                    P_minus = bc_info[f].values[0];
                    if (bc_info[f].tag == DAMPED_OUTLET){
                        rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x(e_minus).dot(x);

                        P_minus = (1.0 - damping_coefficient)*P_minus
                                 + damping_coefficient*driver->fluid_material->getPressure(job,driver,
                                                                                           rho_minus,
                                                                                           p_minus,
                                                                                           rhoE_minus,n_q(q_list[q]));
                    }

                    if (u_minus.dot(normal) > 0) {
                        //flow out of A
                        flux = w_q(q_list[q]) * (rho_minus * u_minus * u_minus.dot(normal)
                                                 + P_minus * normal);
                    } else {
                        //flow into A
                        flux = w_q(q_list[q]) * (rho_bar * u_minus * u_minus.dot(normal)
                                                 + P_minus * normal);
                    }


                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;

                    //tau_minus = 0
                    P_plus = bc_info[f].values[0];
                    if (bc_info[f].tag == DAMPED_OUTLET){
                        rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x(e_plus).dot(x);

                        P_plus = (1.0 - damping_coefficient)*P_plus
                                 + damping_coefficient*driver->fluid_material->getPressure(job,driver,
                                                                                           rho_plus,
                                                                                           p_plus,
                                                                                           rhoE_plus,n_q(q_list[q]));
                    }

                    if (u_plus.dot(normal) < 0) {
                        //flow out of B
                        flux = w_q(q_list[q]) * (rho_plus * u_plus * u_plus.dot(normal)
                                                 + P_plus * normal);
                    } else {
                        //flow into B
                        flux = w_q(q_list[q]) * (rho_bar * u_plus * u_plus.dot(normal)
                                                 + P_plus * normal);
                    }

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == ADIABATIC_WALL ||
                   bc_info[f].tag == THERMAL_WALL){
            //reconstruct pressure and traction at boundary
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = bc_info[f].vector * rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (P_minus*normal
                                             - KinematicVector(tau_minus*normal, job->JOB_TYPE));

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = bc_info[f].vector*rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (P_plus*normal
                                             - KinematicVector(tau_plus*normal, job->JOB_TYPE));

                    result(e_plus) += flux;
                }
            }
        } else if(bc_info[f].tag == SYMMETRIC_WALL) {
            //reconstruct pressure and traction at boundary
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = bc_info[f].vector * rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);

                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (P_minus*normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = bc_info[f].vector*rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);

                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (P_plus*normal);

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == SUPERSONIC_INLET) {
            //density and velocity given (and temperature if thermal)
            for (int q = 0; q < q_list.size(); q++) {
                //neglect shear stress at inlet
                P_minus = driver->fluid_material->getPressureFromDensityAndTemperature(job, driver,
                                                                                       bc_info[f].values[0],
                                                                                       bc_info[f].values[1],
                                                                                       n_q(q_list[q]));

                flux = w_q(q_list[q]) * (bc_info[f].values[0] * bc_info[f].vector * bc_info[f].vector.dot(normal)
                                         + P_minus * normal);

                //check which side of face has valid element
                if (e_minus > -1) {
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == SUPERSONIC_OUTLET) {
            //reconstruct flux (CAREFUL! back-flow unstable)
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p(e_minus) + driver->fluid_body->p_x[e_minus] * x;
                    u_minus = p_minus/rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (p_minus * u_minus.dot(normal)
                                             + P_minus * normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p(e_plus) + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus/rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (p_plus * u_plus.dot(normal)
                                             + P_plus * normal);
                    result(e_plus) += flux;
                }
            }
        } else {
            std::cerr << "ERROR! FVMGridBase does not have flux defined for bc tag " << bc_info[f].tag << "!"
                      << std::endl;
            //don't exit, without flux defined, essentially a wall
        }
    }
    return result;
}


//functions to compute element momentum fluxes (including tractions)
Eigen::VectorXd FVMGridBase::calculateElementEnergyFluxes(Job* job, FiniteVolumeDriver* driver){
    //don't calculate for ISOTHERMAL
    if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL){
        return Eigen::VectorXd::Zero(element_count);
    }

    //flux rates
    Eigen::VectorXd result = Eigen::VectorXd(element_count);
    result.setZero();

    int e_plus, e_minus;    //elements
    double area, flux;      //face values

    //quadrature
    KinematicVector x = KinematicVector(job->JOB_TYPE);
    std::vector<int> q_list;

    //isothermal flux function
    KinematicVector p_plus, p_minus;
    KinematicVector u_plus, u_minus, normal, u_bar; //velocity
    double rho_plus, rho_minus, rho_bar;
    double lambda_1, lambda_2, a_1, a_2, lambda_4, c;
    KinematicVector a_4, s;

    //thermal flux function
    double rhoE_plus, rhoE_minus, rhoE_bar;
    double E_plus, E_minus, H_plus, H_minus, P_plus, P_minus, H_bar;
    double lambda_3, a_3;
    double w_1_bar, u_bar_squared;
    double theta_plus, theta_minus, theta_bar;
    KinematicVector theta_x, heat_flux;

    //traction calculatons
    KinematicTensorArray L = getVelocityGradients(job, driver);
    KinematicTensor L_tmp = KinematicTensor(job->JOB_TYPE);
    MaterialTensor tau_plus, tau_minus;

    //loop over faces and use quadrature to reconstruct flux integral
    u_plus = KinematicVector(job->JOB_TYPE);
    u_minus = KinematicVector(job->JOB_TYPE);
    u_bar = KinematicVector(job->JOB_TYPE);

    for (int f = 0; f < face_count; f++) {
        //face dimensions
        area = getFaceArea(f);
        normal = getFaceNormal(job, f);
        q_list = getFaceQuadraturePoints(f);

        //flux calculation depends on whether face is on boundary
        if (bc_info[f].tag == -1 || (bc_info[f].tag == PERIODIC && face_elements[f][0] > -1)) {
            //face is interior to domain; no BCs (or periodic)
            e_minus = face_elements[f][0];
            e_plus = face_elements[f][1];

            //different flux functions for different simulation TYPES
            if (driver->TYPE == FiniteVolumeDriver::INCOMPRESSIBLE) {
                std::cerr << "FVMGridBase does not have INCOMPRESSIBLE flux function implemented. Exiting."
                          << std::endl;
                exit(0);
            } else if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
                //nothing to do here
            } else if (driver->TYPE == FiniteVolumeDriver::THERMAL) {
                //loop over quadrature points
                for (int q = 0; q < q_list.size(); q++) {
                    //relative position to centroid of A
                    if (bc_info[f].tag == PERIODIC) {
                        //by convention, centroid of A will be on other side of domain
                        //so use relative distance to B but flip normal direction
                        x = x_q[q_list[q]] - x_e[e_plus];
                        x -= 2.0 * (x.dot(normal)) * normal;
                    } else {
                        x = x_q(q_list[q]) - x_e(e_minus);
                    }

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x;
                    u_minus = p_minus / rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    E_minus = rhoE_minus / rho_minus;
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus,
                                                                  n_q(q_list[q]));
                    H_minus = E_minus + P_minus / rho_minus;
                    tau_minus = driver->fluid_material->getShearStress(job, driver, L[e_minus], rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus)  + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    E_plus = rhoE_plus / rho_plus;
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus,
                                                                 n_q(q_list[q]));
                    H_minus = E_plus + P_plus / rho_plus;
                    tau_plus = driver->fluid_material->getShearStress(job, driver, L[e_plus], rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    //approximate Roe advective rate
                    w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
                    rho_bar = w_1_bar * w_1_bar;
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) / (2.0 * w_1_bar);
                    u_bar_squared = u_bar.dot(u_bar);
                    H_bar = (std::sqrt(rho_plus) * H_plus + std::sqrt(rho_minus) * H_minus) / (2.0 * w_1_bar);
                    c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                                                                n_q(q_list[q]));

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    lambda_3 = std::abs(u_bar.dot(normal));
                    lambda_4 = std::abs(u_bar.dot(normal));

                    //Harten's entropy fix
                    if (lambda_1 < delta * c) {
                        lambda_1 = 0.5 * ((lambda_1 * lambda_1) / (delta * c) + (delta * c));
                    }
                    if (lambda_2 < delta * c) {
                        lambda_2 = 0.5 * ((lambda_2 * lambda_2) / (delta * c) + (delta * c));
                    }

                    //calculate Roe eigenvector coefficients
                    a_3 = ((H_bar - u_bar_squared) * (rho_plus - rho_minus)
                           + u_bar.dot(p_plus - p_minus)
                           - (rhoE_plus - rhoE_minus)) / (H_bar - 0.5 * u_bar_squared);

                    a_1 = 0.5 * ((rho_plus - rho_minus)
                                 - a_3
                                 - (p_plus - p_minus - u_bar * (rho_plus - rho_minus)).dot(normal) / c);

                    a_2 = (rho_plus - rho_minus) - a_3 - a_1;

                    a_4 = (p_plus - p_minus) - u_bar*(rho_plus - rho_minus);
                    a_4 -= a_4.dot(normal)*normal;                              //remove normal component of vector

                    //tangent component of u_bar
                    s = u_bar - u_bar.dot(normal)*normal;

                    //flux in n direction
                    flux = w_q(q_list[q]) * 0.5 *((rhoE_plus + P_plus)*u_plus.dot(normal)
                                                  + (rhoE_minus + P_minus)*u_minus.dot(normal)
                                                  - KinematicVector((tau_plus*u_plus),job->JOB_TYPE).dot(normal)
                                                  - KinematicVector((tau_minus*u_minus),job->JOB_TYPE).dot(normal)
                                                  - a_1*lambda_1*(H_bar - u_bar.dot(normal)*c)
                                                  - a_2*lambda_2*(H_bar + u_bar.dot(normal)*c)
                                                  - a_3*lambda_3*(0.5*u_bar_squared)
                                                  - a_4.dot(s)*lambda_4);

                    //calculate heat flux
                    theta_plus = driver->fluid_material->getTemperature(job, driver,
                                                                        driver->fluid_body->rho(e_plus),
                                                                        driver->fluid_body->p[e_plus],
                                                                        driver->fluid_body->rhoE(e_plus),
                                                                        n_q(q_list[q]));
                    theta_minus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_minus),
                                                                         driver->fluid_body->p[e_minus],
                                                                         driver->fluid_body->rhoE(e_minus),
                                                                         n_q(q_list[q]));
                    //relative position of centroid B from A
                    if (bc_info[f].tag == PERIODIC) {
                        //by convention, centroid of A will be on other side of domain
                        //so use twice relative distance to B in normal direction
                        x = x_q[q_list[q]] - x_e[e_plus];
                        x = -2.0 * (x.dot(normal)) * normal;
                    } else {
                        x = x_e(e_plus) - x_e(e_minus);
                    }
                    //estimate thermal gradient
                    theta_x = (theta_plus - theta_minus)*x/(x.dot(x));
                    theta_bar = theta_plus + theta_x.dot(x_q[q_list[q]] - x_e[e_plus]); //theta = theta_0 + dtheta/dx * dx
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_bar, theta_bar, theta_x, n_q(q_list[q]));
                    // add to flux calculation
                    flux += w_q(q_list[q])*heat_flux.dot(normal);

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else {
                std::cerr << "FVMGridBase does not have flux function implemented for TYPE " << driver->TYPE
                          << "! Exiting." << std::endl;
                exit(0);
            }
        } else if (bc_info[f].tag == VELOCITY_INLET) {
            //face has prescribed velocity
            //reconstruct density and traction field at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = bc_info[f].vector * rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rhoE_minus*bc_info[f].vector.dot(normal)
                                             + P_minus*bc_info[f].vector.dot(normal)
                                             - KinematicVector(tau_minus*bc_info[f].vector, job->JOB_TYPE)).dot(normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = bc_info[f].vector*rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rhoE_plus*bc_info[f].vector.dot(normal)
                                             + P_plus*bc_info[f].vector.dot(normal)
                                             - KinematicVector(tau_plus*bc_info[f].vector, job->JOB_TYPE)).dot(normal);


                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == VELOCITY_DENSITY_INLET) {
            //density and velocity given
            //reconstruct traction field at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    p_minus = bc_info[f].vector*bc_info[f].values[0]; //p = u * rho
                    rho_minus = bc_info[f].values[0];
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rhoE_minus*bc_info[f].vector.dot(normal)
                                             + P_minus*bc_info[f].vector.dot(normal)
                                             - KinematicVector(tau_minus*bc_info[f].vector, job->JOB_TYPE)).dot(normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    p_plus = bc_info[f].vector * bc_info[f].values[0]; //p = rho * u
                    rho_plus = bc_info[f].values[0];
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rhoE_plus*bc_info[f].vector.dot(normal)
                                             + P_plus*bc_info[f].vector.dot(normal)
                                             - KinematicVector(tau_plus*bc_info[f].vector, job->JOB_TYPE)).dot(normal);

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == VELOCITY_TEMP_INLET) {
            //velocity and temperature given
            //reconstruct traction field at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = rho_minus * bc_info[f].vector;
                    P_minus = driver->fluid_material->getPressureFromDensityAndTemperature(job, driver, rho_minus, bc_info[f].values[0], n_q(q_list[q]));

                    //rhoE = rho * eps + 0.5*rho*u^2
                    rhoE_minus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver, P_minus, bc_info[f].values[0], n_q(q_list[q]));
                    rhoE_minus += 0.5*rho_minus*bc_info[f].vector.dot(bc_info[f].vector);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rhoE_minus*bc_info[f].vector.dot(normal)
                                             + P_minus*bc_info[f].vector.dot(normal)
                                             - (tau_minus*bc_info[f].vector).dot(normal));

                    //calculate heat flux
                    theta_minus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_minus),
                                                                         driver->fluid_body->p[e_minus],
                                                                         driver->fluid_body->rhoE(e_minus),
                                                                         n_q(q_list[q]));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_minus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_bar, theta_bar, theta_x, n_q(q_list[q]));
                    // add to flux calculation
                    flux += w_q(q_list[q])*heat_flux.dot(normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = rho_plus * bc_info[f].vector;
                    P_plus = driver->fluid_material->getPressureFromDensityAndTemperature(job, driver, rho_plus, bc_info[f].values[0], n_q(q_list[q]));

                    //rhoE = rho * eps + 0.5*rho*u^2
                    rhoE_plus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver, P_plus, bc_info[f].values[0], n_q(q_list[q]));
                    rhoE_plus += 0.5*rho_plus*bc_info[f].vector.dot(bc_info[f].vector);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rhoE_plus*bc_info[f].vector.dot(normal)
                                             + P_plus*bc_info[f].vector.dot(normal)
                                             - (tau_plus*bc_info[f].vector).dot(normal));

                    //calculate heat flux
                    theta_plus = driver->fluid_material->getTemperature(job, driver,
                                                                        driver->fluid_body->rho(e_plus),
                                                                        driver->fluid_body->p[e_plus],
                                                                        driver->fluid_body->rhoE(e_plus),
                                                                        n_q(q_list[q]));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_plus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_bar, theta_bar, theta_x, n_q(q_list[q]));
                    // add to flux calculation
                    flux += w_q(q_list[q])*heat_flux.dot(normal);

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PRESSURE_INLET) {
            //face has prescribed pressure (and temperature)
            //reconstruct velocity field at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                       bc_info[f].values[0],
                                                                                       bc_info[f].values[1],
                                                                                       n_q(q_list[q]));

                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;
                    u_minus = p_minus/rho_minus;

                    //tau_minus = 0
                    P_minus = bc_info[f].values[0];

                    rhoE_minus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                     bc_info[f].values[0],
                                                                                                     bc_info[f].values[1],
                                                                                                     n_q(q_list[q]));
                    rhoE_minus += 0.5*rho_bar*u_minus.dot(u_minus);

                    flux = w_q(q_list[q]) * (rhoE_minus*u_minus.dot(normal)
                                             + P_minus*u_minus.dot(normal));

                    //calculate heat flux
                    theta_minus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_minus),
                                                                         driver->fluid_body->p[e_minus],
                                                                         driver->fluid_body->rhoE(e_minus),
                                                                         n_q(q_list[q]));
                    theta_bar = bc_info[f].values[1];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_minus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_bar, theta_bar, theta_x, n_q(q_list[q]));
                    // add to flux calculation
                    flux += w_q(q_list[q])*heat_flux.dot(normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;
                    u_plus = p_plus/rho_plus;

                    //tau_minus = 0
                    P_plus = bc_info[f].values[0];

                    rhoE_plus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                     bc_info[f].values[0],
                                                                                                     bc_info[f].values[1],
                                                                                                     n_q(q_list[q]));
                    rhoE_plus += 0.5*rho_bar*u_plus.dot(u_plus);

                    flux = w_q(q_list[q]) * (rhoE_plus*u_plus.dot(normal)
                                             + P_plus*u_plus.dot(normal));

                    //calculate heat flux
                    theta_plus = driver->fluid_material->getTemperature(job, driver,
                                                                        driver->fluid_body->rho(e_plus),
                                                                        driver->fluid_body->p[e_plus],
                                                                        driver->fluid_body->rhoE(e_plus),
                                                                        n_q(q_list[q]));
                    theta_bar = bc_info[f].values[1];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_plus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_bar, theta_bar, theta_x, n_q(q_list[q]));
                    // add to flux calculation
                    flux += w_q(q_list[q])*heat_flux.dot(normal);

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PRESSURE_OUTLET || bc_info[f].tag == DAMPED_OUTLET) {
            //face has prescribed pressure (and temperature)
            //reconstruct density and velocity (but adjust if flow direction changes)
            for (int q = 0; q < q_list.size(); q++) {
                rhoE_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                       bc_info[f].values[0],
                                                                                       bc_info[f].values[1],
                                                                                       n_q(q_list[q]));

                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x;
                    u_minus = p_minus / rho_minus;

                    //reconstruct energy for outflow
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x(e_minus).dot(x);

                    //tau_minus = 0
                    P_minus = bc_info[f].values[0];
                    if (bc_info[f].tag == DAMPED_OUTLET){
                        P_minus = (1.0 - damping_coefficient)*P_minus
                                  + damping_coefficient*driver->fluid_material->getPressure(job,driver,
                                                                                            rho_minus,
                                                                                            p_minus,
                                                                                            rhoE_minus,
                                                                                            n_q(q_list[q]));
                    }

                    //calculate energy for inflow
                    rhoE_bar = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                   P_minus,
                                                                                                   bc_info[f].values[1],
                                                                                                   n_q(q_list[q]));
                    rhoE_bar += 0.5*rho_bar*u_minus.dot(u_minus);

                    if (u_minus.dot(normal) > 0) {
                        //flow out of A
                        flux = w_q(q_list[q]) * (rhoE_minus * u_minus.dot(normal)
                                                 + P_minus * u_minus.dot(normal));
                    } else {
                        //flow into A
                        flux = w_q(q_list[q]) * (rhoE_bar * u_minus.dot(normal)
                                                 + P_minus * u_minus.dot(normal));
                    }


                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;

                    //reconstruct energy for outflow
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x(e_plus).dot(x);

                    //tau_minus = 0
                    P_plus = bc_info[f].values[0];
                    if (bc_info[f].tag == DAMPED_OUTLET){
                        P_plus = (1.0 - damping_coefficient)*P_plus
                                  + damping_coefficient*driver->fluid_material->getPressure(job,driver,
                                                                                            rho_plus,
                                                                                            p_plus,
                                                                                            rhoE_plus,
                                                                                            n_q(q_list[q]));
                    }

                    //calculate energy for inflow
                    rhoE_bar = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                   P_plus,
                                                                                                   bc_info[f].values[1],
                                                                                                   n_q(q_list[q]));
                    rhoE_bar += 0.5*rho_bar*u_plus.dot(u_plus);

                    if (u_plus.dot(normal) < 0) {
                        //flow out of B
                        flux = w_q(q_list[q]) * (rhoE_plus * u_plus.dot(normal)
                                                 + P_plus * u_plus.dot(normal));
                    } else {
                        //flow into A
                        flux = w_q(q_list[q]) * (rhoE_bar * u_plus.dot(normal)
                                                 + P_plus * u_plus.dot(normal));
                    }

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == ADIABATIC_WALL || bc_info[f].tag == SYMMETRIC_WALL){
            //nothing to do right now
        } else if (bc_info[f].tag == THERMAL_WALL){
            //face has prescribed temperature
            //calculate heat flux at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);

                    //calculate heat flux
                    theta_minus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_minus),
                                                                         driver->fluid_body->p[e_minus],
                                                                         driver->fluid_body->rhoE(e_minus),
                                                                         n_q(q_list[q]));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_minus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_bar, theta_bar, theta_x, n_q(q_list[q]));
                    // add to flux calculation
                    flux += w_q(q_list[q])*heat_flux.dot(normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;
                    u_plus = p_plus/rho_plus;

                    //tau_minus = 0
                    P_plus = bc_info[f].values[0];

                    rhoE_plus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                    bc_info[f].values[0],
                                                                                                    bc_info[f].values[1],
                                                                                                    n_q(q_list[q]));
                    rhoE_plus += 0.5*rho_bar*u_plus.dot(u_plus);

                    flux = w_q(q_list[q]) * (rhoE_plus*u_plus.dot(normal)
                                             + P_plus*u_plus.dot(normal));

                    //calculate heat flux
                    theta_plus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_plus),
                                                                         driver->fluid_body->p[e_plus],
                                                                         driver->fluid_body->rhoE(e_plus),
                                                                         n_q(q_list[q]));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_plus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_bar, theta_bar, theta_x, n_q(q_list[q]));
                    // add to flux calculation
                    flux += w_q(q_list[q])*heat_flux.dot(normal);

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == SUPERSONIC_INLET) {
            //density and velocity given (and temperature if thermal)
            //don't worry about heat flux for now
            for (int q = 0; q < q_list.size(); q++) {
                //neglect shear stress at inlet
                P_minus = driver->fluid_material->getPressureFromDensityAndTemperature(job, driver,
                                                                                       bc_info[f].values[0],
                                                                                       bc_info[f].values[1],
                                                                                       n_q(q_list[q]));

                rhoE_bar = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                               P_minus,
                                                                                               bc_info[f].values[1],
                                                                                               n_q(q_list[q]));
                rhoE_bar += 0.5*bc_info[f].values[0]*bc_info[f].vector.dot(bc_info[f].vector);


                flux = w_q(q_list[q]) * (rhoE_bar * bc_info[f].vector.dot(normal)
                                         + P_minus * bc_info[f].vector.dot(normal));

                //check which side of face has valid element
                if (e_minus > -1) {
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == SUPERSONIC_OUTLET) {
            //reconstruct flux (CAREFUL! back-flow unstable)
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p(e_minus) + driver->fluid_body->p_x[e_minus] * x;
                    u_minus = p_minus/rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rhoE_minus * u_minus.dot(normal)
                                             + P_minus * u_minus.dot(normal));
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p(e_plus) + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus/rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rhoE_plus * u_plus.dot(normal)
                                             + P_plus * u_plus.dot(normal));
                    result(e_plus) += flux;
                }
            }
        } else {
            std::cerr << "ERROR! FVMGridBase does not have flux defined for bc tag " << bc_info[f].tag << "!"
                      << std::endl;
            //don't exit, without flux defined, essentially a wall
        }
    }
    return result;
}
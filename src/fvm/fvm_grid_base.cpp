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
//function to generate mapping matrices
void FVMGridBase::generateMappings(Job* job, FiniteVolumeDriver* driver){
    //temporary containers
    std::vector<int> nvec(0);
    std::vector<double> valvec(0);
    KinematicVectorArray gradvec(0,job->JOB_TYPE);
    KinematicVector tmpGrad(job->JOB_TYPE);
    KinematicVector tmpVec(job->JOB_TYPE);

    //generate Q matrix [N_i(x_q)]
    Q = MPMScalarSparseMatrix(job->grid->node_count, int_quad_count + ext_quad_count);

    //loop over quadartature points
    for (int q = 0; q < int_quad_count+ext_quad_count; q++) {
        nvec.resize(0);
        valvec.resize(0);

        //for each quadrature point, get MPM basis function values
        tmpVec = x_q[q];

        if (job->grid->whichElement(job, tmpVec) > -1){
            //do nothing
        } else {
            tmpVec *= (1.0 - 1e-14); //adjust tmpVec slightly
        }

        job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);

        //add basis functions to sparse matrix map
        for (int j = 0; j < nvec.size(); j++) {
            Q.push_back(nvec[j], q, valvec[j]); //node, point, value
        }
    }

    //generate gradQ matrix [gradN_i(x_q)]
    gradQ = KinematicVectorSparseMatrix(job->grid->node_count, int_quad_count + ext_quad_count, job->JOB_TYPE);

    //loop over quadrature points
    for (int q = 0; q < int_quad_count+ext_quad_count; q++) {
        nvec.resize(0);
        gradvec.resize(0);

        //for each quadrature point, get MPM basis function gradients
        tmpVec = x_q[q];

        if (job->grid->whichElement(job, tmpVec) > -1){
            //do nothing
        } else {
            tmpVec *= (1.0 - 1e-14); //adjust tmpVec slightly
        }

        job->grid->evaluateBasisFnGradient(job, tmpVec, nvec, gradvec);
        for (int j = 0; j < nvec.size(); j++) {
            gradQ.push_back(nvec[j], q, gradvec[j]); //node, point, value
        }
    }

    //generate M_ei matrix [int{M_e(x)N_i(x)]
    M = MPMScalarSparseMatrix(element_count, job->grid->node_count);

    //loop over elements
    for (int e=0; e<element_count; e++){
        for (int q=0; q<qpe; q++){
            nvec.resize(0);
            valvec.resize(0);

            //for each interior quadrature point, get MPM basis function value
            tmpVec = x_q[e*qpe + q];

            if (job->grid->whichElement(job, tmpVec) > -1){
                //do nothing
            } else {
                tmpVec *= (1.0 - 1e-14); //adjust tmpVec slightly
            }

            job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);

            //add all weighted MPM function values to matrix
            for (int j = 0; j < nvec.size(); j++) {
                M.push_back(e, nvec[j], w_q(e*qpe + q)*valvec[j]); //element, node, quad weight
            }
        }
    }

    return;
}

/*----------------------------------------------------------------------------*/
//functions for solving system of equations

void FVMGridBase::mapMixturePropertiesToQuadraturePoints(Job* job, FiniteVolumeDriver* driver){
    //need to map porosity and solid velocity to quadrature points for flux calculations

    //ask material to update body porosity field
    if (driver->fluid_material->calculatePorosity(job, driver) == 1){
        //successful calculation of porosity
        n_q = Q.operate(driver->fluid_body->n, MPMSparseMatrixBase::TRANSPOSED);          //n_q = Q_iq * n_i
        gradn_q = gradQ.operate(driver->fluid_body->n, MPMSparseMatrixBase::TRANSPOSED);  //gradn_q = gradQ_iq * n_i

        //integrate porosity field over element volumes
        driver->fluid_body->n_e = driver->fluid_grid->M * driver->fluid_body->n;
        for (int e=0; e<driver->fluid_grid->element_count; e++) {
            //average porosity over element volume
            driver->fluid_body->n_e(e) /= driver->fluid_grid->getElementVolume(e);
        }
    } else {
        //no porosity calculation performed
        //n_q.setConstant(1.0);
        //gradn_q.setZero();
    }

    //ask material to update solid velocity field
    if (driver->fluid_material->updateSolidPhaseVelocity(job, driver) == 1){
        //successful update of solid velocity
        v_sq = Q.operate(driver->fluid_body->v_s, MPMSparseMatrixBase::TRANSPOSED);       //v_sq = Q_iq * v_si
    } else {
        //no porosity calculation performed
        //v_sq.setZero();
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
    double n_plus, n_minus, n_bar;

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

        e_minus = face_elements[f][0];
        e_plus = face_elements[f][1];

        //flux calculation depends on whether face is on boundary
        if (bc_info[f].tag == -1 || (bc_info[f].tag == PERIODIC && face_elements[f][0] > -1)) {
            //face is interior to domain; no BCs (or periodic)

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x) / rho_plus;
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) /
                            (std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = std::sqrt(rho_plus * rho_minus);
                    rhoE_minus = driver->fluid_body->rhoE(e_minus);
                    //c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                    //                                            n_q(q_list[q]));
                    c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                                                                std::sqrt(n_minus * n_plus));

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_minus); //n_q(q_list[q]));
                    H_minus = E_minus + P_minus/rho_minus;

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    E_plus = rhoE_plus/rho_plus;
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_plus); //n_q(q_list[q]));
                    H_plus = E_plus + P_plus/rho_plus;

                    //approximate Roe advective rate
                    w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus))/2.0;
                    rho_bar = w_1_bar*w_1_bar;
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) / (2.0 * w_1_bar);
                    u_bar_squared = u_bar.dot(u_bar);
                    H_bar = (std::sqrt(rho_plus) * H_plus + std::sqrt(rho_minus) * H_minus) / (2.0 * w_1_bar);
                    n_bar = (std::sqrt(rho_plus) * n_plus + std::sqrt(rho_minus) * n_minus) / (2.0 * w_1_bar);
                    c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                                                                n_bar);  //n_q(q_list[q]));

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
                                                                                       n_q(q_list[q]));         //use quadrature porosity here
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
                                                                                       n_q(q_list[q]));      //use quadrature porosity
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
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
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
        } else if (bc_info[f].tag == PERIODIC){
            //if you get here, this face is ill-defined. do nothing
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
    double n_plus, n_minus, n_bar;
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

        e_minus = face_elements[f][0];
        e_plus = face_elements[f][1];

        //flux calculation depends on whether face is on boundary
        if (bc_info[f].tag == -1 || (bc_info[f].tag == PERIODIC && face_elements[f][0] > -1)) {
            //face is interior to domain; no BCs (or periodic)

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) /
                            (std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = std::sqrt(rho_plus * rho_minus);
                    rhoE_minus = driver->fluid_body->rhoE(e_minus);
                    n_bar = std::sqrt(n_plus*n_minus);
                    c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                                                                n_bar); //n_q(q_list[q]));

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
                    //P_plus = c*c*(rho_plus - rho_bar) + driver->fluid_material->getPressure(job, driver, rho_bar, p_minus, rhoE_minus, n_q(q_list[q]));
                    //P_minus = c*c*(rho_minus - rho_plus) + P_plus;
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_bar, p_plus, rhoE_minus, n_plus);
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_bar, p_minus, rhoE_minus, n_minus);

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus,
                                                                  n_minus); //n_q(q_list[q]));
                    H_minus = E_minus + P_minus / rho_minus;

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus)  + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    E_plus = rhoE_plus / rho_plus;
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus,
                                                                 n_plus); //n_q(q_list[q]));
                    H_plus = E_plus + P_plus / rho_plus;

                    //approximate Roe advective rate
                    w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
                    rho_bar = w_1_bar * w_1_bar;
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) / (2.0 * w_1_bar);
                    u_bar_squared = u_bar.dot(u_bar);
                    H_bar = (std::sqrt(rho_plus) * H_plus + std::sqrt(rho_minus) * H_minus) / (2.0 * w_1_bar);
                    n_bar = (std::sqrt(rho_plus) * n_plus + std::sqrt(rho_minus) * n_minus) / (2.0 * w_1_bar);
                    c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                                                                n_bar); //n_q(q_list[q]));

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_minus); //n_q(q_list[q]));

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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_plus); //n_q(q_list[q]));

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_minus); //n_q(q_list[q]));

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
                    n_plus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_plus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_plus); //n_q(q_list[q]));

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);
                    P_minus = driver->fluid_material->getPressureFromDensityAndTemperature(job,
                                                                                           driver,
                                                                                           rho_minus,
                                                                                           bc_info[f].values[0],
                                                                                           n_minus); //n_q(q_list[q]));

                    //rhoE = rho * eps + 0.5*rho*u^2
                    rhoE_minus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job,
                                                                                                     driver,
                                                                                                     P_minus,
                                                                                                     bc_info[f].values[0],
                                                                                                     n_minus); //n_q(q_list[q]));
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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);
                    P_plus = driver->fluid_material->getPressureFromDensityAndTemperature(job,
                                                                                          driver,
                                                                                          rho_plus,
                                                                                          bc_info[f].values[0],
                                                                                          n_plus); //n_q(q_list[q]));

                    //rhoE = rho * eps + 0.5*rho*u^2
                    rhoE_plus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job,
                                                                                                    driver,
                                                                                                    P_plus,
                                                                                                    bc_info[f].values[0],
                                                                                                    n_plus); //n_q(q_list[q]));
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
                                                                                       n_q(q_list[q]));         //ok to use quad value here

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
                                                                                       n_q(q_list[q]));         //ok to use quad value here

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
                        n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

                        P_minus = (1.0 - damping_coefficient)*P_minus
                                 + damping_coefficient*driver->fluid_material->getPressure(job,driver,
                                                                                           rho_minus,
                                                                                           p_minus,
                                                                                           rhoE_minus, n_minus); //n_q(q_list[q])); //ok to use quad value here
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
                        n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);

                        P_plus = (1.0 - damping_coefficient)*P_plus
                                 + damping_coefficient*driver->fluid_material->getPressure(job,driver,
                                                                                           rho_plus,
                                                                                           p_plus,
                                                                                           rhoE_plus, n_plus); //n_q(q_list[q])); //ok to use porosity at quad point
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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_minus); //n_q(q_list[q]));

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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_plus); // n_q(q_list[q]));

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_minus); // n_q(q_list[q]));

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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);

                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_plus); //n_q(q_list[q]));

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
                                                                                       1.0); //n_q(q_list[q]));

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_minus); //n_q(q_list[q]));

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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_plus); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (p_plus * u_plus.dot(normal)
                                             + P_plus * normal);
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PERIODIC){
            //if you get here, this face is ill-defined. do nothing
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
    double n_plus, n_minus, n_bar;
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

        //get oriented faces
        e_minus = face_elements[f][0];
        e_plus = face_elements[f][1];

        //flux calculation depends on whether face is on boundary
        if (bc_info[f].tag == -1 || (bc_info[f].tag == PERIODIC && face_elements[f][0] > -1)) {
            //face is interior to domain; no BCs (or periodic)

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);
                    E_minus = rhoE_minus / rho_minus;
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus,
                                                                  n_minus); //n_q(q_list[q]));
                    H_minus = E_minus + P_minus / rho_minus;
                    tau_minus = driver->fluid_material->getShearStress(job, driver, L[e_minus], rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus)  + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);
                    E_plus = rhoE_plus / rho_plus;
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus,
                                                                 n_plus); //n_q(q_list[q]));
                    H_plus = E_plus + P_plus / rho_plus;
                    tau_plus = driver->fluid_material->getShearStress(job, driver, L[e_plus], rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));

                    //approximate Roe advective rate
                    w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
                    rho_bar = w_1_bar * w_1_bar;
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) / (2.0 * w_1_bar);
                    u_bar_squared = u_bar.dot(u_bar);
                    H_bar = (std::sqrt(rho_plus) * H_plus + std::sqrt(rho_minus) * H_minus) / (2.0 * w_1_bar);
                    n_bar = (std::sqrt(rho_plus) * n_plus + std::sqrt(rho_minus) * n_minus) / (2.0 * w_1_bar);
                    c = driver->fluid_material->getSpeedOfSound(job, driver, rho_bar, rho_bar * u_bar, rhoE_minus,
                                                                n_bar); //n_q(q_list[q]));

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_minus); //n_q(q_list[q]));

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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_plus); //n_q(q_list[q]));

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_minus, p_minus, rhoE_minus, n_q(q_list[q]));
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_minus); //n_q(q_list[q]));

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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, L_tmp, rho_plus, p_plus, rhoE_plus, n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_plus); //n_q(q_list[q]));

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);
                    P_minus = driver->fluid_material->getPressureFromDensityAndTemperature(job,
                                                                                           driver,
                                                                                           rho_minus,
                                                                                           bc_info[f].values[0],
                                                                                           n_minus); //n_q(q_list[q]));

                    //rhoE = rho * eps + 0.5*rho*u^2
                    rhoE_minus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job,
                                                                                                     driver,
                                                                                                     P_minus,
                                                                                                     bc_info[f].values[0],
                                                                                                     n_minus); //n_q(q_list[q]));

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
                                                                         driver->fluid_body->n_e(e_minus));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_minus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_minus, theta_bar, theta_x, n_q(q_list[q]));
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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);
                    P_plus = driver->fluid_material->getPressureFromDensityAndTemperature(job,
                                                                                          driver,
                                                                                          rho_plus,
                                                                                          bc_info[f].values[0],
                                                                                          n_plus); //n_q(q_list[q]));

                    //rhoE = rho * eps + 0.5*rho*u^2
                    rhoE_plus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job,
                                                                                                    driver,
                                                                                                    P_plus,
                                                                                                    bc_info[f].values[0],
                                                                                                    n_plus); //n_q(q_list[q]));
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
                                                                        driver->fluid_body->n_e(e_plus));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_plus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_plus, theta_bar, theta_x, n_q(q_list[q]));
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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

                    //tau_minus = 0
                    P_minus = bc_info[f].values[0];

                    rhoE_minus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                     bc_info[f].values[0],
                                                                                                     bc_info[f].values[1],
                                                                                                     n_minus); //n_q(q_list[q]));
                    rhoE_minus += 0.5*rho_bar*u_minus.dot(u_minus);

                    flux = w_q(q_list[q]) * (rhoE_minus*u_minus.dot(normal)
                                             + P_minus*u_minus.dot(normal));

                    //calculate heat flux
                    theta_minus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_minus),
                                                                         driver->fluid_body->p[e_minus],
                                                                         driver->fluid_body->rhoE(e_minus),
                                                                         driver->fluid_body->n_e(e_minus));
                    theta_bar = bc_info[f].values[1];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_minus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_minus, theta_bar, theta_x, n_q(q_list[q]));
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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);

                    //tau_minus = 0
                    P_plus = bc_info[f].values[0];

                    rhoE_plus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                     bc_info[f].values[0],
                                                                                                     bc_info[f].values[1],
                                                                                                     n_plus); //n_q(q_list[q]));
                    rhoE_plus += 0.5*rho_bar*u_plus.dot(u_plus);

                    flux = w_q(q_list[q]) * (rhoE_plus*u_plus.dot(normal)
                                             + P_plus*u_plus.dot(normal));

                    //calculate heat flux
                    theta_plus = driver->fluid_material->getTemperature(job, driver,
                                                                        driver->fluid_body->rho(e_plus),
                                                                        driver->fluid_body->p[e_plus],
                                                                        driver->fluid_body->rhoE(e_plus),
                                                                        driver->fluid_body->n_e(e_plus));
                    theta_bar = bc_info[f].values[1];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_plus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_plus, theta_bar, theta_x, n_q(q_list[q]));
                    // add to flux calculation
                    flux += w_q(q_list[q])*heat_flux.dot(normal);

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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);

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
                                                                                            n_minus); //n_q(q_list[q]));
                    }

                    //calculate energy for inflow
                    rhoE_bar = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                   P_minus,
                                                                                                   bc_info[f].values[1],
                                                                                                   n_minus); //n_q(q_list[q]));
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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);

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
                                                                                            n_plus); //n_q(q_list[q]));
                    }

                    //calculate energy for inflow
                    rhoE_bar = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                   P_plus,
                                                                                                   bc_info[f].values[1],
                                                                                                   n_plus); //n_q(q_list[q]));
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

                    //calculate heat flux
                    theta_minus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_minus),
                                                                         driver->fluid_body->p[e_minus],
                                                                         driver->fluid_body->rhoE(e_minus),
                                                                         driver->fluid_body->n_e(e_minus));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_minus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_minus, theta_bar, theta_x, n_q(q_list[q]));
                    // add to flux calculation
                    flux = w_q(q_list[q])*heat_flux.dot(normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);

                    //calculate heat flux
                    theta_plus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_plus),
                                                                         driver->fluid_body->p[e_plus],
                                                                         driver->fluid_body->rhoE(e_plus),
                                                                         driver->fluid_body->n_e(e_plus));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_plus)*x/(x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, rho_plus, theta_bar, theta_x, n_q(q_list[q]));
                    // add to flux calculation
                    flux = w_q(q_list[q])*heat_flux.dot(normal);

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
                                                                                       1.0);

                rhoE_bar = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                               P_minus,
                                                                                               bc_info[f].values[1],
                                                                                               1.0);
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
                    n_minus = driver->fluid_body->n_e(e_minus) + driver->fluid_body->n_e_x[e_minus].dot(x);
                    P_minus = driver->fluid_material->getPressure(job, driver, rho_minus, p_minus, rhoE_minus, n_minus); //n_q(q_list[q]));

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
                    n_plus = driver->fluid_body->n_e(e_plus) + driver->fluid_body->n_e_x[e_plus].dot(x);
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_plus); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (rhoE_plus * u_plus.dot(normal)
                                             + P_plus * u_plus.dot(normal));
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PERIODIC){
            //if you get here, this face is ill-defined. do nothing
        } else {
            std::cerr << "ERROR! FVMGridBase does not have flux defined for bc tag " << bc_info[f].tag << "!"
                      << std::endl;
            //don't exit, without flux defined, essentially a wall
        }
    }
    return result;
}

/*----------------------------------------------------------------------------*/
//function to calculate interphase force
KinematicVectorArray FVMGridBase::calculateInterphaseForces(Job* job, FiniteVolumeDriver* driver){

    //result defines force on solid by fluid
    KinematicVectorArray result = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray tmp_result = result;

    //declare temporary containers for quadrature point fields
    Eigen::VectorXd P_q = Eigen::VectorXd(int_quad_count + ext_quad_count); //pressure
    KinematicVectorArray f_d = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE); //drag
    KinematicVectorArray n_b = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE); //outward pointing normals
    KinematicVectorArray tmpVecArray = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE);
    Eigen::VectorXd tmpArray = Eigen::VectorXd(int_quad_count + ext_quad_count);

    P_q.setZero();
    f_d.setZero();
    tmpVecArray.setZero();
    tmpArray.setZero();

    double rho, rhoE;
    KinematicVector p, normal;

    //loop over elements
    for (int e=0; e<element_count; e++){
        for (int q=0; q<qpe; q++){
            //approximate local density, momentum, and energy
            rho = driver->fluid_body->rho(e) + driver->fluid_body->rho_x[e].dot(x_q[e*qpe + q] - x_e[e]);
            rhoE = driver->fluid_body->rhoE(e) + driver->fluid_body->rhoE_x[e].dot(x_q[e*qpe + q] - x_e[e]);
            p = driver->fluid_body->p[e] + driver->fluid_body->p_x[e]*(x_q[e*qpe + q] - x_e[e]);

            //fill in drag force
            f_d[e * qpe + q] = driver->fluid_material->getInterphaseDrag(job,
                                                                         driver,
                                                                         rho,
                                                                         (p/rho),
                                                                         v_sq[e*qpe + q],
                                                                         n_q(e*qpe + q));

            //fill in pressure
            P_q(e * qpe + q) = driver->fluid_material->getPressure(job,
                                                                   driver,
                                                                   rho,
                                                                   p,
                                                                   rhoE,
                                                                   n_q(e*qpe + q));

            //fill in tmpVecArray and tmpArray
            //tmpVecArray[e*qpe + q] = (-f_d[e*qpe + q] - P_q(e*qpe + q)*gradn_q[e*qpe+q]) * getQuadratureWeight(e*qpe + q);
            //tmpArray(e*qpe + q) = P_q(e*qpe + q) * (1.0 - n_q(e*qpe+q)) * getQuadratureWeight(e*qpe + q);
            tmpVecArray[e*qpe + q] = -f_d[e*qpe + q] * getQuadratureWeight(e*qpe + q);
            tmpArray(e*qpe + q) = P_q(e*qpe + q) * getQuadratureWeight(e*qpe + q);
        }
    }

    //integrate expression (pg 117, nb 6)
    result = Q*tmpVecArray;
    tmp_result = gradQ*tmpArray;
    for (int i=0; i<result.size(); i++){
        result[i] += (1.0 - driver->fluid_body->n(i))*tmp_result(i);
    }

    //zero temporary vectors
    tmp_result.setZero();
    tmpVecArray.setZero();

    int tmp_e = -1;
    int tmp_q = -1;
    //loop over faces
    for (int f=0; f<face_count; f++){
        for (int q=0; q<qpf; q++){
            //only add contributions from domain boundaries
            if (q_b[int_quad_count + f*qpf + q]){
                //get normal
                normal = getFaceNormal(job, f);
                if (face_elements[f][0] > -1){
                    //A is interior to domain
                    tmp_e = face_elements[f][0];
                } else {
                    //B is interior to domain
                    normal *= -1.0;
                    tmp_e = face_elements[f][1];
                }

                //approximate local density, momentum, and energy
                tmp_q = int_quad_count + f*qpf + q;
                rho = driver->fluid_body->rho(tmp_e) + driver->fluid_body->rho_x[tmp_e].dot(x_q[tmp_q] - x_e[tmp_e]);
                rhoE = driver->fluid_body->rhoE(tmp_e) + driver->fluid_body->rhoE_x[tmp_e].dot(x_q[tmp_q] - x_e[tmp_e]);
                p = driver->fluid_body->p[tmp_e] + driver->fluid_body->p_x[tmp_e]*(x_q[tmp_q] - x_e[tmp_e]);

                //fill in pressure
                P_q(int_quad_count + f * qpf + q) = driver->fluid_material->getPressure(job,
                                                                                       driver,
                                                                                       rho,
                                                                                       p,
                                                                                       rhoE,
                                                                                       n_q(tmp_q));

                //fill in tmpArray
                //tmpVecArray[tmp_q] = -P_q(tmp_q) * (1.0 - n_q(tmp_q)) * normal * getQuadratureWeight(tmp_q);
                tmpVecArray[tmp_q] = -P_q(tmp_q) * normal * getQuadratureWeight(tmp_q);
            }
        }
    }

    //integrate expression (pg 117, nb 6)
    tmp_result = Q*tmpVecArray;
    for (int i=0; i<result.size(); i++){
        result[i] += (1.0 - driver->fluid_body->n(i))*tmp_result(i);
    }

    //divide result by associated nodal volumes
    for (int i=0; i<job->grid->node_count; i++){
        result(i) /= job->grid->nodeVolume(job, i);
    }

    return result;
}
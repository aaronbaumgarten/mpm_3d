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

#include "objects/solvers/solvers.hpp"

static const bool USE_OLD_DRAG_METHOD = false;

/*----------------------------------------------------------------------------*/
void FVMGridBase::init(Job* job, FiniteVolumeDriver* driver){
    //do some basic set-up of parameters

    //get pointer to job->threadPool
    jobThreadPool = &job->threadPool;

    //get number of threads
    num_threads = job->thread_count;

    //initialize memory unit
    memoryUnits.push_back(parallelMemoryUnit());

    //loop over str-props and assign relevant flags
    std::vector<std::string> options = {"USE_REDUCED_QUADRATURE",
                                        "USE_LOCAL_GRADIENT_CORRECTION",
                                        "USE_CUSTOM_DAMPING",
                                        "USE_GAUSS_GREEN",
                                        "USE_LOCAL_POROSITY_CORRECTION",
                                        "USE_DECAYING_DAMPING",
                                        "USE_ENHANCED_M_MATRIX",
                                        "USE_ENHANCED_QUADRATURE",
                                        "USE_HYDROSTATIC_CORRECTION",
                                        "USE_MACH_LIMITER",
                                        "USE_CORRECT_INTERPHASE_ENERGY_FLUX"};
    for (int i=0; i<str_props.size(); i++){
        switch (Parser::findStringID(options, str_props[i])){
            case 0:
                //USE_REDUCED_QUADRATURE
                USE_REDUCED_QUADRATURE = true;
                std::cout << "FiniteVolumeGrid using reduced quadrature." << std::endl;
                break;
            case 1:
                //USE_LOCAL_GRADIENT_CORRECTION
                USE_LOCAL_GRADIENT_CORRECTION = true;
                std::cout << "FiniteVolumeGrid using local gradient correction for shear stresses." << std::endl;
                break;
            case 2:
                //USE_CUSTOM_DAMPING
                damping_coefficient = fp64_props[fp64_props.size()-1];
                std::cout << "FiniteVolumeGrid using custom damping coefficient: " << damping_coefficient << std::endl;
                break;
            case 3:
                //USE_GAUSS_GREEN_POROSITY
                USE_GAUSS_GREEN = true;
                std::cout << "FiniteVolumeGrid using Gauss-Green reconstruction of MPM gradients." << std::endl;
                break;
            case 4:
                //USE_LOCAL_POROSITY_CORRECTION
                USE_LOCAL_POROSITY_CORRECTION = true;
                std::cout << "FiniteVolumeGrid using local porosity correction during flux reconstruction." << std::endl;
                break;
            case 5:
                //USE_DECAYING_DAMPING
                USE_DECAYING_DAMPING = true;
                damping_at_t0 = fp64_props[fp64_props.size()-2];
                damping_time_cst = fp64_props[fp64_props.size()-1];
                std::cout << "FiniteVolumeGrid using decaying damping coefficient." << std::endl;
                //do not use with USE_CUSTOM_DAMPING flag!!!
                break;
            case 6:
                //USE_ENHANCED_M_MATRIX
                USE_ENHANCED_M_MATRIX = true;
                hx_enhanced_quad = fp64_props[fp64_props.size()-1];
                std::cout << "FiniteVolumeGrid using enhanced M matrix calculation." << std::endl;
                //do not use with USE_CUSTOM_DAMPING flag!!!
                break;
            case 7:
                //USE_ENHANCED_QUADRATURE
                USE_ENHANCED_QUADRATURE = true;
                lqpe_enhanced = int_props[int_props.size()-1];
                std::cout << "FiniteVolumeGrid using enhanced quadrature calculation." << std::endl;
                //do not use with USE_CUSTOM_DAMPING flag!!!
                break;
            case 8:
                //USE_HYDROSTATIC_CORRECTION
                USE_HYDROSTATIC_CORRECTION = true;
                std::cout << "FiniteVolumeGrid using hydrostatically corrected pressure BCs." << std::endl;
                //do not use with USE_CUSTOM_DAMPING flag!!!
                break;
            case 9:
                //USE_MACH_LIMIER
                USE_MACH_LIMITER = true;
                std::cout << "FiniteVolumeGrid using Mach limiter for stagnation BCs." << std::endl;
                //do not use with USE_CUSTOM_DAMPING flag!!!
                break;
            case 10:
                //USE_CORRECT_INTERPHASE_ENERGY_FLUX
                USE_CORRECT_INTERPHASE_ENERGY_FLUX = true;
                std::cout << "FiniteVolumeGrid using correct inter-phaseergy flux for fluid (from thesis)." << std::endl;
                break;
            default:
                //do nothing
                break;
        }
    }

    return;
}

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

    //generate gradQ matrix [gradN_i(x_q)]
    gradQ = KinematicVectorSparseMatrix(job->grid->node_count, int_quad_count + ext_quad_count, job->JOB_TYPE);

    //loop over quadartature points
    for (int q = 0; q < int_quad_count+ext_quad_count; q++) {
        nvec.resize(0);
        valvec.resize(0);

        //for each quadrature point, get MPM basis function values
        tmpVec = x_q[q];

        if (job->grid->whichElement(job, tmpVec) > -1){
            //do nothing
        } else {
            tmpVec *= (1.0 - 1e-12); //adjust tmpVec slightly
        }

        job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);

        //add basis functions to sparse matrix map
        for (int j = 0; j < nvec.size(); j++) {
            Q.push_back(nvec[j], q, valvec[j]); //node, point, value
        }
    }



    int e_plus, e_minus, tmp_q, tmp_qq;
    KinematicVector normal;
    if (USE_GAUSS_GREEN){
        //loop over faces and estimate gradN_i(x_q) using N_i(x_q)
        for (int f=0; f<face_count; f++){
            e_minus = face_elements[f][0];
            e_plus = face_elements[f][1];
            normal = getFaceNormal(job, f);
            //loop over boundary quadrature points
            for (int q=0; q<qpf; q++){
                nvec.resize(0);
                valvec.resize(0);

                //get quadrature point index
                tmp_q = int_quad_count + f*qpf + q;

                //for each quadrature point, get MPM basis function values
                tmpVec = x_q[tmp_q];

                if (job->grid->whichElement(job, tmpVec) > -1){
                    //do nothing
                } else {
                    tmpVec *= (1.0 - 1e-14); //adjust tmpVec slightly
                }

                job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);

                //add boundary influence on interior quad points (let gradients be uniform in cells)
                //int(grad(n) dv) = int(n * normal da)
                if (e_minus > -1){
                    for (int qq = 0; qq < qpe; qq++){
                        tmp_qq = e_minus * qpe + qq;
                        for (int j=0; j<nvec.size(); j++){
                            gradQ.push_back(nvec[j], tmp_qq, w_q(tmp_q)*valvec[j]*normal/getElementVolume(e_minus));
                        }
                    }
                }

                if (e_plus > -1){
                    for (int qq = 0; qq < qpe; qq++){
                        tmp_qq = e_plus * qpe + qq;
                        for (int j=0; j<nvec.size(); j++){
                            gradQ.push_back(nvec[j], tmp_qq, -w_q(tmp_q)*valvec[j]*normal/getElementVolume(e_plus));
                        }
                    }
                }

            }
        }
    } else {
        //loop over quadrature points
        for (int q = 0; q < int_quad_count + ext_quad_count; q++) {
            nvec.resize(0);
            gradvec.resize(0);

            //for each quadrature point, get MPM basis function gradients
            tmpVec = x_q[q];

            if (job->grid->whichElement(job, tmpVec) > -1) {
                //do nothing
            } else {
                tmpVec *= (1.0 - 1e-14); //adjust tmpVec slightly
            }

            job->grid->evaluateBasisFnGradient(job, tmpVec, nvec, gradvec);
            for (int j = 0; j < nvec.size(); j++) {
                gradQ.push_back(nvec[j], q, gradvec[j]); //node, point, value
            }
        }
    }

    //generate M_ei matrix [int{M_e(x)N_i(x)]
    M = MPMScalarSparseMatrix(element_count, job->grid->node_count);

    //loop over elements
    if (!USE_ENHANCED_M_MATRIX) {
        for (int e = 0; e < element_count; e++) {
            for (int q = 0; q < qpe; q++) {
                nvec.resize(0);
                valvec.resize(0);

                //for each interior quadrature point, get MPM basis function value
                tmpVec = x_q[e * qpe + q];

                if (job->grid->whichElement(job, tmpVec) > -1) {
                    //do nothing
                } else {
                    tmpVec *= (1.0 - 1e-14); //adjust tmpVec slightly
                }

                job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);

                //add all weighted MPM function values to matrix
                for (int j = 0; j < nvec.size(); j++) {
                    M.push_back(e, nvec[j], w_q(e * qpe + q) * valvec[j]); //element, node, quad weight
                }
            }
        }
    } else {
        //use enhance mapping calculation (for h_i < h_e)
        //create point cartesian point cloud using hx as length scale
        std::vector<double> min_max = getBoundingBox();
        std::vector<int> Nx = std::vector<int>(GRID_DIM);
        for (int i=0; i<GRID_DIM; i++){
            Nx[i] = (int)std::floor((min_max[2*i+1] - min_max[2*i])/hx_enhanced_quad);
        }

        //loop over points in cloud and get cross-contributions to M matrix
        double w = 1;
        for (int i=0; i<GRID_DIM; i++){
            w *= hx_enhanced_quad;
        }
        int e;
        KinematicVector x = KinematicVector(job->JOB_TYPE);
        for (int i=0; i<Nx[0]; i++){
            x[0] = hx_enhanced_quad*(i + 0.5);

            if (GRID_DIM > 1){
                for (int j=0; j<Nx[1]; j++){
                    x[1] = hx_enhanced_quad*(j + 0.5);

                    if (GRID_DIM > 2){
                        for (int k=0; k<Nx[2]; k++){
                            x[2] = hx_enhanced_quad*(k + 0.5);

                            e = whichElement(job, x);

                            nvec.resize(0);
                            valvec.resize(0);

                            if (job->grid->whichElement(job, x) > -1) {
                                //do nothing
                            } else {
                                x *= (1.0 - 1e-14); //adjust tmpVec slightly
                            }

                            job->grid->evaluateBasisFnValue(job, x, nvec, valvec);

                            //add all weighted MPM function values to matrix
                            for (int jj = 0; jj < nvec.size(); jj++) {
                                M.push_back(e, nvec[jj], w * valvec[jj]); //element, node, quad weight
                            }
                        }
                    } else {
                        e = whichElement(job, x);

                        nvec.resize(0);
                        valvec.resize(0);

                        if (job->grid->whichElement(job, x) > -1) {
                            //do nothing
                        } else {
                            x *= (1.0 - 1e-14); //adjust tmpVec slightly
                        }

                        job->grid->evaluateBasisFnValue(job, x, nvec, valvec);

                        //add all weighted MPM function values to matrix
                        for (int jj = 0; jj < nvec.size(); jj++) {
                            M.push_back(e, nvec[jj], w * valvec[jj]); //element, node, quad weight
                        }
                    }
                }
            } else {
                e = whichElement(job, x);

                nvec.resize(0);
                valvec.resize(0);

                if (job->grid->whichElement(job, x) > -1) {
                    //do nothing
                } else {
                    x *= (1.0 - 1e-14); //adjust tmpVec slightly
                }

                job->grid->evaluateBasisFnValue(job, x, nvec, valvec);

                //add all weighted MPM function values to matrix
                for (int jj = 0; jj < nvec.size(); jj++) {
                    M.push_back(e, nvec[jj], w * valvec[jj]); //element, node, quad weight
                }
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
        if (num_threads > 1){
            parallelMultiply(Q, driver->fluid_body->n, n_q, MPMSparseMatrixBase::TRANSPOSED, true);
            parallelMultiply(gradQ, driver->fluid_body->n, gradn_q, MPMSparseMatrixBase::TRANSPOSED, true);
        } else {
            n_q = Q.operate(driver->fluid_body->n, MPMSparseMatrixBase::TRANSPOSED);          //n_q = Q_iq * n_i
            gradn_q = gradQ.operate(driver->fluid_body->n, MPMSparseMatrixBase::TRANSPOSED);  //gradn_q = gradQ_iq * n_i
        }

        //integrate porosity field over element volumes
        if (num_threads > 1){
            parallelMultiply(M, driver->fluid_body->n, driver->fluid_body->n_e, MPMSparseMatrixBase::NORMAL, true);
        } else {
            driver->fluid_body->n_e = driver->fluid_grid->M * driver->fluid_body->n;
        }
        for (int e=0; e<driver->fluid_grid->element_count; e++) {
            //average porosity over element volume
            driver->fluid_body->n_e(e) /= driver->fluid_grid->getElementVolume(e);
        }

        //integrate gradient over element volume
        if (USE_LOCAL_POROSITY_CORRECTION){
            driver->fluid_body->n_e_x.setZero();
            for (int e=0; e<element_count; e++){
                for (int q=0; q<qpe; q++){
                    //average element porosity gradient
                    driver->fluid_body->n_e_x(e) += w_q(e*qpe+q)*gradn_q[e*qpe+q]/v_e(e);
                }
            }
        }

    } else {
        //no porosity calculation performed
        //n_q.setConstant(1.0);
        //gradn_q.setZero();
    }

    //ask material to update solid velocity field
    if (driver->fluid_material->updateSolidPhaseVelocity(job, driver) == 1){
        //successful update of solid velocity
        if (num_threads > 1){
            parallelMultiply(Q, driver->fluid_body->v_s, v_sq, MPMSparseMatrixBase::TRANSPOSED, true);
            if (USE_CORRECT_INTERPHASE_ENERGY_FLUX) {
                parallelMultiply(gradQ, driver->fluid_body->v_s, div_v_sq, MPMSparseMatrixBase::TRANSPOSED, true);
            }
        } else {
            v_sq = Q.operate(driver->fluid_body->v_s, MPMSparseMatrixBase::TRANSPOSED);       //v_sq = Q_iq * v_si
            if (USE_CORRECT_INTERPHASE_ENERGY_FLUX) {
                div_v_sq = gradQ.dot(driver->fluid_body->v_s, MPMSparseMatrixBase::TRANSPOSED);
            }
        }
    } else {
        //no porosity calculation performed
        //v_sq.setZero();
    }

    return;
}

void FVMGridBase::constructVelocityField(Job* job, FiniteVolumeDriver* driver){
    /*
    //reconstruct velocity field
    driver->fluid_body->L.setZero();
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
                driver->fluid_body->L(e, dir, pos) = sol(pos);
            }
        }

    }
     */

    if (num_threads > 1){
        //determine number of threads for element gradient calculation
        int thread_count = 0;
        if (element_count >= num_threads){
            thread_count = num_threads;
        } else {
            thread_count = element_count;
        }

        //boolean of completion status
        //volatile bool firstTaskComplete[thread_count] = {false};
        volatile bool *firstTaskComplete = new volatile bool[thread_count];
        for (int t=0; t<thread_count; t++){
            firstTaskComplete[t] = false;
        }

        //choose interval size
        int k_max = element_count - 1;
        int k_interval = (element_count/thread_count) + 1;
        int k_begin, k_end;

        for (int t=0; t<thread_count; t++) {
            //set interval
            k_begin = t * k_interval;
            k_end = k_begin + k_interval - 1;
            if (k_end > k_max){
                k_end = k_max;
            }

            //send job to threadpool
            jobThreadPool->doJob(std::bind(constructSubsetOfVelocityField,
                                           job,
                                           driver,
                                           this,
                                           k_begin,
                                           k_end,
                                           std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        delete[] firstTaskComplete;
    } else {
        //call construction function
        cSOVF(job, driver, 0, element_count - 1);
    }

    return;
}

void FVMGridBase::cSOVF(Job* job, FiniteVolumeDriver* driver, int k_begin, int k_end){
    //reconstruct velocity field
    KinematicVector u = KinematicVector(job->JOB_TYPE);
    KinematicVector p = KinematicVector(job->JOB_TYPE);
    double rho, rho_0, tmp_dif;

    Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
    KinematicVector x_0, x;
    KinematicVector u_0;

    for (int e=k_begin; e<=k_end; e++){
        driver->fluid_body->L[e].setZero();

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
                driver->fluid_body->L(e, dir, pos) = sol(pos);
            }
        }

    }
    return;
}

void FVMGridBase::constructSubsetOfVelocityField(Job* job,
                                           FiniteVolumeDriver* driver,
                                           FVMGridBase* grid,
                                           int k_begin, int k_end,
                                           volatile bool &done){
    grid->cSOVF(job, driver, k_begin, k_end);
    done = true;
    return;
}

/*----------------------------------------------------------------------------*/
//functions to compute element-wise mass flux
Eigen::VectorXd FVMGridBase::calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver){
    if (num_threads > 1){

        //initialize lhs vector
        Eigen::VectorXd lhs = Eigen::VectorXd(element_count);

        //determine number of threads for face flux evaluation
        int thread_count = 0;
        if (face_count >= num_threads){
            thread_count = num_threads;
        } else {
            thread_count = face_count;
        }

        //intermediate storage vector
        //check that memUnitID is valid
        if (memoryUnits[0].s.size() != thread_count){
            memoryUnits[0].s.resize(thread_count);
        }

        //boolean of completion status
        //volatile bool firstTaskComplete[thread_count] = {false};
        volatile bool *firstTaskComplete = new volatile bool[thread_count];
        for (int t=0; t<thread_count; t++){
            firstTaskComplete[t] = false;
        }

        //choose interval size
        int k_max = face_count - 1;
        int k_interval = (face_count/thread_count) + 1;
        int k_begin, k_end;

        for (int t=0; t<thread_count; t++) {
            //initialize output vector
            if (memoryUnits[0].s[t].size() != element_count){
                memoryUnits[0].s[t] = Eigen::VectorXd(element_count);
            }

            //set interval
            k_begin = t * k_interval;
            k_end = k_begin + k_interval - 1;
            if (k_end > k_max){
                k_end = k_max;
            }

            //send job to threadpool
            jobThreadPool->doJob(std::bind(calcMassFluxes,
                                           job,
                                           driver,
                                           this,
                                           std::ref(memoryUnits[0].s[t]),
                                           k_begin,
                                           k_end,
                                           std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        //determine number of threads for addition
        if (element_count > num_threads){
            thread_count = num_threads;
        } else {
            thread_count = element_count;
        }

        //boolean of completion status
        //volatile bool secondTaskComplete[thread_count] = {false};
        volatile bool *secondTaskComplete = new volatile bool[thread_count];
        for (int t=0; t<thread_count; t++){
            secondTaskComplete[t] = false;
        }

        //choose interval size
        int i_max = element_count - 1;
        int i_interval = (element_count/thread_count) + 1;
        int i_begin, i_end;

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //send job to thread pool
            jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::scalarAddwithFlag,
                                           std::ref(memoryUnits[0].s),
                                           std::ref(lhs),
                                           i_begin,
                                           i_end,
                                           true,
                                           std::ref(secondTaskComplete[t])));
        }

        //join threads
        taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!secondTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        delete[] firstTaskComplete;
        delete[] secondTaskComplete;

        return lhs;

    } else {
        return calculateElementMassFluxes(job, driver, 0, face_count - 1);
    }
}

Eigen::VectorXd FVMGridBase::calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver, int f_start, int f_end){
    //check that limits are within acceptable range
    if (f_start < 0 || f_end >= face_count){
        std::cerr << "ERROR! Start or end indices out of range! " << f_start << " <? 0, " << f_end << " >=? " << face_count << std::endl;
    }

    //flux rates
    Eigen::VectorXd result = Eigen::VectorXd(element_count);
    result.setZero();

    int e_plus, e_minus;    //elements
    double flux, area;      //face values

    //height index
    int H_index = GRID_DIM - 1;

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
    double epsilon_0, epsilon_bar;
    double M, theta_minus, theta_plus;
    double Pt, Tt, t0;

    //loop over faces and use quadrature to reconstruct flux integral
    u_plus = KinematicVector(job->JOB_TYPE);
    u_minus = KinematicVector(job->JOB_TYPE);
    u_bar = KinematicVector(job->JOB_TYPE);

    for (int f = f_start; f <= f_end; f++) {
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
                    n_minus = driver->fluid_body->n_e(e_minus);

                    //scale A properties
                    rho_minus /= n_minus;

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x) / rho_plus;
                    n_plus = driver->fluid_body->n_e(e_plus);

                    //scale B properties
                    rho_plus /= n_plus;

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) /
                            (std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = 0.25 * (std::sqrt(rho_plus) + std::sqrt(rho_minus))*(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    n_bar = n_q(q_list[q]);


                    c = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                n_bar*rho_bar,
                                                                n_bar*rho_bar*u_bar,
                                                                0.0,
                                                                n_bar);

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
                    flux = w_q(q_list[q]) * 0.5 * n_bar * (rho_plus * u_plus.dot(normal) + rho_minus * u_minus.dot(normal)
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

                    //find porosity
                    n_bar = n_q(q_list[q]);

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x;
                    u_minus = p_minus / rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    E_minus = rhoE_minus/rho_minus;
                    n_minus = driver->fluid_body->n_e(e_minus);

                    //scale A properties
                    rho_minus /= n_minus;
                    p_minus /= n_minus;
                    rhoE_minus /= n_minus;

                    P_minus = driver->fluid_material->getPressure(job, driver,
                                                                  n_bar*rho_minus,
                                                                  n_bar*p_minus,
                                                                  n_bar*rhoE_minus,
                                                                  n_bar); //n_q(q_list[q]));
                    H_minus = E_minus + P_minus/rho_minus;

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    E_plus = rhoE_plus/rho_plus;
                    n_plus = driver->fluid_body->n_e(e_plus);

                    //scale B properties
                    rho_plus /= n_plus;
                    p_plus /= n_plus;
                    rhoE_plus /= n_plus;

                    P_plus = driver->fluid_material->getPressure(job, driver,
                                                                 n_bar*rho_plus,
                                                                 n_bar*p_plus,
                                                                 n_bar*rhoE_plus,
                                                                 n_bar); //n_q(q_list[q]));
                    H_plus = E_plus + P_plus/rho_plus;

                    //approximate Roe advective rate
                    w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus))/2.0;
                    rho_bar = w_1_bar*w_1_bar;
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) / (2.0 * w_1_bar);
                    u_bar_squared = u_bar.dot(u_bar);
                    H_bar = (std::sqrt(rho_plus) * H_plus + std::sqrt(rho_minus) * H_minus) / (2.0 * w_1_bar);

                    c = driver->fluid_material->getSpeedOfSoundFromEnthalpy(job,
                                                                            driver,
                                                                            n_bar*rho_bar,
                                                                            n_bar*rho_bar*u_bar,
                                                                            n_bar*rho_bar*H_bar,
                                                                            n_bar);

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

                    //a_1 = 0.5*((rho_plus - rho_minus)
                    //            - a_3
                    //            - (p_plus - p_minus - u_bar*(rho_plus - rho_minus)).dot(normal) / c);
                    a_2 = 0.5 / c * ((p_plus - p_minus).dot(normal) - (rho_plus - rho_minus)*(u_bar.dot(normal) - c) - c*a_3);

                    a_1 = (rho_plus - rho_minus) - a_3 - a_2;

                    //flux in n direction
                    flux = w_q(q_list[q]) * 0.5 * n_bar * (rho_plus * u_plus.dot(normal) + rho_minus * u_minus.dot(normal)
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
                    rho_minus /= driver->fluid_body->n_e(e_minus);
                    flux = n_q(q_list[q]) * rho_minus * w_q(q_list[q]) * bc_info[f].vector.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    rho_plus /= driver->fluid_body->n_e(e_plus);
                    flux = n_q(q_list[q]) * rho_plus * w_q(q_list[q]) * bc_info[f].vector.dot(normal);
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == VELOCITY_DENSITY_INLET){
            //density and velocity given
            for (int q=0; q<q_list.size(); q++){
                //flux = rho * u.dot(n)
                flux = n_q(q_list[q]) * bc_info[f].values[0] * w_q(q_list[q]) * bc_info[f].vector.dot(normal);

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
                    rho_minus /= driver->fluid_body->n_e(e_minus);
                    flux = n_q(q_list[q]) * rho_minus * w_q(q_list[q]) * bc_info[f].vector.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    rho_plus /= driver->fluid_body->n_e(e_plus);
                    flux = n_q(q_list[q]) * rho_plus * w_q(q_list[q]) * bc_info[f].vector.dot(normal);
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

                //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                if (USE_HYDROSTATIC_CORRECTION){
                    rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                           bc_info[f].values[0]
                                                                                           + rho_bar*driver->gravity.dot(x_q[q_list[q]]),
                                                                                           bc_info[f].values[1],
                                                                                           n_q(q_list[q]));
                }

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

                //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                if (USE_HYDROSTATIC_CORRECTION){
                    rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                           bc_info[f].values[0]
                                                                                           + rho_bar*driver->gravity.dot(x_q[q_list[q]]),
                                                                                           bc_info[f].values[1],
                                                                                           n_q(q_list[q]));
                }

                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    u_minus = driver->fluid_body->p(e_minus) + driver->fluid_body->p_x[e_minus]*x;
                    u_minus /= rho_minus;
                    rho_minus /= driver->fluid_body->n_e(e_minus);
                    if (u_minus.dot(normal) < 0){
                        //flow direction is into A
                        flux = rho_bar * w_q(q_list[q]) * u_minus.dot(normal);
                    } else {
                        //flow direction is out of A
                        flux = n_q(q_list[q]) * rho_minus * w_q(q_list[q]) * u_minus.dot(normal);
                    }
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = driver->fluid_body->p(e_plus) + driver->fluid_body->p_x[e_plus]*x;
                    u_plus /= rho_plus;
                    rho_plus /= driver->fluid_body->n_e(e_plus);
                    if (u_plus.dot(normal) > 0){
                        //flow direction is into B
                        flux = rho_bar * w_q(q_list[q]) * u_plus.dot(normal);
                    } else {
                        //flow direction is out of B
                        flux = n_q(q_list[q]) * rho_plus * w_q(q_list[q]) * u_plus.dot(normal);
                    }
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == ADIABATIC_WALL ||
                    bc_info[f].tag == THERMAL_WALL ||
                    bc_info[f].tag == SYMMETRIC_WALL ||
                    bc_info[f].tag == DAMPED_WALL){
            //do nothing, no flux across this boundary
        } else if (bc_info[f].tag == SUPERSONIC_INLET){
            //density and velocity given
            for (int q=0; q<q_list.size(); q++){
                //flux = rho * u.dot(n)
                flux = n_q(q_list[q]) * bc_info[f].values[0] * w_q(q_list[q]) * bc_info[f].vector.dot(normal);

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
                    rho_minus /= driver->fluid_body->n_e(e_minus);
                    flux = n_q(q_list[q]) * rho_minus * w_q(q_list[q]) * u_minus.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = driver->fluid_body->p(e_plus) + driver->fluid_body->p_x[e_plus]*x;
                    u_plus /= rho_plus;
                    rho_plus /= driver->fluid_body->n_e(e_plus);
                    flux = n_q(q_list[q]) * rho_plus * w_q(q_list[q]) * u_plus.dot(normal);
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PERIODIC){
            //if you get here, this face is ill-defined. do nothing
        } else if (bc_info[f].tag == STAGNATION_INLET
                    || bc_info[f].tag == PULSE_STAGNATION_INLET
                    || bc_info[f].tag == RAMP_STAGNATION_INLET){

            //check pulse stagnation condition to switch BCs
            if (bc_info[f].tag == PULSE_STAGNATION_INLET && std::abs(job->t - bc_info[f].vector[2]) < job->dt){
                bc_info[f].values[0] = bc_info[f].vector[0];
                bc_info[f].values[1] = bc_info[f].vector[1];
            }

            //assign BCs
            if (bc_info[f].tag == RAMP_STAGNATION_INLET && job->t < bc_info[f].vector[2]){
                t0 = (job->t / bc_info[f].vector[2]);
                Pt = bc_info[f].vector[0] * (1.0 - t0) + bc_info[f].values[0] * (t0);
                Tt = bc_info[f].vector[1] * (1.0 - t0) + bc_info[f].values[1] * (t0);
            } else {
                Pt = bc_info[f].values[0];
                Tt = bc_info[f].values[1];
            }

            //face has prescribed pressure and temperature
            //reconstruct density and velocity (but adjust if flow direction changes)
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p(e_minus) + driver->fluid_body->p_x[e_minus]*x;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);

                    u_minus = p_minus/rho_minus;

                    if (u_minus.dot(normal) < 0){
                        //flow direction is into A
                        //estimate Mach number
                        c = driver->fluid_material->getSpeedOfSound(job,driver,rho_minus,p_minus,rhoE_minus,n_minus);
                        M = u_minus.norm()/c;

                        //estimate local pressure
                        P_minus = driver->fluid_material->getPressureFromStagnationProperties(job,
                                                                                              driver,
                                                                                              Pt,//bc_info[f].values[0],
                                                                                              M);

                        //estimate local temperature
                        theta_minus = driver->fluid_material->getTemperatureFromStagnationProperties(job,
                                                                                                     driver,
                                                                                                     Tt,//bc_info[f].values[1],
                                                                                                     M);
                        //estimate local density
                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job,
                                                                                               driver,
                                                                                               P_minus,
                                                                                               theta_minus,
                                                                                               n_q(q_list[q]));

                        //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                        if (USE_HYDROSTATIC_CORRECTION){
                            rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                                   P_minus + rho_bar*driver->gravity.dot(x_q(q_list[q])),
                                                                                                   theta_minus,
                                                                                                   n_q(q_list[q]));
                        }

                        //check mach limiter
                        if (USE_MACH_LIMITER && M>1.0){
                            u_minus *= mach_limit/M;
                            M = mach_limit;
                        }

                        flux = rho_bar * w_q(q_list[q]) * u_minus.dot(normal);
                    } else {
                        //flow direction is out of A
                        if (USE_MACH_LIMITER){
                            flux = 0; //not an outlet
                        } else {
                            flux = n_q(q_list[q]) * rho_minus / n_minus * w_q(q_list[q]) * u_minus.dot(normal);
                        }
                    }

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p(e_plus) + driver->fluid_body->p_x[e_plus]*x;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);

                    u_plus = p_plus/rho_plus;

                    if (u_plus.dot(normal) > 0){
                        //flow direction is into B
                        //estimate Mach number
                        c = driver->fluid_material->getSpeedOfSound(job,driver,rho_plus,p_plus,rhoE_plus,n_plus);
                        M = u_plus.norm()/c;

                        //estimate local pressure
                        P_plus = driver->fluid_material->getPressureFromStagnationProperties(job,
                                                                                              driver,
                                                                                              Pt,//bc_info[f].values[0],
                                                                                              M);

                        //estimate local temperature
                        theta_plus = driver->fluid_material->getTemperatureFromStagnationProperties(job,
                                                                                                     driver,
                                                                                                     Tt,//bc_info[f].values[1],
                                                                                                     M);
                        //estimate local density
                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job,
                                                                                               driver,
                                                                                               P_plus,
                                                                                               theta_plus,
                                                                                               n_q(q_list[q]));

                        //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                        if (USE_HYDROSTATIC_CORRECTION){
                            rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                                   P_plus + rho_bar*driver->gravity.dot(x_q(q_list[q])),
                                                                                                   theta_plus,
                                                                                                   n_q(q_list[q]));
                        }

                        //check mach limiter
                        if (USE_MACH_LIMITER && M>1.0){
                            u_plus *= mach_limit/M;
                            M = mach_limit;
                        }

                        flux = rho_bar * w_q(q_list[q]) * u_plus.dot(normal);
                    } else {
                        //flow direction is out of B
                        if (USE_MACH_LIMITER){
                            flux = 0; //not an outlet
                        } else {
                            flux = n_q(q_list[q]) * rho_plus / n_plus * w_q(q_list[q]) * u_plus.dot(normal);
                        }
                    }

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
KinematicVectorArray FVMGridBase::calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver) {
    if (num_threads > 1) {

        //initialize lhs vector
        KinematicVectorArray lhs = KinematicVectorArray(element_count, job->JOB_TYPE);

        //determine number of threads for face flux evaluation
        int thread_count = 0;
        if (face_count >= num_threads) {
            thread_count = num_threads;
        } else {
            thread_count = face_count;
        }

        //intermediate storage vector
        //check that memUnitID is valid
        if (memoryUnits[0].kv.size() != thread_count) {
            memoryUnits[0].kv.resize(thread_count);
        }

        //boolean of completion status
        //volatile bool firstTaskComplete[thread_count] = {false};
        volatile bool *firstTaskComplete = new volatile bool[thread_count];
        for (int t=0; t<thread_count; t++){
            firstTaskComplete[t] = false;
        }

        //choose interval size
        int k_max = face_count - 1;
        int k_interval = (face_count / thread_count) + 1;
        int k_begin, k_end;

        for (int t = 0; t < thread_count; t++) {
            //initialize output vector
            if (memoryUnits[0].kv[t].size() != element_count) {
                memoryUnits[0].kv[t] = KinematicVectorArray(element_count, job->JOB_TYPE);
            }

            //set interval
            k_begin = t * k_interval;
            k_end = k_begin + k_interval - 1;
            if (k_end > k_max) {
                k_end = k_max;
            }

            //send job to threadpool
            jobThreadPool->doJob(std::bind(calcMomentumFluxes,
                                           job,
                                           driver,
                                           this,
                                           std::ref(memoryUnits[0].kv[t]),
                                           k_begin,
                                           k_end,
                                           std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone) {
            //set flag to true
            taskDone = true;
            for (int t = 0; t < thread_count; t++) {
                if (!firstTaskComplete[t]) {
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        //determine number of threads for addition
        if (element_count > num_threads) {
            thread_count = num_threads;
        } else {
            thread_count = element_count;
        }

        //boolean of completion status
        //volatile bool secondTaskComplete[thread_count] = {false};
        volatile bool *secondTaskComplete = new volatile bool[thread_count];
        for (int t=0; t<thread_count; t++){
            secondTaskComplete[t] = false;
        }

        //choose interval size
        int i_max = element_count - 1;
        int i_interval = (element_count / thread_count) + 1;
        int i_begin, i_end;

        for (int t = 0; t < thread_count; t++) {
            //set interval
            i_begin = t * i_interval;
            i_end = i_begin + i_interval - 1;
            if (i_end > i_max) {
                i_end = i_max;
            }
            //send job to thread pool
            jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::vectorAddKwithFlag,
                                           std::ref(memoryUnits[0].kv),
                                           std::ref(lhs),
                                           i_begin,
                                           i_end,
                                           true,
                                           std::ref(secondTaskComplete[t])));
        }

        //join threads
        taskDone = false;
        //wait for task to complete
        while (!taskDone) {
            //set flag to true
            taskDone = true;
            for (int t = 0; t < thread_count; t++) {
                if (!secondTaskComplete[t]) {
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        delete[] firstTaskComplete;
        delete[] secondTaskComplete;

        return lhs;

    } else {
        return calculateElementMomentumFluxes(job, driver, 0, face_count - 1);
    }
}

KinematicVectorArray FVMGridBase::calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver, int f_start, int f_end){

    //check that limits are within acceptable range
    if (f_start < 0 || f_end >= face_count){
        std::cerr << "ERROR! Start or end indices out of range! " << f_start << " <? 0, " << f_end << " >=? " << face_count << std::endl;
    }

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

    double epsilon_0, epsilon_bar;
    double M, theta_minus, theta_plus;
    double Pt, Tt, t0;

    //traction calculatons
    KinematicTensorArray L = driver->fluid_body->L; //getVelocityGradients(job, driver);
    KinematicTensor L_tmp = KinematicTensor(job->JOB_TYPE);
    MaterialTensor tau_plus, tau_minus;
    KinematicVector dx = x;
    KinematicVector du = x;

    //assign damping coefficient
    if (USE_DECAYING_DAMPING){
        damping_coefficient = 1.0 - (1.0 - damping_at_t0) * std::exp(job->t/damping_time_cst)
                              / (damping_at_t0 + (1.0 - damping_at_t0)*std::exp(job->t/damping_time_cst));
    }

    //loop over faces and use quadrature to reconstruct flux integral
    u_plus = KinematicVector(job->JOB_TYPE);
    u_minus = KinematicVector(job->JOB_TYPE);
    u_bar = KinematicVector(job->JOB_TYPE);

    for (int f = f_start; f <= f_end; f++) {
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
                    n_minus = driver->fluid_body->n_e(e_minus);

                    //scale A properties
                    rho_minus /= n_minus;
                    p_minus /= n_minus;

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;
                    n_plus = driver->fluid_body->n_e(e_plus);

                    //scale B properties
                    rho_plus /= n_plus;
                    p_plus /= n_plus;

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) /
                            (std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = 0.25 * (std::sqrt(rho_plus) + std::sqrt(rho_minus)) * (std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    n_bar = n_q(q_list[q]);

                    c = driver->fluid_material->getSpeedOfSound(job, driver,
                                                                n_bar * rho_bar,
                                                                n_bar * rho_bar * u_bar,
                                                                0.0,
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
                    flux = w_q(q_list[q]) * 0.5 * n_bar *(p_plus*u_plus.dot(normal) + p_minus*u_minus.dot(normal)
                                                          - a_1*lambda_1*(u_bar - c*normal)
                                                          - a_2*lambda_2*(u_bar + c*normal)
                                                          - a_4*lambda_4);

                    //add tractions to momentum flux
                    if (!USE_LOCAL_GRADIENT_CORRECTION) {
                        //use standard estimate of traction
                        tau_minus = driver->fluid_material->getShearStress(job, driver,
                                                                           L[e_minus],
                                                                           n_bar*rho_minus,
                                                                           n_bar*p_minus,
                                                                           0.0,
                                                                           n_bar);
                        tau_plus = driver->fluid_material->getShearStress(job, driver,
                                                                          L[e_plus],
                                                                          n_bar*rho_plus,
                                                                          n_bar*p_plus,
                                                                          0.0,
                                                                          n_bar);

                        P_plus = c * c * (rho_plus - rho_bar) +
                                 driver->fluid_material->getPressure(job,
                                                                     driver,
                                                                     n_bar * rho_bar,
                                                                     n_bar * rho_bar * u_bar,
                                                                     0.0,
                                                                     n_bar); //n_q(q_list[q]));

                        P_minus = c * c * (rho_minus - rho_plus) + P_plus;

                        flux += w_q(q_list[q]) * 0.5 * (n_bar*(P_plus + P_minus) * normal
                                                        - KinematicVector((tau_plus + tau_minus) * normal, job->JOB_TYPE));
                    } else {
                        P_plus = c * c * (rho_plus - rho_bar) +
                                 driver->fluid_material->getPressure(job,
                                                                     driver,
                                                                     n_bar*rho_bar,
                                                                     n_bar*rho_bar*u_bar,
                                                                     0.0,
                                                                     n_bar); //n_q(q_list[q]));

                        P_minus = c * c * (rho_minus - rho_plus) + P_plus;

                        //use 'corrected' stencil for du/dx
                        L_tmp = 0.5*(L[e_minus] + L[e_plus]);

                        //find relative position of e_plus and e_minus
                        dx = x_e[e_plus] - x_e[e_minus];

                        //find change in velocity
                        du = u_plus - u_minus;

                        //'correct' L using centroid difference
                        L_tmp += (du - L_tmp*dx).tensor(dx)/(dx.norm()*dx.norm());

                        //calculate tau
                        tau_plus = driver->fluid_material->getShearStress(job,
                                                                          driver,
                                                                          L_tmp,
                                                                          n_bar*rho_bar,
                                                                          n_bar*rho_bar*u_bar,
                                                                          0.0,
                                                                          n_bar);


                        flux += w_q(q_list[q]) * (0.5*n_bar*(P_plus + P_minus) * normal
                                                        - KinematicVector(tau_plus * normal, job->JOB_TYPE));
                    }

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

                    //get quadrature porosity
                    n_bar = n_q(q_list[q]);

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x;
                    u_minus = p_minus / rho_minus;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    E_minus = rhoE_minus / rho_minus;
                    n_minus = driver->fluid_body->n_e(e_minus);

                    //scale A properties
                    rho_minus /= n_minus;
                    p_minus /= n_minus;
                    rhoE_minus /= n_minus;

                    P_minus = driver->fluid_material->getPressure(job, driver,
                                                                  n_bar*rho_minus,
                                                                  n_bar*p_minus,
                                                                  n_bar*rhoE_minus,
                                                                  n_bar); //n_q(q_list[q]));
                    H_minus = E_minus + P_minus / rho_minus;

                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    u_plus = p_plus / rho_plus;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus)  + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    E_plus = rhoE_plus / rho_plus;
                    n_plus = driver->fluid_body->n_e(e_plus);

                    //scale B properties
                    rho_plus /= n_plus;
                    p_plus /= n_plus;
                    rhoE_plus /= n_plus;

                    P_plus = driver->fluid_material->getPressure(job, driver,
                                                                 n_bar*rho_plus,
                                                                 n_bar*p_plus,
                                                                 n_bar*rhoE_plus,
                                                                 n_bar); //n_q(q_list[q]));
                    H_plus = E_plus + P_plus / rho_plus;

                    //approximate Roe advective rate
                    w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
                    rho_bar = w_1_bar * w_1_bar;
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) / (2.0 * w_1_bar);
                    u_bar_squared = u_bar.dot(u_bar);
                    H_bar = (std::sqrt(rho_plus) * H_plus + std::sqrt(rho_minus) * H_minus) / (2.0 * w_1_bar);

                    c = driver->fluid_material->getSpeedOfSoundFromEnthalpy(job,
                                                                            driver,
                                                                            n_bar*rho_bar,
                                                                            n_bar*rho_bar*u_bar,
                                                                            n_bar*rho_bar*H_bar,
                                                                            n_bar);

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

                    //a_1 = 0.5 * ((rho_plus - rho_minus)
                    //             - a_3
                    //             - (p_plus - p_minus - u_bar * (rho_plus - rho_minus)).dot(normal) / c);
                    a_2 = 0.5 / c * ((p_plus - p_minus).dot(normal) - (rho_plus - rho_minus)*(u_bar.dot(normal) - c) - c*a_3);

                    a_1 = (rho_plus - rho_minus) - a_3 - a_2;

                    a_4 = (p_plus - p_minus) - u_bar*(rho_plus - rho_minus);
                    a_4 -= a_4.dot(normal)*normal;                              //remove normal component of vector


                    //flux in n direction
                    flux = w_q(q_list[q]) * 0.5 * n_bar * (p_plus*u_plus.dot(normal) + p_minus*u_minus.dot(normal)
                                                           - a_1*lambda_1*(u_bar - c*normal)
                                                           - a_2*lambda_2*(u_bar + c*normal)
                                                           - a_3*lambda_3*(u_bar)
                                                           - a_4*lambda_4);

                    //add tractions to momentum flux
                    if(!USE_LOCAL_GRADIENT_CORRECTION) {
                        tau_minus = driver->fluid_material->getShearStress(job, driver,
                                                                           L[e_minus],
                                                                           n_bar*rho_minus,
                                                                           n_bar*p_minus,
                                                                           n_bar*rhoE_minus,
                                                                           n_bar);
                        tau_plus = driver->fluid_material->getShearStress(job, driver,
                                                                          L[e_plus],
                                                                          n_bar*rho_plus,
                                                                          n_bar*p_plus,
                                                                          n_bar*rhoE_plus,
                                                                          n_bar);
                    } else {
                        //use 'corrected' stencil for du/dx
                        L_tmp = 0.5*(L[e_minus] + L[e_plus]);

                        //find relative position of e_plus and e_minus
                        dx = x_e[e_plus] - x_e[e_minus];

                        //find change in velocity
                        du = u_plus - u_minus;

                        //'correct' L using centroid difference
                        L_tmp += (du - L_tmp*dx).tensor(dx)/(dx.norm()*dx.norm());

                        //average density and energy
                        rho_bar = (rho_plus + rho_minus)/2.0;
                        rhoE_bar = (rhoE_plus + rhoE_minus)/2.0;

                        //calculate tau
                        tau_plus = driver->fluid_material->getShearStress(job,
                                                                          driver,
                                                                          L_tmp,
                                                                          n_bar*rho_bar,
                                                                          (n_bar*rho_bar*u_bar),
                                                                          n_bar*rhoE_bar,
                                                                          n_bar);
                        tau_minus = tau_plus;
                    }

                    flux += w_q(q_list[q]) * 0.5 * (n_bar*(P_plus + P_minus)*normal - KinematicVector((tau_plus + tau_minus)*normal, job->JOB_TYPE));

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
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);
                    n_bar = n_q(q_list[q]);

                    p_minus = bc_info[f].vector * rho_minus;

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver,
                                                                       L_tmp,
                                                                       rho_minus*n_bar/n_minus,
                                                                       n_bar*p_minus/n_minus,
                                                                       n_bar*rhoE_minus/n_minus,
                                                                       n_bar);
                    P_minus = driver->fluid_material->getPressure(job, driver,
                                                                  n_bar*rho_minus/n_minus,
                                                                  n_bar*p_minus/n_minus,
                                                                  n_bar*rhoE_minus/n_minus,
                                                                  n_bar); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar*rho_minus/n_minus*bc_info[f].vector*bc_info[f].vector.dot(normal)
                                             + n_bar*P_minus*normal
                                             - KinematicVector(tau_minus*normal, job->JOB_TYPE));

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);

                    p_plus = bc_info[f].vector*rho_plus;
                    n_bar = n_q(q_list[q]);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job,
                                                                      driver,
                                                                      L_tmp,
                                                                      n_bar*rho_plus/n_plus,
                                                                      n_bar*p_plus/n_plus,
                                                                      n_bar*rhoE_plus/n_plus,
                                                                      n_bar);

                    P_plus = driver->fluid_material->getPressure(job, driver,
                                                                 n_bar*rho_plus/n_plus,
                                                                 n_bar*p_plus/n_plus,
                                                                 n_bar*rhoE_plus/n_plus,
                                                                 n_bar); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar*rho_plus/n_plus*bc_info[f].vector*bc_info[f].vector.dot(normal)
                                             + n_bar*P_plus*normal
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
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);
                    p_minus = bc_info[f].vector*bc_info[f].values[0]; //p = u * rho
                    rho_minus = bc_info[f].values[0];

                    n_bar = n_q(q_list[q]);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver,
                                                                       L_tmp,
                                                                       n_bar*rho_minus/n_minus,
                                                                       n_bar*p_minus/n_minus,
                                                                       n_bar*rhoE_minus/n_minus,
                                                                       n_bar);
                    P_minus = driver->fluid_material->getPressure(job, driver,
                                                                  n_bar*rho_minus/n_minus,
                                                                  n_bar*p_minus/n_minus,
                                                                  n_bar*rhoE_minus/n_minus,
                                                                  n_bar); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar*p_minus/n_minus*bc_info[f].vector.dot(normal)
                                             + n_bar*P_minus*normal
                                             - KinematicVector(tau_minus*normal, job->JOB_TYPE));

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);
                    p_plus = bc_info[f].vector * bc_info[f].values[0]; //p = rho * u
                    rho_plus = bc_info[f].values[0];

                    n_bar = n_q(q_list[q]);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver,
                                                                      L_tmp,
                                                                      n_bar*rho_plus/n_plus,
                                                                      n_bar*p_plus/n_plus,
                                                                      n_bar*rhoE_plus/n_plus,
                                                                      n_bar);
                    P_plus = driver->fluid_material->getPressure(job, driver, rho_plus, p_plus, rhoE_plus, n_plus); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar*p_plus/n_plus*bc_info[f].vector.dot(normal)
                                             + n_bar*P_plus*normal
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
                    n_minus = driver->fluid_body->n_e(e_minus);
                    p_minus = rho_minus * bc_info[f].vector;

                    n_bar = n_q(q_list[q]);

                    P_minus = driver->fluid_material->getPressureFromDensityAndTemperature(job,
                                                                                           driver,
                                                                                           n_bar*rho_minus/n_minus,
                                                                                           bc_info[f].values[0],
                                                                                           n_bar); //n_q(q_list[q]));

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

                    tau_minus = driver->fluid_material->getShearStress(job, driver,
                                                                       L_tmp,
                                                                       n_bar*rho_minus/n_minus,
                                                                       n_bar*p_minus/n_minus,
                                                                       n_bar*rhoE_minus/n_minus,
                                                                       n_bar);

                    flux = w_q(q_list[q]) * (n_bar*p_minus/n_minus*bc_info[f].vector.dot(normal)
                                             + n_bar*P_minus*normal
                                             - KinematicVector(tau_minus*normal, job->JOB_TYPE));

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);

                    p_plus = rho_plus * bc_info[f].vector;

                    n_bar = n_q(q_list[q]);

                    P_plus = driver->fluid_material->getPressureFromDensityAndTemperature(job,
                                                                                          driver,
                                                                                          n_bar*rho_plus/n_plus,
                                                                                          bc_info[f].values[0],
                                                                                          n_bar); //n_q(q_list[q]));

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

                    tau_plus = driver->fluid_material->getShearStress(job, driver,
                                                                      L_tmp,
                                                                      n_bar*rho_plus/n_plus,
                                                                      n_bar*p_plus/n_plus,
                                                                      n_bar*rhoE_plus/n_plus,
                                                                      n_bar); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar*p_plus/n_plus*bc_info[f].vector.dot(normal)
                                             + n_bar*P_plus*normal
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

                    //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                    if (USE_HYDROSTATIC_CORRECTION){
                        P_minus += rho_bar*driver->gravity.dot(x_q[q_list[q]]);

                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                               P_minus,
                                                                                               bc_info[f].values[1],
                                                                                               n_q(q_list[q]));
                    }

                    flux = w_q(q_list[q]) * (rho_bar*u_minus*u_minus.dot(normal)
                                             + n_q(q_list[q])*P_minus*normal);

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

                    //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                    if (USE_HYDROSTATIC_CORRECTION){
                        P_plus += rho_bar*driver->gravity.dot(x_q[q_list[q]]);

                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                               P_plus,
                                                                                               bc_info[f].values[1],
                                                                                               n_q(q_list[q]));
                    }

                    flux = w_q(q_list[q]) * (rho_bar*u_plus*u_plus.dot(normal)
                                             + n_q(q_list[q])*P_plus*normal);

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
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x(e_minus).dot(x);
                    u_minus = p_minus / rho_minus;
                    n_minus = driver->fluid_body->n_e(e_minus);

                    //tau_minus = 0
                    P_minus = bc_info[f].values[0];

                    //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                    if (USE_HYDROSTATIC_CORRECTION){
                        P_minus += rho_bar*driver->gravity.dot(x_q[q_list[q]]);

                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                               P_minus,
                                                                                               bc_info[f].values[1],
                                                                                               n_q(q_list[q]));
                    }

                    if (bc_info[f].tag == DAMPED_OUTLET && u_minus.dot(normal) > 0){
                        //only adjust P for outflow

                        P_minus = (1.0 - damping_coefficient)*P_minus
                                 + damping_coefficient*(driver->fluid_material->getPressure(job,driver,
                                                                                            n_q(q_list[q])*rho_minus/n_minus,
                                                                                            n_q(q_list[q])*p_minus/n_minus,
                                                                                            n_q(q_list[q])*rhoE_minus/n_minus,
                                                                                            n_q(q_list[q])));
                    }

                    if (u_minus.dot(normal) > 0) {
                        //flow out of A
                        flux = w_q(q_list[q]) * (n_q(q_list[q]) * rho_minus / n_minus * u_minus * u_minus.dot(normal)
                                                 + n_q(q_list[q])*P_minus * normal);
                    } else {
                        //flow direction is into A
                        //estimate Mach number
                        c = driver->fluid_material->getSpeedOfSound(job,driver,rho_minus,p_minus,rhoE_minus,n_minus);
                        M = u_minus.norm()/c;

                        //estimate local pressure
                        P_minus = driver->fluid_material->getPressureFromStagnationProperties(job,
                                                                                              driver,
                                                                                              P_minus,
                                                                                              M);

                        //estimate local temperature
                        theta_minus = driver->fluid_material->getTemperatureFromStagnationProperties(job,
                                                                                                     driver,
                                                                                                     bc_info[f].values[1],
                                                                                                     M);
                        //estimate local density
                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job,
                                                                                               driver,
                                                                                               P_minus,
                                                                                               theta_minus,
                                                                                               n_q(q_list[q]));

                        //flow into A
                        flux = w_q(q_list[q]) * (rho_bar * u_minus * u_minus.dot(normal)
                                                 + n_q(q_list[q])*P_minus * normal);
                    }

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x(e_plus).dot(x);
                    u_plus = p_plus / rho_plus;
                    n_plus = driver->fluid_body->n_e(e_plus);

                    //tau_minus = 0
                    P_plus = bc_info[f].values[0];

                    //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                    if (USE_HYDROSTATIC_CORRECTION){
                        P_plus += rho_bar*driver->gravity.dot(x_q[q_list[q]]);

                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                               P_plus,
                                                                                               bc_info[f].values[1],
                                                                                               n_q(q_list[q]));
                    }

                    if (bc_info[f].tag == DAMPED_OUTLET && u_plus.dot(normal) < 0){
                        //only adjust P for outflow
                        n_plus = driver->fluid_body->n_e(e_plus);

                        P_plus = (1.0 - damping_coefficient)*P_plus
                                 + damping_coefficient*driver->fluid_material->getPressure(job,driver,
                                                                                           n_q(q_list[q])*rho_plus/n_plus,
                                                                                           n_q(q_list[q])*p_plus/n_plus,
                                                                                           n_q(q_list[q])*rhoE_plus/n_plus,
                                                                                           n_q(q_list[q])); //n_q(q_list[q])); //ok to use porosity at quad point
                    }

                    if (u_plus.dot(normal) < 0) {
                        //flow out of B
                        flux = w_q(q_list[q]) * (n_q(q_list[q])*rho_plus/n_plus * u_plus * u_plus.dot(normal)
                                                 + n_q(q_list[q])*P_plus * normal);
                    } else {
                        //flow into B
                        //estimate Mach number
                        c = driver->fluid_material->getSpeedOfSound(job,driver,rho_plus,p_plus,rhoE_plus,n_plus);
                        M = u_plus.norm()/c;

                        //estimate local pressure
                        P_plus = driver->fluid_material->getPressureFromStagnationProperties(job,
                                                                                              driver,
                                                                                              P_plus,
                                                                                              M);

                        //estimate local temperature
                        theta_plus = driver->fluid_material->getTemperatureFromStagnationProperties(job,
                                                                                                     driver,
                                                                                                     bc_info[f].values[1],
                                                                                                     M);
                        //estimate local density
                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job,
                                                                                               driver,
                                                                                               P_plus,
                                                                                               theta_plus,
                                                                                               n_q(q_list[q]));

                        flux = w_q(q_list[q]) * (rho_bar * u_plus * u_plus.dot(normal)
                                                 + n_q(q_list[q])*P_plus * normal);
                    }

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == ADIABATIC_WALL ||
                   bc_info[f].tag == THERMAL_WALL ||
                   bc_info[f].tag == DAMPED_WALL){
            //reconstruct pressure and traction at boundary
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);

                    p_minus = bc_info[f].vector * rho_minus;

                    n_bar = n_q(q_list[q]);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                           /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver,
                                                                       L_tmp,
                                                                       n_bar*rho_minus/n_minus,
                                                                       n_bar*p_minus/n_minus,
                                                                       n_bar*rhoE_minus/n_minus,
                                                                       n_bar);
                    P_minus = driver->fluid_material->getPressure(job, driver,
                                                                  n_bar*rho_minus/n_minus,
                                                                  n_bar*p_minus/n_minus,
                                                                  n_bar*rhoE_minus/n_minus,
                                                                  n_bar); //n_q(q_list[q]));

                    if (bc_info[f].tag == DAMPED_WALL){
                        //add damping term to pressure
                        P_minus -= bc_info[f].values[0]*L_tmp.trace();
                    }

                    flux = w_q(q_list[q]) * (n_bar*P_minus*normal
                                             - KinematicVector(tau_minus*normal, job->JOB_TYPE));

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);

                    n_plus = driver->fluid_body->n_e(e_plus);

                    p_plus = bc_info[f].vector*rho_plus;

                    n_bar = n_q(q_list[q]);

                    //estimate L
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_plus][ii]/driver->fluid_body->rho(e_plus))
                                           /(x_f[f] - x_e[e_plus]).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver,
                                                                      L_tmp,
                                                                      n_bar*rho_plus/n_plus,
                                                                      n_bar*p_plus/n_plus,
                                                                      n_bar*rhoE_plus/n_plus,
                                                                      n_bar);
                    P_plus = driver->fluid_material->getPressure(job, driver,
                                                                 n_bar*rho_plus/n_plus,
                                                                 n_bar*p_plus/n_plus,
                                                                 n_bar*rhoE_plus/n_plus,
                                                                 n_bar); // n_q(q_list[q]));

                    if (bc_info[f].tag == DAMPED_WALL){
                        //add damping term to pressure
                        P_plus -= bc_info[f].values[0]*L_tmp.trace();
                    }

                    flux = w_q(q_list[q]) * (n_bar*P_plus*normal
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
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);

                    //p_minus = p_minus - p_minus.dot(normal)*normal; //bc_info[f].vector * rho_minus;

                    P_minus = driver->fluid_material->getPressure(job, driver,
                                                                  n_q(q_list[q])*rho_minus/n_minus,
                                                                  n_q(q_list[q])*p_minus/n_minus,
                                                                  n_q(q_list[q])*rhoE_minus/n_minus,
                                                                  n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_q(q_list[q])*P_minus*normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);

                    //p_plus = p_plus - p_plus.dot(normal)*normal; //bc_info[f].vector*rho_plus;

                    P_plus = driver->fluid_material->getPressure(job, driver,
                                                                 n_q(q_list[q])*rho_plus/n_plus,
                                                                 n_q(q_list[q])*p_plus/n_plus,
                                                                 n_q(q_list[q])*rhoE_plus/n_plus,
                                                                 n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_q(q_list[q])*P_plus*normal);

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

                flux = w_q(q_list[q]) * (n_q(q_list[q]) * bc_info[f].values[0] * bc_info[f].vector * bc_info[f].vector.dot(normal)
                                         + n_q(q_list[q])*P_minus * normal);

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
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);

                    u_minus = p_minus/rho_minus;

                    P_minus = driver->fluid_material->getPressure(job, driver,
                                                                  n_q(q_list[q])*rho_minus/n_minus,
                                                                  n_q(q_list[q])*p_minus/n_minus,
                                                                  n_q(q_list[q])*rhoE_minus/n_minus,
                                                                  n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_q(q_list[q]) * p_minus / n_minus * u_minus.dot(normal)
                                             + n_q(q_list[q])*P_minus * normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);

                    u_plus = p_plus/rho_plus;

                    P_plus = driver->fluid_material->getPressure(job, driver,
                                                                 n_q(q_list[q])*rho_plus/n_plus,
                                                                 n_q(q_list[q])*p_plus/n_plus,
                                                                 n_q(q_list[q])*rhoE_plus/n_plus,
                                                                 n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_q(q_list[q])*p_plus/n_plus * u_plus.dot(normal)
                                             + n_q(q_list[q])*P_plus * normal);
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PERIODIC){
            //if you get here, this face is ill-defined. do nothing
        } else if (bc_info[f].tag == STAGNATION_INLET
                   || bc_info[f].tag == PULSE_STAGNATION_INLET
                   || bc_info[f].tag == RAMP_STAGNATION_INLET){

            //check pulse stagnation condition to switch BCs
            if (bc_info[f].tag == PULSE_STAGNATION_INLET && std::abs(job->t - bc_info[f].vector[2]) < job->dt){
                bc_info[f].values[0] = bc_info[f].vector[0];
                bc_info[f].values[1] = bc_info[f].vector[1];
            }

            //assign BCs
            if (bc_info[f].tag == RAMP_STAGNATION_INLET && job->t < bc_info[f].vector[2]){
                t0 = (job->t / bc_info[f].vector[2]);
                Pt = bc_info[f].vector[0] * (1.0 - t0) + bc_info[f].values[0] * (t0);
                Tt = bc_info[f].vector[1] * (1.0 - t0) + bc_info[f].values[1] * (t0);
            } else {
                Pt = bc_info[f].values[0];
                Tt = bc_info[f].values[1];
            }

            //face has prescribed pressure and temperature
            //reconstruct density and velocity (but adjust if flow direction changes)
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p(e_minus) + driver->fluid_body->p_x[e_minus]*x;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);

                    u_minus = p_minus/rho_minus;

                    if (u_minus.dot(normal) < 0){
                        //flow direction is into A
                        //estimate Mach number
                        c = driver->fluid_material->getSpeedOfSound(job,driver,rho_minus,p_minus,rhoE_minus,n_minus);
                        M = u_minus.norm()/c;

                        //estimate local pressure
                        P_minus = driver->fluid_material->getPressureFromStagnationProperties(job,
                                                                                              driver,
                                                                                              Pt,//bc_info[f].values[0],
                                                                                              M);

                        //estimate local temperature
                        theta_minus = driver->fluid_material->getTemperatureFromStagnationProperties(job,
                                                                                                     driver,
                                                                                                     Tt,//bc_info[f].values[1],
                                                                                                     M);
                        //estimate local density
                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job,
                                                                                               driver,
                                                                                               P_minus,
                                                                                               theta_minus,
                                                                                               n_q(q_list[q]));
                        //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                        if (USE_HYDROSTATIC_CORRECTION){
                            P_minus += rho_bar*driver->gravity.dot(x_q(q_list[q]));


                            rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                                   P_minus,
                                                                                                   theta_minus,
                                                                                                   n_q(q_list[q]));
                        }

                        //check mach limiter
                        if (USE_MACH_LIMITER && M>1.0){
                            u_minus *= mach_limit/M;
                            M = mach_limit;
                        }

                        flux = w_q(q_list[q]) * (rho_bar * u_minus * u_minus.dot(normal)
                                                 + n_q(q_list[q])*P_minus*normal);
                    } else {
                        if (USE_MACH_LIMITER){
                            //no flux, use wall BC
                            p_minus = p_minus - p_minus.dot(normal)*normal; //bc_info[f].vector * rho_minus;

                            P_minus = driver->fluid_material->getPressure(job, driver,
                                                                          n_q(q_list[q])*rho_minus/n_minus,
                                                                          n_q(q_list[q])*p_minus/n_minus,
                                                                          n_q(q_list[q])*rhoE_minus/n_minus,
                                                                          n_q(q_list[q]));

                            flux = w_q(q_list[q]) * (n_q(q_list[q])*P_minus*normal);
                        } else {
                            P_minus = Pt; //bc_info[f].values[0];

                            //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                            if (USE_HYDROSTATIC_CORRECTION) {
                                rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                                       Pt,//bc_info[f].values[0],
                                                                                                       Tt,//bc_info[f].values[1],
                                                                                                       n_q(q_list[q]));

                                P_minus += rho_bar * driver->gravity.dot(x_q(q_list[q]));

                            }

                            //flow direction is out of A
                            flux = w_q(q_list[q]) * (n_q(q_list[q]) * p_minus / n_minus * u_minus.dot(normal)
                                                     + n_q(q_list[q]) * P_minus * normal);
                        }
                    }

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p(e_plus) + driver->fluid_body->p_x[e_plus]*x;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);

                    u_plus = p_plus/rho_plus;

                    if (u_plus.dot(normal) > 0){
                        //flow direction is into B
                        //estimate Mach number
                        c = driver->fluid_material->getSpeedOfSound(job,driver,rho_plus,p_plus,rhoE_plus,n_plus);
                        M = u_plus.norm()/c;

                        //estimate local pressure
                        P_plus = driver->fluid_material->getPressureFromStagnationProperties(job,
                                                                                             driver,
                                                                                             Pt,//bc_info[f].values[0],
                                                                                             M);

                        //estimate local temperature
                        theta_plus = driver->fluid_material->getTemperatureFromStagnationProperties(job,
                                                                                                    driver,
                                                                                                    Tt,//bc_info[f].values[1],
                                                                                                    M);
                        //estimate local density
                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job,
                                                                                               driver,
                                                                                               P_plus,
                                                                                               theta_plus,
                                                                                               n_q(q_list[q]));

                        //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                        if (USE_HYDROSTATIC_CORRECTION){
                            P_plus += rho_bar*driver->gravity.dot(x_q(q_list[q]));


                            rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                                   P_plus,
                                                                                                   theta_plus,
                                                                                                   n_q(q_list[q]));
                        }

                        //check mach limiter
                        if (USE_MACH_LIMITER && M>1.0){
                            u_plus *= mach_limit/M;
                            M = mach_limit;
                        }

                        flux = w_q(q_list[q]) * (rho_bar*u_plus*u_plus.dot(normal)
                                                 + n_q(q_list[q])*P_plus*normal);
                    } else {
                        if (USE_MACH_LIMITER){
                            //no flux, use wall BC
                            p_plus = p_plus - p_plus.dot(normal)*normal; //bc_info[f].vector * rho_minus;

                            P_plus = driver->fluid_material->getPressure(job, driver,
                                                                          n_q(q_list[q])*rho_plus/n_plus,
                                                                          n_q(q_list[q])*p_plus/n_plus,
                                                                          n_q(q_list[q])*rhoE_plus/n_plus,
                                                                          n_q(q_list[q]));

                            flux = w_q(q_list[q]) * (n_q(q_list[q])*P_plus*normal);
                        } else {
                            P_plus = Pt; //bc_info[f].values[0];

                            //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                            if (USE_HYDROSTATIC_CORRECTION) {
                                rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                                       Pt,//bc_info[f].values[0],
                                                                                                       Tt,//bc_info[f].values[1],
                                                                                                       n_q(q_list[q]));

                                P_plus += rho_bar * driver->gravity.dot(x_q(q_list[q]));

                            }

                            //flow direction is out of B
                            flux = w_q(q_list[q]) * n_q(q_list[q]) * (p_plus / n_plus * u_minus.dot(normal)
                                                                      + P_plus * normal);
                        }
                    }

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
    if (num_threads > 1){

        //initialize lhs vector
        Eigen::VectorXd lhs = Eigen::VectorXd(element_count);

        //determine number of threads for face flux evaluation
        int thread_count = 0;
        if (face_count >= num_threads){
            thread_count = num_threads;
        } else {
            thread_count = face_count;
        }

        //intermediate storage vector
        //check that memUnitID is valid
        if (memoryUnits[0].s.size() != thread_count){
            memoryUnits[0].s.resize(thread_count);
        }

        //boolean of completion status
        //volatile bool firstTaskComplete[thread_count] = {false};
        volatile bool *firstTaskComplete = new volatile bool[thread_count];
        for (int t=0; t<thread_count; t++){
            firstTaskComplete[t] = false;
        }

        //choose interval size
        int k_max = face_count - 1;
        int k_interval = (face_count/thread_count) + 1;
        int k_begin, k_end;

        for (int t=0; t<thread_count; t++) {
            //initialize output vector
            if (memoryUnits[0].s[t].size() != element_count){
                memoryUnits[0].s[t] = Eigen::VectorXd(element_count);
            }

            //set interval
            k_begin = t * k_interval;
            k_end = k_begin + k_interval - 1;
            if (k_end > k_max){
                k_end = k_max;
            }
            //send job to threadpool
            jobThreadPool->doJob(std::bind(calcEnergyFluxes,
                                           job,
                                           driver,
                                           this,
                                           std::ref(memoryUnits[0].s[t]),
                                           k_begin,
                                           k_end,
                                           std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        //determine number of threads for addition
        if (element_count > num_threads){
            thread_count = num_threads;
        } else {
            thread_count = element_count;
        }

        //boolean of completion status
        //volatile bool secondTaskComplete[thread_count] = {false};
        volatile bool *secondTaskComplete = new volatile bool[thread_count];
        for (int t=0; t<thread_count; t++){
            secondTaskComplete[t] = false;
        }

        //choose interval size
        int i_max = element_count - 1;
        int i_interval = (element_count/thread_count) + 1;
        int i_begin, i_end;

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //send job to thread pool
            jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::scalarAddwithFlag,
                                           std::ref(memoryUnits[0].s),
                                           std::ref(lhs),
                                           i_begin,
                                           i_end,
                                           true,
                                           std::ref(secondTaskComplete[t])));
        }

        //join threads
        taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!secondTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        delete[] firstTaskComplete;
        delete[] secondTaskComplete;

        return lhs;

    } else {
        return calculateElementEnergyFluxes(job, driver, 0, face_count - 1);
    }
}

Eigen::VectorXd FVMGridBase::calculateElementEnergyFluxes(Job* job, FiniteVolumeDriver* driver, int f_start, int f_end){

    //check that limits are within acceptable range
    if (f_start < 0 || f_end >= face_count){
        std::cerr << "ERROR! Start or end indices out of range! " << f_start << " <? 0, " << f_end << " >=? " << face_count << std::endl;
    }

    //don't calculate for ISOTHERMAL
    if (driver->TYPE != FiniteVolumeDriver::THERMAL){
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
    double epsilon_0, epsilon_bar;
    double M, Pt, Tt, t0;

    //traction calculatons
    KinematicTensorArray L = driver->fluid_body->L; //getVelocityGradients(job, driver);
    KinematicTensor L_tmp = KinematicTensor(job->JOB_TYPE);
    MaterialTensor tau_plus, tau_minus;
    KinematicVector dx = x;
    KinematicVector du = x;

    //assign damping coefficient
    if (USE_DECAYING_DAMPING){
        damping_coefficient = 1.0 - (1.0 - damping_at_t0) * std::exp(job->t/damping_time_cst)
                                    / (damping_at_t0 + (1.0 - damping_at_t0)*std::exp(job->t/damping_time_cst));
    }

    //loop over faces and use quadrature to reconstruct flux integral
    u_plus = KinematicVector(job->JOB_TYPE);
    u_minus = KinematicVector(job->JOB_TYPE);
    u_bar = KinematicVector(job->JOB_TYPE);

    for (int f = f_start; f <= f_end; f++) {
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
            for (int q=0; q<q_list.size(); q++) {
                //different flux functions for different simulation TYPES
                if (driver->TYPE == FiniteVolumeDriver::INCOMPRESSIBLE) {
                    std::cerr << "FVMGridBase does not have INCOMPRESSIBLE flux function implemented. Exiting."
                              << std::endl;
                    exit(0);
                } else if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
                    //nothing to do here
                } else if (driver->TYPE == FiniteVolumeDriver::THERMAL) {
                    //calculate A side properties
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
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);

                    u_minus = p_minus / rho_minus;
                    E_minus = rhoE_minus / rho_minus;

                    //scale A properties
                    n_bar = n_q(q_list[q]);
                    rho_minus /= n_minus;
                    p_minus /= n_minus;
                    rhoE_minus /= n_minus;

                    P_minus = driver->fluid_material->getPressure(job,
                                                                  driver,
                                                                  n_bar*rho_minus,
                                                                  n_bar*p_minus,
                                                                  n_bar*rhoE_minus,
                                                                  n_bar); //n_q(q_list[q]));
                    H_minus = E_minus + P_minus / rho_minus;
                    tau_minus = driver->fluid_material->getShearStress(job,
                                                                       driver,
                                                                       L[e_minus],
                                                                       n_bar*rho_minus,
                                                                       n_bar*p_minus,
                                                                       n_bar*rhoE_minus,
                                                                       n_bar); //n_q(q_list[q]));

                    //calculate B properties
                    //relative position to centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);

                    u_plus = p_plus / rho_plus;
                    E_plus = rhoE_plus / rho_plus;

                    //scale B properties
                    rho_plus /= n_plus;
                    p_plus /= n_plus;
                    rhoE_plus /= n_plus;

                    P_plus = driver->fluid_material->getPressure(job,
                                                                 driver,
                                                                 n_bar*rho_plus,
                                                                 n_bar*p_plus,
                                                                 n_bar*rhoE_plus,
                                                                 n_bar); //n_q(q_list[q]));
                    H_plus = E_plus + P_plus / rho_plus;
                    tau_plus = driver->fluid_material->getShearStress(job, driver,
                                                                      L[e_plus],
                                                                      n_bar*rho_plus,
                                                                      n_bar*p_plus,
                                                                      n_bar*rhoE_plus,
                                                                      n_bar);

                    //approximate Roe advective rate
                    w_1_bar = (std::sqrt(rho_plus) + std::sqrt(rho_minus)) / 2.0;
                    rho_bar = w_1_bar * w_1_bar;
                    u_bar = (std::sqrt(rho_plus) * u_plus + std::sqrt(rho_minus) * u_minus) / (2.0 * w_1_bar);
                    u_bar_squared = u_bar.dot(u_bar);
                    H_bar = (std::sqrt(rho_plus) * H_plus + std::sqrt(rho_minus) * H_minus) / (2.0 * w_1_bar);

                    c = driver->fluid_material->getSpeedOfSoundFromEnthalpy(job,
                                                                            driver,
                                                                            n_bar*rho_bar,
                                                                            n_bar*rho_bar * u_bar,
                                                                            n_bar*rho_bar * H_bar,
                                                                            n_bar);

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

                    //a_1 = 0.5 * ((rho_plus - rho_minus)
                    //             - a_3
                    //             - (p_plus - p_minus - u_bar * (rho_plus - rho_minus)).dot(normal) / c);
                    a_2 = 0.5 / c * ((p_plus - p_minus).dot(normal) - (rho_plus - rho_minus)*(u_bar.dot(normal) - c) - c*a_3);

                    a_1 = (rho_plus - rho_minus) - a_3 - a_2;

                    a_4 = (p_plus - p_minus) - u_bar*(rho_plus - rho_minus);
                    a_4 -= a_4.dot(normal)*normal;                              //remove normal component of vector

                    //tangent component of u_bar
                    s = u_bar - u_bar.dot(normal) * normal;

                    if (USE_LOCAL_GRADIENT_CORRECTION) {
                        //use 'corrected' stencil for du/dx
                        L_tmp = 0.5 * (L[e_minus] + L[e_plus]);

                        //find relative position of e_plus and e_minus
                        dx = x_e[e_plus] - x_e[e_minus];

                        //find change in velocity
                        du = u_plus - u_minus;

                        //'correct' L using centroid difference
                        L_tmp += (du - L_tmp * dx).tensor(dx) / (dx.norm() * dx.norm());



                        //calculate tau
                        tau_plus = driver->fluid_material->getShearStress(job,
                                                                          driver,
                                                                          L_tmp,
                                                                          n_bar*(rho_plus + rho_minus)/2.0,
                                                                          (n_bar*(rho_plus + rho_minus)/2.0 * u_bar),
                                                                          n_bar*(rhoE_plus + rhoE_minus)/2.0,
                                                                          n_bar);
                        tau_minus = tau_plus;
                    }

                    //flux in n direction
                    flux = w_q(q_list[q]) * 0.5 * (n_bar*(rhoE_plus + P_plus) * u_plus.dot(normal)
                                                   + n_bar*(rhoE_minus + P_minus) * u_minus.dot(normal)
                                                   - KinematicVector((tau_plus * u_plus), job->JOB_TYPE).dot(normal)
                                                   - KinematicVector((tau_minus * u_minus), job->JOB_TYPE).dot(normal)
                                                   - n_bar * a_1 * lambda_1 * (H_bar - u_bar.dot(normal) * c)
                                                   - n_bar * a_2 * lambda_2 * (H_bar + u_bar.dot(normal) * c)
                                                   - n_bar * a_3 * lambda_3 * (0.5 * u_bar_squared)
                                                   - n_bar * a_4.dot(s) * lambda_4);

                    //calculate heat flux
                    theta_plus = driver->fluid_material->getTemperature(job, driver,
                                                                        driver->fluid_body->rho(e_plus),
                                                                        driver->fluid_body->p[e_plus],
                                                                        driver->fluid_body->rhoE(e_plus),
                                                                        n_plus);
                    theta_minus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_minus),
                                                                         driver->fluid_body->p[e_minus],
                                                                         driver->fluid_body->rhoE(e_minus),
                                                                         n_minus);
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
                    theta_x = (theta_plus - theta_minus) * x / (x.dot(x));
                    theta_bar = theta_plus +
                                theta_x.dot(x_q[q_list[q]] - x_e[e_plus]); //theta = theta_0 + dtheta/dx * dx
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver, n_bar*rho_bar, theta_bar, theta_x, n_bar);
                    // add to flux calculation
                    flux += w_q(q_list[q]) * heat_flux.dot(normal);

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                } else {
                    std::cerr << "FVMGridBase does not have flux function implemented for TYPE " << driver->TYPE
                              << "! Exiting." << std::endl;
                    exit(0);
                }
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
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);
                    n_bar = n_q(q_list[q]);

                    //assign boundary p_minus
                    p_minus = bc_info[f].vector * rho_minus;

                    //estimate L
                    for (int ii = 0; ii < GRID_DIM; ii++) {
                        for (int jj = 0; jj < GRID_DIM; jj++) {
                            L_tmp(ii, jj) = (bc_info[f].vector[ii] -
                                             driver->fluid_body->p[e_minus][ii] / driver->fluid_body->rho(e_minus))
                                            / (x_f[f] - x_e[e_minus]).dot(normal) * normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver,
                                                                       L_tmp,
                                                                       n_bar*rho_minus/n_minus,
                                                                       n_bar*p_minus/n_minus,
                                                                       n_bar*rhoE_minus/n_minus,
                                                                       n_bar);

                    P_minus = driver->fluid_material->getPressure(job, driver,
                                                                  n_bar*rho_minus/n_minus,
                                                                  n_bar*p_minus/n_minus,
                                                                  n_bar*rhoE_minus/n_minus,
                                                                  n_bar); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar * rhoE_minus/n_minus * bc_info[f].vector.dot(normal)
                                             + n_bar * P_minus * bc_info[f].vector.dot(normal)
                                             - KinematicVector(tau_minus * bc_info[f].vector, job->JOB_TYPE)).dot(normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);
                    n_bar = n_q(q_list[q]);

                    //assign boundary condition
                    p_plus = bc_info[f].vector * rho_plus;

                    //estimate L
                    for (int ii = 0; ii < GRID_DIM; ii++) {
                        for (int jj = 0; jj < GRID_DIM; jj++) {
                            L_tmp(ii, jj) = (bc_info[f].vector[ii] -
                                             driver->fluid_body->p[e_plus][ii] / driver->fluid_body->rho(e_plus))
                                            / (x_f[f] - x_e[e_plus]).dot(normal) * normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver,
                                                                      L_tmp,
                                                                      n_bar*rho_plus/n_plus,
                                                                      n_bar*p_plus/n_plus,
                                                                      n_bar*rhoE_plus/n_plus,
                                                                      n_bar); //n_q(q_list[q]));
                    P_plus = driver->fluid_material->getPressure(job, driver,
                                                                 n_bar*rho_plus/n_plus,
                                                                 n_bar*p_plus/n_plus,
                                                                 n_bar*rhoE_plus/n_plus,
                                                                 n_bar); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar * rhoE_plus/n_plus * bc_info[f].vector.dot(normal)
                                             + n_bar * P_plus * bc_info[f].vector.dot(normal)
                                             - KinematicVector(tau_plus * bc_info[f].vector, job->JOB_TYPE)).dot(
                            normal);


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
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);
                    n_bar = n_q(q_list[q]);

                    p_minus = n_bar * bc_info[f].vector * bc_info[f].values[0]; //p = u * rho
                    rho_minus = n_bar * bc_info[f].values[0];

                    //estimate L
                    for (int ii = 0; ii < GRID_DIM; ii++) {
                        for (int jj = 0; jj < GRID_DIM; jj++) {
                            L_tmp(ii, jj) = (bc_info[f].vector[ii] -
                                             driver->fluid_body->p[e_minus][ii] / driver->fluid_body->rho(e_minus))
                                            / (x_f[f] - x_e[e_minus]).dot(normal) * normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver,
                                                                       L_tmp,
                                                                       rho_minus,
                                                                       p_minus,
                                                                       n_bar*rhoE_minus/n_minus,
                                                                       n_bar);
                    P_minus = driver->fluid_material->getPressure(job, driver,
                                                                  rho_minus,
                                                                  p_minus,
                                                                  n_bar*rhoE_minus/n_minus,
                                                                  n_bar); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar * rhoE_minus/n_minus * bc_info[f].vector.dot(normal)
                                             + n_bar * P_minus * bc_info[f].vector.dot(normal)
                                             - KinematicVector(tau_minus * bc_info[f].vector, job->JOB_TYPE).dot(
                            normal));
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);
                    n_bar = n_q(q_list[q]);

                    p_plus = n_bar * bc_info[f].vector * bc_info[f].values[0]; //p = rho * u
                    rho_plus = n_bar * bc_info[f].values[0];

                    //estimate L
                    for (int ii = 0; ii < GRID_DIM; ii++) {
                        for (int jj = 0; jj < GRID_DIM; jj++) {
                            L_tmp(ii, jj) = (bc_info[f].vector[ii] -
                                             driver->fluid_body->p[e_plus][ii] / driver->fluid_body->rho(e_plus))
                                            / (x_f[f] - x_e[e_plus]).dot(normal) * normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver,
                                                                      L_tmp,
                                                                      rho_plus,
                                                                      p_plus,
                                                                      n_bar * rhoE_plus / n_plus,
                                                                      n_bar);
                    P_plus = driver->fluid_material->getPressure(job, driver,
                                                                 rho_plus,
                                                                 p_plus,
                                                                 n_bar * rhoE_plus / n_plus,
                                                                 n_bar); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar * rhoE_plus/n_plus * bc_info[f].vector.dot(normal)
                                             + n_bar * P_plus * bc_info[f].vector.dot(normal)
                                             - KinematicVector(tau_plus * bc_info[f].vector, job->JOB_TYPE).dot(
                            normal));

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
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);
                    n_bar = n_q(q_list[q]);

                    p_minus = n_bar * rho_minus / n_minus * bc_info[f].vector;

                    P_minus = driver->fluid_material->getPressureFromDensityAndTemperature(job,
                                                                                           driver,
                                                                                           n_bar * rho_minus / n_minus,
                                                                                           bc_info[f].values[0],
                                                                                           n_bar); //n_q(q_list[q]));

                    //rhoE = rho * eps + 0.5*rho*u^2
                    rhoE_minus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job,
                                                                                                     driver,
                                                                                                     P_minus,
                                                                                                     bc_info[f].values[0],
                                                                                                     n_minus); //n_q(q_list[q]));

                    rhoE_minus += 0.5 * rho_minus * bc_info[f].vector.dot(bc_info[f].vector);

                    //estimate L
                    for (int ii = 0; ii < GRID_DIM; ii++) {
                        for (int jj = 0; jj < GRID_DIM; jj++) {
                            L_tmp(ii, jj) = (bc_info[f].vector[ii] -
                                             driver->fluid_body->p[e_minus][ii] / driver->fluid_body->rho(e_minus))
                                            / (x_f[f] - x_e[e_minus]).dot(normal) * normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver,
                                                                       L_tmp,
                                                                       n_bar * rho_minus / n_minus,
                                                                       p_minus,
                                                                       n_bar * rhoE_minus / n_minus,
                                                                       n_bar);

                    flux = w_q(q_list[q]) * (n_bar * rhoE_minus / n_minus * bc_info[f].vector.dot(normal)
                                             + n_bar * P_minus * bc_info[f].vector.dot(normal)
                                             - (tau_minus * bc_info[f].vector).dot(normal));

                    //calculate heat flux
                    theta_minus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_minus),
                                                                         driver->fluid_body->p[e_minus],
                                                                         driver->fluid_body->rhoE(e_minus),
                                                                         driver->fluid_body->n_e(e_minus));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_minus) * x / (x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver,
                                                                    n_bar * rho_minus/n_minus,
                                                                    theta_bar,
                                                                    theta_x,
                                                                    n_bar);
                    // add to flux calculation
                    flux += w_q(q_list[q]) * heat_flux.dot(normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);
                    n_bar = n_q(q_list[q]);

                    p_plus = n_bar * rho_plus / n_plus * bc_info[f].vector;

                    P_plus = driver->fluid_material->getPressureFromDensityAndTemperature(job,
                                                                                          driver,
                                                                                          n_bar*rho_plus/n_plus,
                                                                                          bc_info[f].values[0],
                                                                                          n_bar); //n_q(q_list[q]));

                    //rhoE = rho * eps + 0.5*rho*u^2
                    rhoE_plus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job,
                                                                                                    driver,
                                                                                                    P_plus,
                                                                                                    bc_info[f].values[0],
                                                                                                    n_plus); //n_q(q_list[q]));
                    rhoE_plus += 0.5 * rho_plus * bc_info[f].vector.dot(bc_info[f].vector);

                    //estimate L
                    for (int ii = 0; ii < GRID_DIM; ii++) {
                        for (int jj = 0; jj < GRID_DIM; jj++) {
                            L_tmp(ii, jj) = (bc_info[f].vector[ii] -
                                             driver->fluid_body->p[e_plus][ii] / driver->fluid_body->rho(e_plus))
                                            / (x_f[f] - x_e[e_plus]).dot(normal) * normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver,
                                                                      L_tmp,
                                                                      n_bar*rho_plus/n_plus,
                                                                      p_plus,
                                                                      n_bar*rhoE_plus/n_plus,
                                                                      n_bar);

                    flux = w_q(q_list[q]) * (n_bar * rhoE_plus / n_plus * bc_info[f].vector.dot(normal)
                                             + n_bar * P_plus * bc_info[f].vector.dot(normal)
                                             - (tau_plus * bc_info[f].vector).dot(normal));

                    //calculate heat flux
                    theta_plus = driver->fluid_material->getTemperature(job, driver,
                                                                        driver->fluid_body->rho(e_plus),
                                                                        driver->fluid_body->p[e_plus],
                                                                        driver->fluid_body->rhoE(e_plus),
                                                                        driver->fluid_body->n_e(e_plus));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_plus) * x / (x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver,
                                                                    n_bar * rho_plus / n_plus,
                                                                    theta_bar,
                                                                    theta_x,
                                                                    n_bar);
                    // add to flux calculation
                    flux += w_q(q_list[q]) * heat_flux.dot(normal);

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
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);
                    n_bar = n_q(q_list[q]);

                    u_minus = p_minus / rho_minus;

                    //tau_minus = 0
                    P_minus = bc_info[f].values[0];

                    //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                    if (USE_HYDROSTATIC_CORRECTION){
                        P_minus += rho_bar*driver->gravity.dot(x_q[q_list[q]]);
                    }

                    rhoE_minus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                     P_minus,
                                                                                                     bc_info[f].values[1],
                                                                                                     n_minus); //n_q(q_list[q]));
                    rhoE_minus += 0.5 * rho_minus * u_minus.dot(u_minus);

                    flux = w_q(q_list[q]) * (n_bar * rhoE_minus / n_minus * u_minus.dot(normal)
                                             + n_bar * P_minus * u_minus.dot(normal));

                    //calculate heat flux
                    theta_minus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_minus),
                                                                         driver->fluid_body->p[e_minus],
                                                                         driver->fluid_body->rhoE(e_minus),
                                                                         driver->fluid_body->n_e(e_minus));
                    theta_bar = bc_info[f].values[1];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_minus) * x / (x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver,
                                                                    n_bar * rho_minus / n_minus,
                                                                    theta_bar,
                                                                    theta_x,
                                                                    n_bar);
                    // add to flux calculation
                    flux += w_q(q_list[q]) * heat_flux.dot(normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);
                    u_plus = p_plus / rho_plus;

                    n_bar = n_q(q_list[q]);

                    //tau_minus = 0
                    P_plus = bc_info[f].values[0];

                    //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                    if (USE_HYDROSTATIC_CORRECTION){
                        P_plus += rho_bar*driver->gravity.dot(x_q[q_list[q]]);
                    }

                    rhoE_plus = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                    P_plus,
                                                                                                    bc_info[f].values[1],
                                                                                                    n_plus); //n_q(q_list[q]));
                    rhoE_plus += 0.5 * rho_plus * u_plus.dot(u_plus);

                    flux = w_q(q_list[q]) * (n_bar * rhoE_plus / n_plus * u_plus.dot(normal)
                                             + n_bar * P_plus * u_plus.dot(normal));

                    //calculate heat flux
                    theta_plus = driver->fluid_material->getTemperature(job, driver,
                                                                        driver->fluid_body->rho(e_plus),
                                                                        driver->fluid_body->p[e_plus],
                                                                        driver->fluid_body->rhoE(e_plus),
                                                                        driver->fluid_body->n_e(e_plus));
                    theta_bar = bc_info[f].values[1];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_plus) * x / (x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver,
                                                                    n_bar * rho_plus / n_plus,
                                                                    theta_bar,
                                                                    theta_x,
                                                                    n_bar);
                    // add to flux calculation
                    flux += w_q(q_list[q]) * heat_flux.dot(normal);

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
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);
                    n_bar = n_q(q_list[q]);

                    u_minus = p_minus / rho_minus;

                    //tau_minus = 0
                    P_minus = bc_info[f].values[0];

                    //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                    if (USE_HYDROSTATIC_CORRECTION){
                        P_minus += rho_bar*driver->gravity.dot(x_q[q_list[q]]);
                    }

                    if (bc_info[f].tag == DAMPED_OUTLET && u_minus.dot(normal) > 0) {
                        //only adjust if outflow
                        P_minus = (1.0 - damping_coefficient) * P_minus
                                  + damping_coefficient * driver->fluid_material->getPressure(job, driver,
                                                                                              n_bar*rho_minus/n_minus,
                                                                                              n_bar*p_minus/n_minus,
                                                                                              n_bar*rhoE_minus/n_minus,
                                                                                              n_bar); //n_q(q_list[q]));
                    }

                    //calculate energy for inflow
                    rhoE_bar = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                   P_minus,
                                                                                                   bc_info[f].values[1],
                                                                                                   n_minus); //n_q(q_list[q]));


                    rhoE_bar += 0.5 * rho_minus * u_minus.dot(u_minus);

                    if (u_minus.dot(normal) > 0) {
                        //flow out of A
                        flux = w_q(q_list[q]) * (n_bar * rhoE_minus/n_minus * u_minus.dot(normal)
                                                 + n_bar * P_minus * u_minus.dot(normal));
                    } else {
                        //flow into A
                        flux = w_q(q_list[q]) * (n_bar * rhoE_bar/n_minus * u_minus.dot(normal)
                                                 + n_bar * P_minus * u_minus.dot(normal));
                    }

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);
                    n_bar = n_q(q_list[q]);
                    u_plus = p_plus / rho_plus;

                    //tau_minus = 0
                    P_plus = bc_info[f].values[0];

                    //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                    if (USE_HYDROSTATIC_CORRECTION){
                        P_plus += rho_bar*driver->gravity.dot(x_q[q_list[q]]);
                    }

                    if (bc_info[f].tag == DAMPED_OUTLET && u_plus.dot(normal) < 0) {
                        //only adjust P for outflow
                        P_plus = (1.0 - damping_coefficient) * P_plus
                                 + damping_coefficient * driver->fluid_material->getPressure(job, driver,
                                                                                             n_bar * rho_plus / n_plus,
                                                                                             n_bar * p_plus / n_plus,
                                                                                             n_bar * rhoE_plus / n_plus,
                                                                                             n_bar); //n_q(q_list[q]));
                    }

                    //calculate energy for inflow
                    rhoE_bar = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver,
                                                                                                   P_plus,
                                                                                                   bc_info[f].values[1],
                                                                                                   n_plus); //n_q(q_list[q]));
                    rhoE_bar += 0.5 * rho_plus * u_plus.dot(u_plus);

                    if (u_plus.dot(normal) < 0) {
                        //flow out of B
                        flux = w_q(q_list[q]) * (n_bar * rhoE_plus / n_plus * u_plus.dot(normal)
                                                 + n_bar * P_plus * u_plus.dot(normal));
                    } else {
                        //flow into A
                        flux = w_q(q_list[q]) * (n_bar * rhoE_bar / n_plus * u_plus.dot(normal)
                                                 + n_bar * P_plus * u_plus.dot(normal));
                    }

                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == ADIABATIC_WALL || bc_info[f].tag == SYMMETRIC_WALL ||
                   bc_info[f].tag == DAMPED_WALL) {
            //nothing to do right now
        } else if (bc_info[f].tag == THERMAL_WALL) {
            //face has prescribed temperature
            //calculate heat flux at quadrature points
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);

                    //calculate heat flux
                    theta_minus = driver->fluid_material->getTemperature(job, driver,
                                                                         driver->fluid_body->rho(e_minus),
                                                                         driver->fluid_body->p[e_minus],
                                                                         driver->fluid_body->rhoE(e_minus),
                                                                         driver->fluid_body->n_e(e_minus));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_minus) * x / (x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver,
                                                                    n_q(q_list[q])*rho_minus/n_minus,
                                                                    theta_bar,
                                                                    theta_x,
                                                                    n_q(q_list[q]));
                    // add to flux calculation
                    flux = w_q(q_list[q]) * heat_flux.dot(normal);

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);

                    //calculate heat flux
                    theta_plus = driver->fluid_material->getTemperature(job, driver,
                                                                        driver->fluid_body->rho(e_plus),
                                                                        driver->fluid_body->p[e_plus],
                                                                        driver->fluid_body->rhoE(e_plus),
                                                                        driver->fluid_body->n_e(e_plus));
                    theta_bar = bc_info[f].values[0];
                    //estimate thermal gradient
                    theta_x = (theta_bar - theta_plus) * x / (x.dot(x));
                    //estimate heat flux
                    heat_flux = driver->fluid_material->getHeatFlux(job, driver,
                                                                    n_q(q_list[q])*rho_plus/n_plus,
                                                                    theta_bar,
                                                                    theta_x,
                                                                    n_q(q_list[q]));
                    // add to flux calculation
                    flux = w_q(q_list[q]) * heat_flux.dot(normal);

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
                rhoE_bar += 0.5 * bc_info[f].values[0] * bc_info[f].vector.dot(bc_info[f].vector);


                flux = w_q(q_list[q]) * n_q(q_list[q]) * (rhoE_bar * bc_info[f].vector.dot(normal)
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
                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus] * x;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);
                    n_bar = n_q(q_list[q]);
                    u_minus = p_minus / rho_minus;

                    P_minus = driver->fluid_material->getPressure(job, driver,
                                                                  n_bar * rho_minus / n_minus,
                                                                  n_bar * p_minus / n_minus,
                                                                  n_bar * rhoE_minus / n_minus,
                                                                  n_bar); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar * rhoE_minus/n_minus * u_minus.dot(normal)
                                             + n_bar * P_minus * u_minus.dot(normal));
                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus] * x;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);
                    n_bar = n_q(q_list[q]);
                    u_plus = p_plus / rho_plus;

                    P_plus = driver->fluid_material->getPressure(job, driver,
                                                                 n_bar * rho_plus / n_plus,
                                                                 n_bar * p_plus / n_plus,
                                                                 n_bar * rhoE_plus / n_plus,
                                                                 n_bar); //n_q(q_list[q]));

                    flux = w_q(q_list[q]) * (n_bar * rhoE_plus / n_plus * u_plus.dot(normal)
                                             + n_bar * P_plus * u_plus.dot(normal));
                    result(e_plus) += flux;
                }
            }
        } else if (bc_info[f].tag == PERIODIC) {
            //if you get here, this face is ill-defined. do nothing
        } else if (bc_info[f].tag == STAGNATION_INLET
                   || bc_info[f].tag == PULSE_STAGNATION_INLET
                   || bc_info[f].tag == RAMP_STAGNATION_INLET){

            //check pulse stagnation condition to switch BCs
            if (bc_info[f].tag == PULSE_STAGNATION_INLET && std::abs(job->t - bc_info[f].vector[2]) < job->dt){
                bc_info[f].values[0] = bc_info[f].vector[0];
                bc_info[f].values[1] = bc_info[f].vector[1];
            }

            //assign BCs
            if (bc_info[f].tag == RAMP_STAGNATION_INLET && job->t < bc_info[f].vector[2]){
                t0 = (job->t / bc_info[f].vector[2]);
                Pt = bc_info[f].vector[0] * (1.0 - t0) + bc_info[f].values[0] * (t0);
                Tt = bc_info[f].vector[1] * (1.0 - t0) + bc_info[f].values[1] * (t0);
            } else {
                Pt = bc_info[f].values[0];
                Tt = bc_info[f].values[1];
            }

            //face has prescribed pressure and temperature
            //reconstruct density and velocity (but adjust if flow direction changes)
            for (int q = 0; q < q_list.size(); q++) {
                if (e_minus > -1) {
                    //relative position from centroid of A
                    x = x_q[q_list[q]] - x_e[e_minus];

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p(e_minus) + driver->fluid_body->p_x[e_minus]*x;
                    rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                    n_minus = driver->fluid_body->n_e(e_minus);

                    u_minus = p_minus/rho_minus;

                    if (u_minus.dot(normal) < 0){
                        //flow direction is into A
                        //estimate Mach number
                        c = driver->fluid_material->getSpeedOfSound(job,driver,rho_minus,p_minus,rhoE_minus,n_minus);
                        M = u_minus.norm()/c;

                        //estimate local pressure
                        P_minus = driver->fluid_material->getPressureFromStagnationProperties(job,
                                                                                              driver,
                                                                                              Pt,//bc_info[f].values[0],
                                                                                              M);

                        //estimate local temperature
                        theta_minus = driver->fluid_material->getTemperatureFromStagnationProperties(job,
                                                                                                     driver,
                                                                                                     Tt,//bc_info[f].values[1],
                                                                                                     M);
                        //estimate local density
                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job,
                                                                                               driver,
                                                                                               P_minus,
                                                                                               theta_minus,
                                                                                               n_q(q_list[q]));

                        //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                        if (USE_HYDROSTATIC_CORRECTION){
                            P_minus += rho_bar*driver->gravity.dot(x_q(q_list[q]));


                            rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                                   P_minus,
                                                                                                   theta_minus,
                                                                                                   n_q(q_list[q]));
                        }

                        //estimate local energy
                        rhoE_bar = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job,
                                                                                                       driver,
                                                                                                       P_minus,
                                                                                                       theta_minus,
                                                                                                       n_q(q_list[q]));
                        rhoE_bar += 0.5*rho_bar*u_minus.dot(u_minus);

                        //check mach limiter
                        if (USE_MACH_LIMITER && M>1.0){
                            u_minus *= mach_limit/M;
                            M = mach_limit;
                        }

                        flux = w_q(q_list[q]) * (rhoE_bar * u_minus.dot(normal)
                                                 + n_q(q_list[q])*P_minus*u_minus.dot(normal));
                    } else {
                        //flow direction is out of A
                        if (USE_MACH_LIMITER){
                            flux = 0;
                        } else {
                            flux = w_q(q_list[q]) * (n_q(q_list[q]) * rhoE_minus / n_minus * u_minus.dot(normal)
                                                     + n_q(q_list[q]) * Pt * u_minus.dot(
                                    normal)); //bc_info[f].values[0] * u_minus.dot(normal));
                        }
                    }

                    result(e_minus) -= flux;
                }

                if (e_plus > -1) {
                    //relative position from centroid of B
                    x = x_q[q_list[q]] - x_e[e_plus];

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p(e_plus) + driver->fluid_body->p_x[e_plus]*x;
                    rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x);
                    n_plus = driver->fluid_body->n_e(e_plus);

                    u_plus = p_plus/rho_plus;

                    if (u_plus.dot(normal) > 0){
                        //flow direction is into B
                        //estimate Mach number
                        c = driver->fluid_material->getSpeedOfSound(job,driver,rho_plus,p_plus,rhoE_plus,n_plus);
                        M = u_plus.norm()/c;

                        //estimate local pressure
                        P_plus = driver->fluid_material->getPressureFromStagnationProperties(job,
                                                                                             driver,
                                                                                             Pt,//bc_info[f].values[0],
                                                                                             M);

                        //estimate local temperature
                        theta_plus = driver->fluid_material->getTemperatureFromStagnationProperties(job,
                                                                                                    driver,
                                                                                                    Tt,//bc_info[f].values[1],
                                                                                                    M);
                        //estimate local density
                        rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job,
                                                                                               driver,
                                                                                               P_plus,
                                                                                               theta_plus,
                                                                                               n_q(q_list[q]));

                        //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                        if (USE_HYDROSTATIC_CORRECTION){
                            P_plus += rho_bar*driver->gravity.dot(x_q(q_list[q]));


                            rho_bar = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                                   P_plus,
                                                                                                   theta_plus,
                                                                                                   n_q(q_list[q]));
                        }

                        //estimate local energy
                        rhoE_bar = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job,
                                                                                                       driver,
                                                                                                       P_plus,
                                                                                                       theta_plus,
                                                                                                       n_q(q_list[q]));
                        rhoE_bar += 0.5*rho_bar*u_plus.dot(u_plus);

                        //check mach limiter
                        if (USE_MACH_LIMITER && M>1.0){
                            u_plus *= mach_limit/M;
                            M = mach_limit;
                        }

                        flux = w_q(q_list[q]) * (rhoE_bar*u_plus.dot(normal)
                                                 + n_q(q_list[q])*P_plus*u_plus.dot(normal));
                    } else {
                        //flow direction is out of B
                        if (USE_MACH_LIMITER){
                            flux = 0;
                        } else {
                            flux = w_q(q_list[q]) * n_q(q_list[q]) * (rhoE_plus / n_plus * u_plus.dot(normal)
                                                                      + Pt * u_plus.dot(
                                    normal)); //bc_info[f].values[0] * u_plus.dot(normal));
                        }
                    }

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

/*----------------------------------------------------------------------------*/
//function to calculate interphase force
KinematicVectorArray FVMGridBase::calculateInterphaseForces(Job* job, FiniteVolumeDriver* driver){

    //result defines force on solid by fluid
    KinematicVectorArray result = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray tmp_result = result;

    //declare temporary containers for quadrature point fields
    KinematicVectorArray tmpVecArray = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE);
    Eigen::VectorXd tmpArray = Eigen::VectorXd(int_quad_count + ext_quad_count);

    tmpVecArray.setZero();
    tmpArray.setZero();

    double rho, rhoE, n;
    KinematicVector p, normal;

    //loop over elements
    calculateElementIntegrandsForInterphaseForce(job, driver, tmpVecArray, tmpArray, 0, element_count-1);
    calculateFaceIntegrandsForInterphaseForce(job, driver, tmpVecArray, tmpArray, 0, face_count-1);

    //integrate expression (pg 117, nb 6)
    //result = Q*tmpVecArray;
    if (num_threads > 1){
        //add components to result vector
        parallelMultiply(Q, tmpVecArray, result, Q.NORMAL, false);
        parallelMultiply(gradQ, tmpArray, result, gradQ.NORMAL, false);

        //scale result vector
        for (int i = 0; i < result.size(); i++) {
            result[i] *= (1.0 - driver->fluid_body->n(i));
        }

        //new drag calculation
        if (!USE_OLD_DRAG_METHOD){
            //calculate volume weighted velocity and mass density
            KinematicVectorArray u_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
            KinematicVectorArray u = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
            Eigen::VectorXd m_i = Eigen::VectorXd(job->grid->node_count);
            for (int e=0; e<driver->fluid_grid->element_count; e++){
                u[e] = driver->fluid_body->p(e) / driver->fluid_body->rho(e);
            }
            parallelMultiply(M, u, u_i, M.TRANSPOSED, true);
            for (int i=0; i<job->grid->node_count; i++){
                u_i[i] /= job->grid->nodeVolume(job, i);
            }
            parallelMultiply(M, driver->fluid_body->rho, m_i, M.TRANSPOSED, true);

            //calculate drag contribution
            double rho_f;
            for (int i=0; i<job->grid->node_count; i++){
                rho_f = m_i(i)/job->grid->nodeVolume(job,i);
                result[i] -= job->grid->nodeVolume(job,i)*
                             driver->fluid_material->getInterphaseDrag(job,
                                                                       driver,
                                                                       rho_f,
                                                                       u_i[i],
                                                                       driver->fluid_body->v_s[i],
                                                                       driver->fluid_body->n(i));
            }
        }

    } else {
        tmp_result = Q * tmpVecArray;
        for (int i = 0; i < result.size(); i++) {
            result[i] += (1.0 - driver->fluid_body->n(i)) * tmp_result(i); //ensures nodes with zero porosity have zero force
        }
        tmp_result = gradQ * tmpArray;
        for (int i = 0; i < result.size(); i++) {
            result[i] += (1.0 - driver->fluid_body->n(i)) * tmp_result(i);
        }

        //new drag calculation
        if (!USE_OLD_DRAG_METHOD){
            //calculate volume weighted velocity and mass density
            KinematicVectorArray u_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
            KinematicVectorArray u = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
            Eigen::VectorXd m_i = Eigen::VectorXd(job->grid->node_count);
            for (int e=0; e<driver->fluid_grid->element_count; e++){
                u[e] = driver->fluid_body->p(e) / driver->fluid_body->rho(e);
            }
            u_i = M.operate(u,M.TRANSPOSED);
            for (int i=0; i<job->grid->node_count; i++){
                u_i[i] /= job->grid->nodeVolume(job, i);
            }
            m_i = M.operate(driver->fluid_body->rho, M.TRANSPOSED);

            //calculate drag contribution
            double rho_f;
            for (int i=0; i<job->grid->node_count; i++){
                rho_f = m_i(i)/job->grid->nodeVolume(job,i);
                result[i] -= job->grid->nodeVolume(job,i)*
                             driver->fluid_material->getInterphaseDrag(job,
                                                                       driver,
                                                                       rho_f,
                                                                       u_i[i],
                                                                       driver->fluid_body->v_s[i],
                                                                       driver->fluid_body->n(i));
            }
        }
    }

    //divide result by associated nodal volumes
    for (int i=0; i<job->grid->node_count; i++){
        result(i) /= job->grid->nodeVolume(job, i);
    }

    return result;
}

KinematicVectorArray FVMGridBase::calculateBuoyantForces(Job* job, FiniteVolumeDriver* driver){

    //result defines force on solid by fluid
    KinematicVectorArray result = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray tmp_result = result;

    //declare temporary containers for quadrature point fields
    KinematicVectorArray tmpVecArray = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE);
    Eigen::VectorXd tmpArray = Eigen::VectorXd(int_quad_count + ext_quad_count);

    tmpVecArray.setZero();
    tmpArray.setZero();

    double rho, rhoE, n;
    KinematicVector p, normal;

    //loop over elements
    calculateElementIntegrandsForBuoyantForce(job, driver, tmpArray, 0, element_count-1);
    calculateFaceIntegrandsForBuoyantForce(job, driver, tmpVecArray, 0, face_count-1);

    //integrate expression (pg 117, nb 6)
    //result = Q*tmpVecArray;
    if (num_threads > 1){
        //add components to result vector
        parallelMultiply(Q, tmpVecArray, result, Q.NORMAL, false);
        parallelMultiply(gradQ, tmpArray, result, gradQ.NORMAL, false);

        //scale result vector
        for (int i = 0; i < result.size(); i++) {
            result[i] *= (1.0 - driver->fluid_body->n(i));
        }
    } else {
        tmp_result = Q * tmpVecArray;
        for (int i = 0; i < result.size(); i++) {
            result[i] += (1.0 - driver->fluid_body->n(i)) * tmp_result(i); //ensures nodes with zero porosity have zero force
        }
        tmp_result = gradQ * tmpArray;
        for (int i = 0; i < result.size(); i++) {
            result[i] += (1.0 - driver->fluid_body->n(i)) * tmp_result(i);
        }
    }

    //divide result by associated nodal volumes
    for (int i=0; i<job->grid->node_count; i++){
        result(i) /= job->grid->nodeVolume(job, i);
    }

    return result;
}

KinematicVectorArray FVMGridBase::calculateDragForces(Job* job, FiniteVolumeDriver* driver){

    //result defines force on solid by fluid
    KinematicVectorArray result = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray tmp_result = result;

    //declare temporary containers for quadrature point fields
    KinematicVectorArray tmpVecArray = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE);
    Eigen::VectorXd tmpArray = Eigen::VectorXd(int_quad_count + ext_quad_count);

    tmpVecArray.setZero();
    tmpArray.setZero();

    double rho, rhoE, n;
    KinematicVector p, normal;

    if (USE_OLD_DRAG_METHOD) {
        //loop over elements
        calculateElementIntegrandsForDragForce(job, driver, tmpVecArray, 0, element_count - 1);

        //integrate expression (pg 117, nb 6)
        //result = Q*tmpVecArray;
        if (num_threads > 1) {
            //add components to result vector
            parallelMultiply(Q, tmpVecArray, result, Q.NORMAL, false);

            //scale result vector
            for (int i = 0; i < result.size(); i++) {
                result[i] *= (1.0 - driver->fluid_body->n(i));
            }
        } else {
            tmp_result = Q * tmpVecArray;
            for (int i = 0; i < result.size(); i++) {
                result[i] += (1.0 - driver->fluid_body->n(i)) *
                             tmp_result(i); //ensures nodes with zero porosity have zero force
            }
        }

        //divide result by associated nodal volumes
        for (int i=0; i<job->grid->node_count; i++){
            result(i) /= job->grid->nodeVolume(job, i);
        }

    } else {
        //new drag calculation
        //calculate volume weighted velocity and mass density
        KinematicVectorArray u_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
        KinematicVectorArray u = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
        Eigen::VectorXd m_i = Eigen::VectorXd(job->grid->node_count);
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            u[e] = driver->fluid_body->p(e) / driver->fluid_body->rho(e);
        }
        if (num_threads > 1) {
            parallelMultiply(M, u, u_i, M.TRANSPOSED, true);
        } else {
            u_i = M.operate(u, M.TRANSPOSED);
        }
        for (int i=0; i<job->grid->node_count; i++){
            u_i[i] /= job->grid->nodeVolume(job, i);
        }
        if (num_threads > 1) {
            parallelMultiply(M, driver->fluid_body->rho, m_i, M.TRANSPOSED, true);
        } else {
            m_i = M.operate(driver->fluid_body->rho, M.TRANSPOSED);
        }


        //calculate drag contribution
        double rho_f;
        for (int i=0; i<job->grid->node_count; i++){
            rho_f = m_i(i)/job->grid->nodeVolume(job,i);
            result[i] = -1.0 * driver->fluid_material->getInterphaseDrag(job,
                                                                         driver,
                                                                         rho_f,
                                                                         u_i[i],
                                                                         driver->fluid_body->v_s[i],
                                                                         driver->fluid_body->n(i));
        }
    }

    return result;
}

/*----------------------------------------------------------------------------*/

void FVMGridBase::calculateElementIntegrandsForInterphaseForce(Job *job,
                                                     FiniteVolumeDriver *driver,
                                                     KinematicVectorArray& kv,
                                                     Eigen::VectorXd& v,
                                                     int e_start, int e_end){

    //check that limits are within acceptable range
    if (e_start < 0 || e_end >= element_count){
        std::cerr << "ERROR! Start or end indices out of range! " << e_start << " <? 0, " << e_end << " >=? " << element_count << std::endl;
    }

    double rho, rhoE, n;
    KinematicVector p, normal, f_d;
    double P;
    double epsilon_0, epsilon_bar;

    //loop over elements
    for (int e=e_start; e<=e_end; e++){
        for (int q=0; q<qpe; q++){
            //approximate local density, momentum, and energy
            rho = driver->fluid_body->rho(e) + driver->fluid_body->rho_x[e].dot(x_q[e*qpe + q] - x_e[e]);
            rhoE = driver->fluid_body->rhoE(e) + driver->fluid_body->rhoE_x[e].dot(x_q[e*qpe + q] - x_e[e]);
            p = driver->fluid_body->p[e] + driver->fluid_body->p_x[e]*(x_q[e*qpe + q] - x_e[e]);
            n = driver->fluid_body->n_e(e);

            //fill in drag force
            if (USE_OLD_DRAG_METHOD) {
                f_d = driver->fluid_material->getInterphaseDrag(job,
                                                                driver,
                                                                rho,
                                                                (p / rho),
                                                                v_sq[e * qpe + q],
                                                                n); //n_q(e*qpe + q));
            }

            //fill in pressure
            P = driver->fluid_material->getPressure(job,
                                                    driver,
                                                    rho,
                                                    p,
                                                    rhoE,
                                                    n); //n_q(e*qpe + q));

            //fill in tmpVecArray and tmpArray
            if ((1.0 - n) > 1e-10 && USE_OLD_DRAG_METHOD) {
                //ensures nodes with zero porosity have zero force
                kv[e * qpe + q] = -f_d * getQuadratureWeight(e * qpe + q) / (1.0 - n);
            } else {
                //zero
                kv[e * qpe + q].setZero();
            }
            v(e*qpe + q) = P * getQuadratureWeight(e*qpe + q);
        }
    }
    return;
}

void FVMGridBase::calculateFaceIntegrandsForInterphaseForce(Job *job,
                                                  FiniteVolumeDriver *driver,
                                                  KinematicVectorArray& kv,
                                                  Eigen::VectorXd& v,
                                                  int f_start, int f_end){

    //check that limits are within acceptable range
    if (f_start < 0 || f_end >= face_count){
        std::cerr << "ERROR! Start or end indices out of range! " << f_start << " <? 0, " << f_end << " >=? " << face_count << std::endl;
    }

    double rho, rhoE, n;
    KinematicVector p, normal;
    double P;
    double epsilon_0, epsilon_bar;

    int tmp_e = -1;
    int tmp_q = -1;
    //loop over faces
    for (int f=f_start; f<=f_end; f++){
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
                n = driver->fluid_body->n_e(tmp_e);

                //fill in pressure
                P = driver->fluid_material->getPressure(job,
                                                        driver,
                                                        rho,
                                                        p,
                                                        rhoE,
                                                        n); //n_q(tmp_q));

                //fill in tmpArray
                kv[tmp_q] = -P * normal * getQuadratureWeight(tmp_q);
            }
        }
    }
    return;
}

void FVMGridBase::calculateElementIntegrandsForBuoyantForce(Job *job,
                                                               FiniteVolumeDriver *driver,
                                                               Eigen::VectorXd& v,
                                                               int e_start, int e_end){

    //check that limits are within acceptable range
    if (e_start < 0 || e_end >= element_count){
        std::cerr << "ERROR! Start or end indices out of range! " << e_start << " <? 0, " << e_end << " >=? " << element_count << std::endl;
    }

    double rho, rhoE, n;
    KinematicVector p, normal, f_d;
    double P;
    double epsilon_0, epsilon_bar;

    //loop over elements
    for (int e=e_start; e<=e_end; e++){
        for (int q=0; q<qpe; q++){
            //approximate local density, momentum, and energy
            rho = driver->fluid_body->rho(e) + driver->fluid_body->rho_x[e].dot(x_q[e*qpe + q] - x_e[e]);
            rhoE = driver->fluid_body->rhoE(e) + driver->fluid_body->rhoE_x[e].dot(x_q[e*qpe + q] - x_e[e]);
            p = driver->fluid_body->p[e] + driver->fluid_body->p_x[e]*(x_q[e*qpe + q] - x_e[e]);
            n = driver->fluid_body->n_e(e);

            //fill in pressure
            P = driver->fluid_material->getPressure(job,
                                                    driver,
                                                    rho,
                                                    p,
                                                    rhoE,
                                                    n); //n_q(e*qpe + q));

            v(e*qpe + q) = P * getQuadratureWeight(e*qpe + q);
        }
    }
    return;
}

void FVMGridBase::calculateFaceIntegrandsForBuoyantForce(Job *job,
                                                            FiniteVolumeDriver *driver,
                                                            KinematicVectorArray& kv,
                                                            int f_start, int f_end,
                                                            bool BOUNDARY_ONLY){

    //check that limits are within acceptable range
    if (f_start < 0 || f_end >= face_count){
        std::cerr << "ERROR! Start or end indices out of range! " << f_start << " <? 0, " << f_end << " >=? " << face_count << std::endl;
    }

    double rho, rhoE, n, n_bar;
    KinematicVector p, normal;
    double P;
    double P_plus, P_minus;

    double rho_plus, rhoE_plus, rho_minus, rhoE_minus;
    double n_plus, n_minus;
    KinematicVector p_plus, p_minus, x;
    double epsilon_0, epsilon_bar;

    KinematicTensor L_tmp = KinematicTensor(job->JOB_TYPE);

    //assign damping coefficient
    if (USE_DECAYING_DAMPING){
        damping_coefficient = 1.0 - (1.0 - damping_at_t0) * std::exp(job->t/damping_time_cst)
                                    / (damping_at_t0 + (1.0 - damping_at_t0)*std::exp(job->t/damping_time_cst));
    }

    int tmp_e = -1;
    int tmp_q = -1;
    int e_plus = -1;
    int e_minus = -1;
    //loop over faces
    for (int f=f_start; f<=f_end; f++){
        for (int q=0; q<qpf; q++){
            //only add contributions from domain boundaries (unless BOUNDARY_ONLY flag is false)
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
                n = driver->fluid_body->n_e(tmp_e);
                n_bar = n_q(tmp_q);

                //initialize P
                P = 0;

                if (bc_info[f].tag == VELOCITY_INLET) {
                    //fill in pressure
                    P = driver->fluid_material->getPressure(job,
                                                            driver,
                                                            rho,
                                                            p,
                                                            rhoE,
                                                            n); //n_q(tmp_q));

                } else if (bc_info[f].tag == VELOCITY_DENSITY_INLET) {
                    //density and velocity given
                    rho = n_bar*bc_info[f].values[0];
                    p = n_bar*bc_info[f].values[0] * bc_info[f].vector;

                    //fill in pressure
                    P = driver->fluid_material->getPressure(job,
                                                            driver,
                                                            rho,
                                                            p,
                                                            n_bar*rhoE/n,
                                                            n_bar); //n_q(tmp_q));

                } else if (bc_info[f].tag == VELOCITY_TEMP_INLET) {
                    //velocity and temperature given
                    p = rho*bc_info[f].vector;

                    //fill in pressure
                    P = driver->fluid_material->getPressureFromDensityAndTemperature(job, driver, rho, bc_info[f].values[0], n);

                } else if (bc_info[f].tag == PRESSURE_INLET) {
                    //face has prescribed pressure (and temperature)
                    P = bc_info[f].values[0];

                    //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                    if (USE_HYDROSTATIC_CORRECTION){
                        rho = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                              P,
                                                                                              bc_info[f].values[1],
                                                                                              n_q(tmp_q));

                        P += rho*driver->gravity.dot(x_q[tmp_q]);
                    }

                } else if (bc_info[f].tag == PRESSURE_OUTLET || bc_info[f].tag == DAMPED_OUTLET) {
                    //face has prescribed pressure
                    P = bc_info[f].values[0];

                    //if hydrostatic flag is set, use vertical position of quadrature point to adjust pressure and density
                    if (USE_HYDROSTATIC_CORRECTION){
                        rho = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver,
                                                                                           P,
                                                                                           bc_info[f].values[1],
                                                                                           n_q(tmp_q));

                        P += rho*driver->gravity.dot(x_q[tmp_q]);
                    }

                    if (bc_info[f].tag == DAMPED_OUTLET && p.dot(normal) > 0) {
                        //only adjust if outflow
                        P = (1.0 - damping_coefficient) * P
                            + damping_coefficient * driver->fluid_material->getPressure(job, driver,
                                                                                        rho,
                                                                                        p,
                                                                                        rhoE,
                                                                                        n); //n_q(q_list[q]));
                    }
                } else if (bc_info[f].tag == ADIABATIC_WALL ||
                           bc_info[f].tag == DAMPED_WALL ||
                           bc_info[f].tag == SYMMETRIC_WALL ||
                           bc_info[f].tag == THERMAL_WALL){
                    P = driver->fluid_material->getPressure(job, driver, rho, p, rhoE, n); //n_q(q_list[q]));

                    if (bc_info[f].tag == DAMPED_WALL){
                        //add damping term to pressure
                        //estimate L
                        for (int ii=0; ii<GRID_DIM; ii++){
                            for (int jj=0; jj<GRID_DIM; jj++){
                                L_tmp(ii,jj) = (bc_info[f].vector[ii] - driver->fluid_body->p[e_minus][ii]/driver->fluid_body->rho(e_minus))
                                               /(x_f[f] - x_e[e_minus]).dot(normal)*normal[jj];
                            }
                        }
                        P -= bc_info[f].values[0]*L_tmp.trace();
                    }
                } else if (bc_info[f].tag == SUPERSONIC_INLET) {
                    n = n_q[tmp_q];
                    rho = n*bc_info[f].values[0];
                    p = rho*bc_info[f].vector;
                    P = driver->fluid_material->getPressureFromDensityAndTemperature(job,driver,rho,bc_info[f].values[1],n);
                } else if (bc_info[f].tag == SUPERSONIC_OUTLET) {
                    //fill in pressure
                    P = driver->fluid_material->getPressure(job,
                                                            driver,
                                                            rho,
                                                            p,
                                                            rhoE,
                                                            n); //n_q(tmp_q));
                }

                //fill in tmpArray
                //tmpVecArray[tmp_q] = -P_q(tmp_q) * (1.0 - n_q(tmp_q)) * normal * getQuadratureWeight(tmp_q);
                kv[tmp_q] = -P * normal * getQuadratureWeight(tmp_q);
            } else if (!BOUNDARY_ONLY && face_elements[f][0] > -1 && face_elements[f][1] > -1) {
                //not a bounding face and well defined
                //get normal and element neighbors
                normal = getFaceNormal(job, f);
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //approximate local density, momentum, and energy
                tmp_q = int_quad_count + f*qpf + q;

                rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x_q[tmp_q] - x_e[e_plus]);
                rhoE_plus = driver->fluid_body->rhoE(e_plus) + driver->fluid_body->rhoE_x[e_plus].dot(x_q[tmp_q] - x_e[e_plus]);
                p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*(x_q[tmp_q] - x_e[e_plus]);
                n_plus = driver->fluid_body->n_e(e_plus);

                //relative position to centroid of A
                if (bc_info[f].tag == PERIODIC) {
                    //by convention, centroid of A will be on other side of domain
                    //so use relative distance to B but flip normal direction
                    x = x_q[tmp_q] - x_e[e_plus];
                    x -= 2.0 * (x.dot(normal)) * normal;
                } else {
                    x = x_q(tmp_q) - x_e(e_minus);
                }

                rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                rhoE_minus = driver->fluid_body->rhoE(e_minus) + driver->fluid_body->rhoE_x[e_minus].dot(x);
                p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*(x);
                n_minus = driver->fluid_body->n_e(e_minus);

                //for consistency with energy equations, use average pressure
                P_plus = driver->fluid_material->getPressure(job,
                                                             driver,
                                                             rho_plus,
                                                             p_plus,
                                                             rhoE_plus,
                                                             n_plus);


                P_minus = driver->fluid_material->getPressure(job,
                                                              driver,
                                                              rho_minus,
                                                              p_minus,
                                                              rhoE_minus,
                                                              n_minus);

                P = 0.5*(P_plus + P_minus);

                //fill in tmpArray
                //tmpVecArray[tmp_q] = -P_q(tmp_q) * (1.0 - n_q(tmp_q)) * normal * getQuadratureWeight(tmp_q);
                kv[tmp_q] = -P * normal * getQuadratureWeight(tmp_q);
            }
        }
    }
    return;
}


void FVMGridBase::calculateElementIntegrandsForDragForce(Job *job,
                                                               FiniteVolumeDriver *driver,
                                                               KinematicVectorArray& kv,
                                                               int e_start, int e_end){

    //check that limits are within acceptable range
    if (e_start < 0 || e_end >= element_count){
        std::cerr << "ERROR! Start or end indices out of range! " << e_start << " <? 0, " << e_end << " >=? " << element_count << std::endl;
    }

    //do nothing

    return;
}

/*----------------------------------------------------------------------------*/
//functions to get drag corrections
KinematicVectorArray FVMGridBase::calculateCorrectedDragForces(Job* job,
                                                          FiniteVolumeDriver* driver,
                                                          const Eigen::VectorXd &K_n){

    //result defines force on solid by fluid
    KinematicVectorArray result = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray tmp_result = result;

    //declare temporary containers for quadrature point fields
    KinematicVectorArray tmpVecArray = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE);
    Eigen::VectorXd tmpArray = Eigen::VectorXd(int_quad_count + ext_quad_count);

    tmpVecArray.setZero();
    tmpArray.setZero();

    double rho, rhoE, n;
    KinematicVector p, normal;

    //new drag calculation
    //calculate volume weighted velocity and mass density
    KinematicVectorArray u_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray u = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    Eigen::VectorXd m_i = Eigen::VectorXd(job->grid->node_count);
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        u[e] = driver->fluid_body->p(e) / driver->fluid_body->rho(e);
    }
    if (num_threads > 1) {
        parallelMultiply(M, u, u_i, M.TRANSPOSED, true);
    } else {
        u_i = M.operate(u, M.TRANSPOSED);
    }
    for (int i=0; i<job->grid->node_count; i++){
        u_i[i] /= job->grid->nodeVolume(job, i);
    }
    if (num_threads > 1) {
        parallelMultiply(M, driver->fluid_body->rho, m_i, M.TRANSPOSED, true);
    } else {
        m_i = M.operate(driver->fluid_body->rho, M.TRANSPOSED);
    }


    //calculate drag contribution
    double rho_f;
    for (int i=0; i<job->grid->node_count; i++){
        rho_f = m_i(i)/job->grid->nodeVolume(job,i);
        result[i] = -K_n(i) * (driver->fluid_body->v_s[i] - u_i[i]);
    }

    return result;
}

Eigen::VectorXd FVMGridBase::getCorrectedDragCoefficients(Job* job, FiniteVolumeDriver* driver) {
    Eigen::VectorXd result;

    double rho, rhoE, n, n_bar;
    KinematicVector p, normal, f_d;
    double P;
    double epsilon_0, epsilon_bar;

    result = Eigen::VectorXd(job->grid->node_count);

    //calculate volume weighted velocity and mass density
    KinematicVectorArray u_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray u = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    Eigen::VectorXd m_i = Eigen::VectorXd(job->grid->node_count);
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        u[e] = driver->fluid_body->p(e) / driver->fluid_body->rho(e);
    }
    if (num_threads > 1) {
        parallelMultiply(M, u, u_i, M.TRANSPOSED, true);
    } else {
        u_i = M.operate(u,M.TRANSPOSED);
    }
    for (int i=0; i<job->grid->node_count; i++){
        u_i[i] /= job->grid->nodeVolume(job, i);
    }
    if (num_threads > 1) {
        parallelMultiply(M, driver->fluid_body->rho, m_i, M.TRANSPOSED, true);
    } else {
        m_i = M.operate(driver->fluid_body->rho,M.TRANSPOSED);
    }


    //calculate drag coefficient
    double rho_f;
    for (int i=0; i<job->grid->node_count; i++){
        rho_f = m_i(i)/job->grid->nodeVolume(job,i);
        result(i) = driver->fluid_material->getInterphaseDragCoefficient(job,
                                                                         driver,
                                                                         rho_f,
                                                                         u_i[i],
                                                                         driver->fluid_body->v_s[i],
                                                                         driver->fluid_body->n(i),
                                                                         FiniteVolumeMaterial::CORRECTED_DRAG);
    }

    return result;
}

void FVMGridBase::calculateElementIntegrandsForCorrectedDragForces(Job *job,
                                                                   FiniteVolumeDriver *driver,
                                                                   KinematicVectorArray &kv,
                                                                   const Eigen::VectorXd &K_n,
                                                                   int e_start, int e_end){

    //check that limits are within acceptable range
    if (e_start < 0 || e_end >= element_count){
        std::cerr << "ERROR! Start or end indices out of range! " << e_start << " <? 0, " << e_end << " >=? " << element_count << std::endl;
    }

    //do nothing

    return;
}

void FVMGridBase::calculateSplitIntegralInterphaseForces(Job* job,
                                                    FiniteVolumeDriver* driver,
                                                    KinematicVectorArray& f_i,
                                                    KinematicVectorArray& f_e){
    KinematicVectorArray tmp_f_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray tmp_f_e = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);

    //get buoyant force
    calculateSplitIntegralBuoyantForces(job, driver, f_i, f_e);

    //get drag force
    calculateSplitIntegralDragForces(job, driver, tmp_f_i, tmp_f_e);

    //add drag to buoyant forces
    for (int i=0; i<job->grid->node_count; i++){
        f_i[i] += tmp_f_i[i];
    }

    for (int e=0; e<driver->fluid_grid->element_count; e++){
        f_e[e] += tmp_f_e[e];
    }

    return;
}

void FVMGridBase::calculateSplitIntegralBuoyantForces(Job* job,
                                                 FiniteVolumeDriver* driver,
                                                 KinematicVectorArray& f_i,
                                                 KinematicVectorArray& f_e){
    //declare temporary containers for quadrature point fields
    KinematicVectorArray tmpVecArray = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE);
    Eigen::VectorXd tmpArray = Eigen::VectorXd(int_quad_count + ext_quad_count);
    KinematicVectorArray tmp_result = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);

    tmpVecArray.setZero();
    tmpArray.setZero();

    double rho, rhoE, n;
    KinematicVector p, normal;

    //get integral values at quadrature points
    calculateElementIntegrandsForBuoyantForce(job, driver, tmpArray, 0, element_count-1);
    calculateFaceIntegrandsForBuoyantForce(job, driver, tmpVecArray, 0, face_count-1); //DO NOT include internal faces

    //integrate expression (pg 143, nb 6) by multiplying phi_j, N_j, gradN_j
    //loop over nodes to find f_i
    if (num_threads > 1){
        //add components to result vector
        parallelMultiply(Q, tmpVecArray, f_i, Q.NORMAL, true);       //clear f_i here
        parallelMultiply(gradQ, tmpArray, f_i, gradQ.NORMAL, false); //do not clear f_i here

        //scale result vector
        for (int i = 0; i < f_i.size(); i++) {
            f_i[i] *= (1.0 - driver->fluid_body->n(i));
        }
    } else {
        tmp_result = Q * tmpVecArray;
        for (int i = 0; i < f_i.size(); i++) {
            f_i[i] = (1.0 - driver->fluid_body->n(i)) * tmp_result(i); //ensures nodes with zero porosity have zero force
        }

        tmp_result = gradQ * tmpArray;
        for (int i = 0; i < f_i.size(); i++) {
            f_i[i] += (1.0 - driver->fluid_body->n(i)) * tmp_result(i);
        }
    }

    //divide result by associated nodal volumes
    for (int i=0; i<job->grid->node_count; i++){
        f_i(i) /= job->grid->nodeVolume(job, i);
    }

    //redo this calculation but DO include internal faces
    //tmpVecArray.setZero();
    //calculateFaceIntegrandsForBuoyantForce(job, driver, tmpVecArray, 0, face_count-1, false); //DO include internal faces

    //loop over elements to find f_e
    int tmp_e = -1;
    int tmp_q = -1;
    f_e.setZero();
    for (int e=0; e<element_count; e++){
        for (int q=0; q<qpe; q++){
            //for each quadrature point in element, multiply integrand by gradphi
            // \nabla phi = -\nabla n, requires partition of unity
            f_e[e] -= tmpArray(e*qpe + q)*gradn_q[e*qpe + q];
        }
    }
    /*
    for (int f=0; f<face_count; f++){
        for (int q=0; q<qpf; q++){
            //for each quadrature point on face, multiply integrand by phi
            tmp_q = int_quad_count + f*qpf + q;
            if (q_b[int_quad_count + f*qpf + q]) {
                //determine which element is interior
                if (face_elements[f][0] > -1) {
                    //A is interior to domain
                    tmp_e = face_elements[f][0];
                } else {
                    //B is interior to domain
                    tmp_e = face_elements[f][1];
                }
                f_e[tmp_e] += tmpVecArray[tmp_q]*(1.0 - n_q(tmp_q));
            } else if (face_elements[f][0] > -1 && face_elements[f][1] > -1){
                //subtract face integral from A and add to B (according to normal direction of face)
                f_e[face_elements[f][0]] += tmpVecArray[tmp_q]*(1.0 - n_q(tmp_q));
                f_e[face_elements[f][1]] -= tmpVecArray[tmp_q]*(1.0 - n_q(tmp_q));
            }
        }
    }
     */
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        //rescale element force
        f_e[e] /= driver->fluid_grid->getElementVolume(e);
    }

    /*
    //check that everything is kosher
    KinematicVector sum_e = KinematicVector(job->JOB_TYPE);
    KinematicVector sum_j = sum_e;
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        sum_e += f_e[e] * driver->fluid_grid->getElementVolume(e);
    }
    for (int j=0; j<job->grid->node_count; j++){
        sum_j += f_i[j] * job->grid->nodeVolume(job, j);
    }
    std::cout << "Element Forces: " << sum_e[0] << ", " << sum_e[1] << " ?= Nodal Forces: " << sum_j[0] << ", " << sum_j[1] << std::endl;
    std::cout << sum_e.norm() / sum_j.norm() << std::endl;
    */


    return;
}

void FVMGridBase::calculateSplitIntegralDragForces(Job* job,
                                              FiniteVolumeDriver* driver,
                                              KinematicVectorArray& f_i,
                                              KinematicVectorArray& f_e){
    //result defines force on solid by fluid
    KinematicVectorArray Kvs_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    Eigen::VectorXd K_i = Eigen::VectorXd(job->grid->node_count);

    //new drag calculation
    //calculate volume weighted velocity and mass density on nodes to estimate K_i
    KinematicVectorArray u_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray u = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    Eigen::VectorXd m_i = Eigen::VectorXd(job->grid->node_count);
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        u[e] = driver->fluid_body->p(e) / driver->fluid_body->rho(e);
    }
    if (num_threads > 1) {
        parallelMultiply(M, u, u_i, M.TRANSPOSED, true);
    } else {
        u_i = M.operate(u, M.TRANSPOSED);
    }
    for (int i=0; i<job->grid->node_count; i++){
        u_i[i] /= job->grid->nodeVolume(job, i);
    }

    if (num_threads > 1) {
        parallelMultiply(M, driver->fluid_body->rho, m_i, M.TRANSPOSED, true);
    } else {
        m_i = M.operate(driver->fluid_body->rho, M.TRANSPOSED);
    }

    //calculate drag coefficient associated with each node
    double rho_f;
    for (int i=0; i<job->grid->node_count; i++){
        rho_f = m_i(i)/job->grid->nodeVolume(job,i);
        K_i(i) = driver->fluid_material->getInterphaseDragCoefficient(job,
                                                             driver,
                                                             rho_f,
                                                             u_i[i],
                                                             driver->fluid_body->v_s[i],
                                                             driver->fluid_body->n(i),
                                                             FiniteVolumeMaterial::REGULAR_DRAG);

        //estimate nodal drag
        f_i[i] = -K_i(i)*(driver->fluid_body->v_s[i] - u_i[i]);
        Kvs_i[i] = K_i(i) * driver->fluid_body->v_s[i];
    }

    //calculate element-wise weighted drag coefficient
    Eigen::VectorXd K_e = Eigen::VectorXd(driver->fluid_grid->element_count);
    if (num_threads > 1){
        parallelMultiply(M, K_i, K_e, M.NORMAL, true);
    } else {
        K_e = M*K_i;
    }

    //calculate element-wise weighted solid drag component velocity
    KinematicVectorArray Kvs_e = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    if (num_threads > 1){
        parallelMultiply(M, Kvs_i, Kvs_e, M.NORMAL, true);
    } else {
        Kvs_e = M*Kvs_i;
    }

    //calculate element estimate for drag
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        f_e[e] = -(Kvs_e[e] - K_e(e)*u[e]) / driver->fluid_grid->getElementVolume(e);
    }

    return;
}

void FVMGridBase::calculateSplitIntegralCorrectedDragForces(Job* job,
                                                       FiniteVolumeDriver* driver,
                                                       KinematicVectorArray& f_i,
                                                       KinematicVectorArray& f_e,
                                                       const Eigen::VectorXd &K_n){
    //result defines force on solid by fluid
    KinematicVectorArray Kvs_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);

    //new drag calculation
    //calculate volume weighted velocity and mass density on nodes to estimate K_i
    KinematicVectorArray u_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray u = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    Eigen::VectorXd m_i = Eigen::VectorXd(job->grid->node_count);
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        u[e] = driver->fluid_body->p(e) / driver->fluid_body->rho(e);
    }
    if (num_threads > 1) {
        parallelMultiply(M, u, u_i, M.TRANSPOSED, true);
    } else {
        u_i = M.operate(u, M.TRANSPOSED);
    }
    for (int i=0; i<job->grid->node_count; i++){
        u_i[i] /= job->grid->nodeVolume(job, i);
    }

    if (num_threads > 1) {
        parallelMultiply(M, driver->fluid_body->rho, m_i, M.TRANSPOSED, true);
    } else {
        m_i = M.operate(driver->fluid_body->rho, M.TRANSPOSED);
    }

    //calculate drag coefficient associated with each node
    double rho_f;
    for (int i=0; i<job->grid->node_count; i++){
        rho_f = m_i(i)/job->grid->nodeVolume(job,i);

        //estimate nodal drag
        f_i[i] = -K_n(i)*(driver->fluid_body->v_s[i] - u_i[i]);
        Kvs_i[i] = K_n(i) * driver->fluid_body->v_s[i];
    }

    //calculate element-wise weighted drag coefficient
    Eigen::VectorXd K_e = Eigen::VectorXd(driver->fluid_grid->element_count);
    if (num_threads > 1){
        parallelMultiply(M, K_n, K_e, M.NORMAL, true);
    } else {
        K_e = M*K_n;
    }

    //calculate element-wise weighted solid drag component velocity
    KinematicVectorArray Kvs_e = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    if (num_threads > 1){
        parallelMultiply(M, Kvs_i, Kvs_e, M.NORMAL, true);
    } else {
        Kvs_e = M*Kvs_i;
    }

    //calculate element estimate for drag
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        f_e[e] = -(Kvs_e[e] - K_e(e)*u[e]) / driver->fluid_grid->getElementVolume(e);
    }

    return;
}


/*----------------------------------------------------------------------------*/
//energy flux terms due to interphase forces
Eigen::VectorXd FVMGridBase::calculateInterphaseEnergyFlux(Job* job, FiniteVolumeDriver* driver){
    //result defines interphase energy flux into fluid elements
    Eigen::VectorXd result = Eigen::VectorXd(element_count);
    result.setZero();

    //don't bother for ISOTHERMAL sims
    if (driver->TYPE != driver->THERMAL){
        return result;
    }

    /*-------------------------*/
    /* BEGIN DRAG CONTRIBUTION */
    /*-------------------------*/

    KinematicVectorArray Kvs_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    Eigen::VectorXd K_i = Eigen::VectorXd(job->grid->node_count);
    Eigen::VectorXd Kvsvs_i = Eigen::VectorXd(job->grid->node_count);

    //new drag calculation
    //calculate volume weighted velocity and mass density on nodes to estimate K_i
    KinematicVectorArray u_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    KinematicVectorArray u = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    Eigen::VectorXd m_i = Eigen::VectorXd(job->grid->node_count);
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        u[e] = driver->fluid_body->p(e) / driver->fluid_body->rho(e);
    }
    if (num_threads > 1) {
        parallelMultiply(M, u, u_i, M.TRANSPOSED, true);
    } else {
        u_i = M.operate(u, M.TRANSPOSED);
    }
    for (int i=0; i<job->grid->node_count; i++){
        u_i[i] /= job->grid->nodeVolume(job, i);
    }

    if (num_threads > 1) {
        parallelMultiply(M, driver->fluid_body->rho, m_i, M.TRANSPOSED, true);
    } else {
        m_i = M.operate(driver->fluid_body->rho, M.TRANSPOSED);
    }

    //calculate drag coefficient associated with each node
    double rho_f;
    for (int i=0; i<job->grid->node_count; i++){
        rho_f = m_i(i)/job->grid->nodeVolume(job,i);
        K_i(i) = driver->fluid_material->getInterphaseDragCoefficient(job,
                                                                      driver,
                                                                      rho_f,
                                                                      u_i[i],
                                                                      driver->fluid_body->v_s[i],
                                                                      driver->fluid_body->n(i),
                                                                      FiniteVolumeMaterial::REGULAR_DRAG);

        //intermediate calculation
        Kvs_i[i] = K_i(i) * driver->fluid_body->v_s[i];
        Kvsvs_i(i) = K_i(i) * driver->fluid_body->v_s[i].dot(driver->fluid_body->v_s[i]);
    }

    //first term in drag component
    if (num_threads > 1){
        parallelMultiply(M, Kvsvs_i, result, M.NORMAL, true);
    } else {
        result = M*Kvsvs_i;
    }

    //calculate element-wise weighted solid drag coefficient
    //Eigen::VectorXd K_e = Eigen::VectorXd(driver->fluid_grid->element_count);
    //if (num_threads > 1){
    //    parallelMultiply(M, K_i, K_e, M.NORMAL, true);
    //} else {
    //    K_e = M*K_i;
    //}

    //calculate element-wise weighted solid drag component velocity
    KinematicVectorArray Kvs_e = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    if (num_threads > 1){
        parallelMultiply(M, Kvs_i, Kvs_e, M.NORMAL, true);
    } else {
        Kvs_e = M*Kvs_i;
    }

    //add drag term
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        result(e) -= Kvs_e[e].dot(u[e]);
        //result(e) = Kvs_e[e].dot(u[e]) - K_e(e)*u[e].dot(u[e]);
    }

    /*-----------------------------*/
    /* BEGIN PRESSURE CONTRIBUTION */
    /*-----------------------------*/

    if (num_threads > 1) {
        //determine number of threads for element flux evaluation
        int thread_count = 0;
        if (element_count >= num_threads) {
            thread_count = num_threads;
        } else {
            thread_count = element_count;
        }

        //intermediate storage vector
        //check that memUnitID is valid
        if (memoryUnits[0].s.size() != thread_count) {
            memoryUnits[0].s.resize(thread_count);
        }

        //boolean of completion status
        //volatile bool firstTaskComplete[thread_count] = {false};
        volatile bool *firstTaskComplete = new volatile bool[thread_count];
        for (int t=0; t<thread_count; t++){
            firstTaskComplete[t] = false;
        }

        //choose interval size
        int k_max = element_count - 1;
        int k_interval = (element_count / thread_count) + 1;
        int k_begin, k_end;

        for (int t = 0; t < thread_count; t++) {
            //initialize output vector
            if (memoryUnits[0].s[t].size() != element_count) {
                memoryUnits[0].s[t] = Eigen::VectorXd(element_count);
            }

            //zero output vector
            memoryUnits[0].s[t].setZero();

            //set interval
            k_begin = t * k_interval;
            k_end = k_begin + k_interval - 1;
            if (k_end > k_max) {
                k_end = k_max;
            }
            //send job to threadpool
            jobThreadPool->doJob(std::bind(calcElementIntegrandsForInterphaseEnergyFlux,
                                           job,
                                           driver,
                                           this,
                                           std::ref(memoryUnits[0].s[t]),
                                           k_begin,
                                           k_end,
                                           std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone) {
            //set flag to true
            taskDone = true;
            for (int t = 0; t < thread_count; t++) {
                if (!firstTaskComplete[t]) {
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        //determine number of threads for addition
        if (element_count > num_threads) {
            thread_count = num_threads;
        } else {
            thread_count = element_count;
        }

        //boolean of completion status
        //volatile bool secondTaskComplete[thread_count] = {false};
        volatile bool *secondTaskComplete = new volatile bool[thread_count];
        for (int t=0; t<thread_count; t++){
            secondTaskComplete[t] = false;
        }

        //choose interval size
        int i_max = element_count - 1;
        int i_interval = (element_count / thread_count) + 1;
        int i_begin, i_end;

        for (int t = 0; t < thread_count; t++) {
            //set interval
            i_begin = t * i_interval;
            i_end = i_begin + i_interval - 1;
            if (i_end > i_max) {
                i_end = i_max;
            }
            //send job to thread pool
            jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::scalarAddwithFlag,
                                           std::ref(memoryUnits[0].s),
                                           std::ref(result),
                                           i_begin,
                                           i_end,
                                           false,
                                           std::ref(secondTaskComplete[t])));
        }

        //join threads
        taskDone = false;
        //wait for task to complete
        while (!taskDone) {
            //set flag to true
            taskDone = true;
            for (int t = 0; t < thread_count; t++) {
                if (!secondTaskComplete[t]) {
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        delete[] firstTaskComplete;
        delete[] secondTaskComplete;

        /*

        //determine number of threads for face flux evaluation
        if (face_count >= num_threads){
            thread_count = num_threads;
        } else {
            thread_count = face_count;
        }

        //intermediate storage vector
        //check that memUnitID is valid
        if (memoryUnits[0].s.size() != thread_count){
            memoryUnits[0].s.resize(thread_count);
        }

        //boolean of completion status
        volatile bool thirdTaskComplete[thread_count] = {false};

        //choose interval size
        k_max = face_count - 1;
        k_interval = (face_count/thread_count) + 1;

        for (int t=0; t<thread_count; t++) {
            //initialize output vector
            if (memoryUnits[0].s[t].size() != element_count){
                memoryUnits[0].s[t] = Eigen::VectorXd(element_count);
            }

            //zero output vector
            memoryUnits[0].s[t].setZero();

            //set interval
            k_begin = t * k_interval;
            k_end = k_begin + k_interval - 1;
            if (k_end > k_max){
                k_end = k_max;
            }
            //send job to threadpool
            jobThreadPool->doJob(std::bind(calcFaceIntegrandsForInterphaseEnergyFlux,
                                           job,
                                           driver,
                                           this,
                                           std::ref(memoryUnits[0].s[t]),
                                           k_begin,
                                           k_end,
                                           std::ref(thirdTaskComplete[t])));
        }

        //join threads
        taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!thirdTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        //determine number of threads for addition
        if (element_count > num_threads){
            thread_count = num_threads;
        } else {
            thread_count = element_count;
        }

        //boolean of completion status
        volatile bool fourthTaskComplete[thread_count] = {false};

        //choose interval size
        i_max = element_count - 1;
        i_interval = (element_count/thread_count) + 1;

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //send job to thread pool
            jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::scalarAddwithFlag,
                                           std::ref(memoryUnits[0].s),
                                           std::ref(result),
                                           i_begin,
                                           i_end,
                                           false,
                                           std::ref(fourthTaskComplete[t])));
        }

        //join threads
        taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!fourthTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
         */
    } else {
        //do calculation serially
        calculateElementIntegrandsForInterphaseEnergyFlux(job, driver, result, 0, element_count - 1);
        calculateFaceIntegrandsForInterphaseEnergyFlux(job, driver, result, 0, face_count - 1);
    }

    return result;
}

Eigen::VectorXd FVMGridBase::calculateInterphaseEnergyFluxUsingElementBasedForce(Job* job,
                                                                                 FiniteVolumeDriver* driver,
                                                                                 const KinematicVectorArray &f_e){
    std::cerr << "ERROR: calculateInterphaseEnergyFluxUsingElementBasedForce() is deprecated! Exiting." << std::endl;
    exit(0);

    Eigen::VectorXd result = Eigen::VectorXd(element_count);
    return result;
}

Eigen::VectorXd FVMGridBase::calculateInterphaseEnergyFluxUsingNodeBasedDrag(Job* job,
                                                                             FiniteVolumeDriver* driver,
                                                                             const KinematicVectorArray &f_d){
    std::cerr << "ERROR: FVMGridBase::calculateInterphseEnergyFluxUsingNodeBasedDrag() not implemented! Exiting." << std::endl;
    exit(0);

    return Eigen::VectorXd(0);
}

Eigen::VectorXd FVMGridBase::calculateElementIntegrandsForInterphaseEnergyFlux(Job* job,
                                                                                FiniteVolumeDriver* driver,
                                                                               Eigen::VectorXd &result,
                                                                                int e_start, int e_end){

    //check that limits are within acceptable range
    if (e_start < 0 || e_end >= element_count){
        std::cerr << "ERROR! Start or end indices out of range! " << e_start << " <? 0, " << e_end << " >=? " << element_count << std::endl;
    }

    double rho, rhoE, n, P, trL;
    KinematicVector p, u;
    double epsilon_0, epsilon_bar;

    //loop over elements
    for (int e = e_start; e <= e_end; e++) {
        for (int q = 0; q < qpe; q++) {
            //approximate local density, momentum, and energy
            rho = driver->fluid_body->rho(e) + driver->fluid_body->rho_x[e].dot(x_q[e * qpe + q] - x_e[e]);
            rhoE = driver->fluid_body->rhoE(e) + driver->fluid_body->rhoE_x[e].dot(x_q[e * qpe + q] - x_e[e]);
            p = driver->fluid_body->p[e] + driver->fluid_body->p_x[e] * (x_q[e * qpe + q] - x_e[e]);
            n = driver->fluid_body->n_e(e);

            //estimate local pressure
            P = driver->fluid_material->getPressure(job, driver, rho, p, rhoE, n);

            //fill in pressure term contribution to energy flux
            u = p/rho;
            //result(e) += getQuadratureWeight(e*qpe + q) * P * u.dot(gradn_q[e*qpe + q]);
            result(e) += getQuadratureWeight(e*qpe + q) * P * v_sq[e*qpe + q].dot(gradn_q[e*qpe + q]);

            //CORRECTION pg 42, nb #9
            if (USE_CORRECT_INTERPHASE_ENERGY_FLUX){
                result(e) -= getQuadratureWeight(e*qpe + q) * P * (1.0 - n_q(e*qpe + q)) * div_v_sq(e*qpe + q);
            }
        }
    }

    return result;

}

Eigen::VectorXd FVMGridBase::calculateFaceIntegrandsForInterphaseEnergyFlux(Job* job,
                                                                            FiniteVolumeDriver* driver,
                                                                            Eigen::VectorXd &result,
                                                                            int f_start, int f_end){

    //check that limits are within acceptable range
    if (f_start < 0 || f_end >= face_count){
        std::cerr << "ERROR! Start or end indices out of range! " << f_start << " <? 0, " << f_end << " >=? " << element_count << std::endl;
    }

    //do nothing, this has already been accounted for in energy fluxes above

    return result;

}


/*----------------------------------------------------------------------------*/
//parallel functions
void FVMGridBase::parallelMultiply(const MPMScalarSparseMatrix &S,
                                   const Eigen::VectorXd &x,
                                   Eigen::VectorXd &lhs,
                                   int SPEC,
                                   bool clear,
                                   int memUnitID){
    //get length of sparse matrix storage
    int k_max = S.size() - 1;

    //get length of output vector
    int i_max = lhs.rows() - 1;

    //determine number of threads for matrix vector mult
    int thread_count = 0;
    if (S.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = S.size();
    }

    //intermediate storage vector
    //check that memUnitID is valid
    if (memUnitID >= memoryUnits.size()){
        std::cerr << "ERROR! Invalid memUnitID passed to parallelMultiply! " << memUnitID << " >= " << memoryUnits.size() << std::endl;
    }
    if (memoryUnits[memUnitID].s.size() != thread_count){
        memoryUnits[memUnitID].s.resize(thread_count);
    }

    //boolean of completion status
    //volatile bool firstTaskComplete[thread_count] = {false};
    volatile bool *firstTaskComplete = new volatile bool[thread_count];
    for (int t=0; t<thread_count; t++){
        firstTaskComplete[t] = false;
    }

    //choose interval size
    int k_interval = (S.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        if (memoryUnits[memUnitID].s[t].rows() != lhs.rows()){
            memoryUnits[memUnitID].s[t] = Eigen::VectorXd(lhs.rows());
        }

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //send job to threadpool
        jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::ssmOperateStoSwithFlag,
                                       std::ref(S),
                                       std::ref(x),
                                       std::ref(memoryUnits[memUnitID].s[t]),
                                       SPEC,
                                       k_begin,
                                       k_end,
                                       std::ref(firstTaskComplete[t])));
    }

    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            //std::cout << firstTaskComplete[t] << "?" << std::endl;
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //determine number of threads for addition
    if (lhs.rows() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.rows();
    }

    //volatile bool secondTaskComplete[thread_count] = {false};
    volatile bool *secondTaskComplete = new volatile bool[thread_count];
    for (int t=0; t<thread_count; t++){
        secondTaskComplete[t] = false;
    }

    //choose interval size
    int i_interval = (lhs.rows()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //send job to threadpool
        jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::scalarAddwithFlag,
                                       std::ref(memoryUnits[memUnitID].s),
                                       std::ref(lhs),
                                       i_begin,
                                       i_end,
                                       clear,
                                       std::ref(secondTaskComplete[t])));
    }

    //join threads
    taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    delete[] firstTaskComplete;
    delete[] secondTaskComplete;

    //should be all done
    return;
}

void FVMGridBase::parallelMultiply(const MPMScalarSparseMatrix &S,
                                   const KinematicVectorArray &x,
                                   KinematicVectorArray &lhs,
                                   int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = S.size() - 1;

    //get length of output vector
    int i_max = lhs.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count = 0;
    if (S.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = S.size();
    }

    //intermediate storage vector
    //check that memUnitID is valid
    if (memUnitID >= memoryUnits.size()){
        std::cerr << "ERROR! Invalid memUnitID passed to parallelMultiply! " << memUnitID << " >= " << memoryUnits.size() << std::endl;
    }
    if (memoryUnits[memUnitID].kv.size() != thread_count){
        memoryUnits[memUnitID].kv.resize(thread_count);
    }

    //boolean of completion status
    //volatile bool firstTaskComplete[thread_count] = {false};
    volatile bool *firstTaskComplete = new volatile bool[thread_count];
    for (int t=0; t<thread_count; t++){
        firstTaskComplete[t] = false;
    }

    //choose interval size
    int k_interval = (S.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        //lhs_vec[t] = KinematicVectorArray(lhs.size(),lhs.VECTOR_TYPE);
        if (memoryUnits[memUnitID].kv[t].size() != lhs.size()){
            memoryUnits[memUnitID].kv[t] = KinematicVectorArray(lhs.size(),lhs.VECTOR_TYPE);
        }

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //send job to threadpool
        jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::ssmOperateVtoVwithFlag,
                                       std::ref(S),
                                       std::ref(x),
                                       std::ref(memoryUnits[memUnitID].kv[t]),
                                       SPEC,
                                       k_begin,
                                       k_end,
                                       std::ref(firstTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //determine number of threads for addition
    if (lhs.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.size();
    }

    //boolean of completion status
    //volatile bool secondTaskComplete[thread_count] = {false};
    volatile bool *secondTaskComplete = new volatile bool[thread_count];
    for (int t=0; t<thread_count; t++){
        secondTaskComplete[t] = false;
    }

    //choose interval size
    int i_interval = (lhs.size()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //send job to thread pool
        jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::vectorAddKwithFlag,
                                       std::ref(memoryUnits[memUnitID].kv),
                                       std::ref(lhs),
                                       i_begin,
                                       i_end,
                                       clear,
                                       std::ref(secondTaskComplete[t])));
    }

    //join threads
    taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    delete[] firstTaskComplete;
    delete[] secondTaskComplete;

    //should be all done
    return;
}

void FVMGridBase::parallelMultiply(const KinematicVectorSparseMatrix &gradS,
                                   const Eigen::VectorXd &x,
                                   KinematicVectorArray &lhs,
                                   int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = gradS.size() - 1;

    //get length of output vector
    int i_max = lhs.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count = 0;
    if (gradS.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = gradS.size();
    }

    //intermediate storage vector
    //check that memUnitID is valid
    if (memUnitID >= memoryUnits.size()){
        std::cerr << "ERROR! Invalid memUnitID passed to parallelMultiply! " << memUnitID << " >= " << memoryUnits.size() << std::endl;
    }
    if (memoryUnits[memUnitID].kv.size() != thread_count){
        memoryUnits[memUnitID].kv.resize(thread_count);
    }

    //boolean of completion status
    //volatile bool firstTaskComplete[thread_count] = {false};
    volatile bool *firstTaskComplete = new volatile bool[thread_count];
    for (int t=0; t<thread_count; t++){
        firstTaskComplete[t] = false;
    }

    //choose interval size
    int k_interval = (gradS.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        //lhs_vec[t] = KinematicVectorArray(lhs.size(),lhs.VECTOR_TYPE);
        if (memoryUnits[memUnitID].kv[t].size() != lhs.size()){
            memoryUnits[memUnitID].kv[t] = KinematicVectorArray(lhs.size(),lhs.VECTOR_TYPE);
        }

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //send job to threadpool
        jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::kvsmOperateStoVwithFlag,
                                       std::ref(gradS),
                                       std::ref(x),
                                       std::ref(memoryUnits[memUnitID].kv[t]),
                                       SPEC,
                                       k_begin,
                                       k_end,
                                       std::ref(firstTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //determine number of threads for addition
    if (lhs.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.size();
    }

    //boolean of completion status
    //volatile bool secondTaskComplete[thread_count] = {false};
    volatile bool *secondTaskComplete = new volatile bool[thread_count];
    for (int t=0; t<thread_count; t++){
        secondTaskComplete[t] = false;
    }

    //choose interval size
    int i_interval = (lhs.size()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //send job to thread pool
        jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::vectorAddKwithFlag,
                                       std::ref(memoryUnits[memUnitID].kv),
                                       std::ref(lhs),
                                       i_begin,
                                       i_end,
                                       clear,
                                       std::ref(secondTaskComplete[t])));
    }

    //join threads
    taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    delete[] firstTaskComplete;
    delete[] secondTaskComplete;

    //should be all done
    return;
}

void FVMGridBase::parallelMultiply(const KinematicVectorSparseMatrix &gradS,
                                   const KinematicVectorArray &v,
                                   Eigen::VectorXd &lhs,
                                   int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = gradS.size() - 1;

    //get length of output vector
    int i_max = lhs.rows() - 1;

    //determine number of threads for matrix vector mult
    int thread_count = 0;
    if (gradS.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = gradS.size();
    }

    //intermediate storage vector
    //check that memUnitID is valid
    if (memUnitID >= memoryUnits.size()){
        std::cerr << "ERROR! Invalid memUnitID passed to parallelMultiply! " << memUnitID << " >= " << memoryUnits.size() << std::endl;
    }
    if (memoryUnits[memUnitID].s.size() != thread_count){
        memoryUnits[memUnitID].s.resize(thread_count);
    }

    //boolean of completion status
    //volatile bool firstTaskComplete[thread_count] = {false};
    volatile bool *firstTaskComplete = new volatile bool[thread_count];
    for (int t=0; t<thread_count; t++){
        firstTaskComplete[t] = false;
    }

    //choose interval size
    int k_interval = (gradS.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        //lhs_vec[t] = KinematicVectorArray(lhs.size(),lhs.VECTOR_TYPE);
        if (memoryUnits[memUnitID].s[t].size() != lhs.rows()){
            memoryUnits[memUnitID].s[t] = Eigen::VectorXd(lhs.rows());
        }

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //send job to threadpool
        jobThreadPool->doJob(std::bind(kvsmOperateVtoS,
                                       std::ref(gradS),
                                       std::ref(v),
                                       std::ref(memoryUnits[memUnitID].s[t]),
                                       SPEC,
                                       k_begin,
                                       k_end,
                                       std::ref(firstTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //determine number of threads for addition
    if (lhs.rows() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.rows();
    }

    //boolean of completion status
    //volatile bool secondTaskComplete[thread_count] = {false};
    volatile bool *secondTaskComplete = new volatile bool[thread_count];
    for (int t=0; t<thread_count; t++){
        secondTaskComplete[t] = false;
    }

    //choose interval size
    int i_interval = (lhs.rows()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //send job to thread pool
        jobThreadPool->doJob(std::bind(ThreadPoolExplicitUSL::scalarAddwithFlag,
                                       std::ref(memoryUnits[memUnitID].s),
                                       std::ref(lhs),
                                       i_begin,
                                       i_end,
                                       clear,
                                       std::ref(secondTaskComplete[t])));
    }

    //join threads
    taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    delete[] firstTaskComplete;
    delete[] secondTaskComplete;

    //should be all done
    return;
}

//kinematic vecotr sparse matrix multiplying scalar vector
void FVMGridBase::kvsmOperateVtoS(const KinematicVectorSparseMatrix &gradS,
                                          const KinematicVectorArray &v,
                                          Eigen::VectorXd &lhs,
                                          int SPEC,
                                          int k_begin, int k_end, volatile bool &done){
    assert(v.size() == gradS.cols(SPEC) && v.VECTOR_TYPE == gradS.VECTOR_TYPE && "Vector sparse matrix multiplication (dot product) failed.");
    lhs.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = gradS.get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = gradS.get_j_index_ref(SPEC);

    for (int k=k_begin;k<=k_end;k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        lhs(i_tmp) += gradS.buffer[k].dot(v[j_tmp]);
    }

    done = true;

    return;
}


/*----------------------------------------------------------------------------*/
//functions for parallel calculation of element fluxes by face
void FVMGridBase::calcMassFluxes(Job *job,
                           FiniteVolumeDriver *driver,
                           FVMGridBase *grid,
                           Eigen::VectorXd &lhs,
                           int k_begin, int k_end,
                           volatile bool &done){

    //be careful that lhs is only passed here
    lhs = grid->calculateElementMassFluxes(job, driver, k_begin, k_end);

    //let everyone know we are done
    done = true;

    return;
}

void FVMGridBase::calcMomentumFluxes(Job *job,
                               FiniteVolumeDriver *driver,
                               FVMGridBase *grid,
                               KinematicVectorArray &lhs,
                               int k_begin, int k_end,
                               volatile bool &done){

    //be careful that lhs is only passed here
    lhs = grid->calculateElementMomentumFluxes(job, driver, k_begin, k_end);

    //let everyone know we are done
    done = true;

    return;
}

void FVMGridBase::calcEnergyFluxes(Job *job,
                             FiniteVolumeDriver *driver,
                             FVMGridBase *grid,
                             Eigen::VectorXd &lhs,
                             int k_begin, int k_end,
                             volatile bool &done){

    //be careful that lhs is only passed here
    lhs = grid->calculateElementEnergyFluxes(job, driver, k_begin, k_end);

    //let everyone know we are done
    done = true;

    return;
}


void FVMGridBase::calcElementIntegrandForInterphaseForces(Job *job,
                                                    FiniteVolumeDriver *driver,
                                                    FVMGridBase *grid,
                                                    KinematicVectorArray &kv,
                                                    Eigen::VectorXd& v,
                                                    int k_begin, int k_end,
                                                    volatile bool &done){

    //be careful that lhs is only passed here
    grid->calculateElementIntegrandsForInterphaseForce(job, driver, kv, v, k_begin, k_end);

    //let everyone know we are done
    done = true;

    return;
}


void FVMGridBase::calcFaceIntegrandForInterphaseForces(Job *job,
                                                 FiniteVolumeDriver *driver,
                                                 FVMGridBase *grid,
                                                 KinematicVectorArray &kv,
                                                 Eigen::VectorXd& v,
                                                 int k_begin, int k_end,
                                                 volatile bool &done){

    //be careful that lhs is only passed here
    grid->calculateFaceIntegrandsForInterphaseForce(job, driver, kv, v, k_begin, k_end);

    //let everyone know we are done
    done = true;

    return;
}

void FVMGridBase::calcElementIntegrandsForInterphaseEnergyFlux(Job* job,
                                                         FiniteVolumeDriver* driver,
                                                         FVMGridBase* grid,
                                                         Eigen::VectorXd& v,
                                                         int k_begin, int k_end,
                                                         volatile bool &done){

    //careful that v is unique!
    grid->calculateElementIntegrandsForInterphaseEnergyFlux(job, driver, v, k_begin, k_end);

    //let everyone know we are done
    done = true;

    return;
}


void FVMGridBase::calcFaceIntegrandsForInterphaseEnergyFlux(Job* job,
                                                      FiniteVolumeDriver* driver,
                                                      FVMGridBase* grid,
                                                      Eigen::VectorXd& v,
                                                      int k_begin, int k_end,
                                                      volatile bool &done){
    //careful that v is unique!
    grid->calculateFaceIntegrandsForInterphaseEnergyFlux(job, driver, v, k_begin, k_end);

    //let everyone know we are done
    done = true;

    return;
}


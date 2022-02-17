//
// Created by aaron on 2/24/20.
// fvm_steady_state_solver.cpp
//

//
// Created by aaron on 12/23/19.
// fvm_default_solver.cpp
//


#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <job.hpp>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "mpm_objects.hpp"
#include "fvm_objects.hpp"
#include "fvm_solvers.hpp"

void FVMSteadyStateSolver::init(Job* job, FiniteVolumeDriver* driver){
    //size flux containers
    density_fluxes = Eigen::VectorXd(driver->fluid_grid->element_count);
    momentum_fluxes = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    energy_fluxes = Eigen::VectorXd(driver->fluid_grid->element_count);

    //size temporary containers
    tmp_df = Eigen::VectorXd(driver->fluid_grid->element_count);
    tmp_mf = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    tmp_ef = Eigen::VectorXd(driver->fluid_grid->element_count);

    //fp64 properties defined, overwrite associated variables
    //nu
    if (fp64_props.size() > 0){
        nu = fp64_props[0];
    }

    //h
    if (fp64_props.size() > 1){
        h = fp64_props[1];
    }

    //REL_TOL
    if (fp64_props.size() > 2){
        REL_TOL = fp64_props[2];
    }

    //ABS_TOL
    if (fp64_props.size() > 3){
        ABS_TOL = fp64_props[3];
    }

    //LineSearch_TOL
    if (fp64_props.size() > 4){
        LineSearch_TOL = fp64_props[4];
    } else {
        LineSearch_TOL = REL_TOL;
    }

    //BiCGSTAB_TOL
    if (fp64_props.size() > 5){
        BiCGSTAB_TOL = fp64_props[5];
    }

    //INNER_SOLVER
    if (int_props.size() > 0){
        INNER_SOLVER = int_props[0];
        if (INNER_SOLVER != BICGSTAB &&
            INNER_SOLVER != DIRECT &&
            INNER_SOLVER != JACOBI){
            std::cerr << "ERROR! FVMSteadyStateSolver does not have INNER_SOLVER defined for input " << INNER_SOLVER << "! Exiting." << std::endl;
            exit(0);
        }
    }

    //max_s (limit iterations of continuation steps)
    if (int_props.size() > 1){
        max_s = int_props[1];
    }

    //max_n (limit iterations of outer loop)
    if (int_props.size() > 2){
        max_n = int_props[2];
    }

    //max_i (limit iterations of inner loop)
    if (int_props.size() > 3){
        max_i = int_props[3];
    }

    //size of vectors
    if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL){
        //rho, rho*u
        vector_size = driver->fluid_grid->element_count*(1 + driver->fluid_grid->GRID_DIM);
    } else if (driver->TYPE == FiniteVolumeDriver::THERMAL){
        //rho, rho*u, rho*E
        vector_size = driver->fluid_grid->element_count*(2 + driver->fluid_grid->GRID_DIM);
    } else {
        std::cerr << "ERROR! FVMSteadyStateSolver not implemented for simulationg TYPE " << driver->TYPE << "! Exiting." << std::endl;
        exit(0);
    }

    //initialize u_0
    convertStateSpaceToVector(job, driver,
                              u_0,
                              driver->fluid_body->rho,
                              driver->fluid_body->p,
                              driver->fluid_body->rhoE);

    std::cout << "FVMSteadyStateSolver properties:" << std::endl;
    std::cout << "   nu      = " << nu << std::endl;
    std::cout << "   h       = " << h << std::endl;
    std::cout << "   REL_TOL = " << REL_TOL << std::endl;
    std::cout << "   ABS_TOL = " << ABS_TOL << std::endl;
    std::cout << "   LineSearch_TOL = " << LineSearch_TOL << std::endl;
    std::cout << "   BiCGSTAB_TOL   = " << BiCGSTAB_TOL << std::endl;
    std::cout << "   max_i   = " << max_i << std::endl;
    std::cout << "   max_n   = " << max_n << std::endl;
    std::cout << "   max_s   = " << max_s << std::endl;

    //done
    std::cout << "FiniteVolumeSolver initialized." << std::endl;
    return;
}


//flux functions
Eigen::VectorXd FVMSteadyStateSolver::F(Job* job, FiniteVolumeDriver* driver,
                                        const Eigen::VectorXd& u,
                                        double nu_i){
    //convert input vector into fluid_body state space
    convertVectorToStateSpace(job, driver,
                              u,
                              driver->fluid_body->rho,
                              driver->fluid_body->p,
                              driver->fluid_body->rhoE);

    //calculate fluxes
    driver->fluid_grid->constructDensityField(job, driver);
    driver->fluid_grid->constructMomentumField(job, driver);
    driver->fluid_grid->constructEnergyField(job, driver);
    driver->fluid_grid->constructPorosityField(job, driver);
    driver->fluid_grid->constructVelocityField(job, driver);


    //fluid fluxes associated with each volume
    density_fluxes = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
    momentum_fluxes = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2
    energy_fluxes = driver->fluid_grid->calculateElementEnergyFluxes(job, driver);      //J/s?

    //add gravity and scale by element volume
    double volume;
    for (int e = 0; e < driver->fluid_grid->element_count; e++) {
        volume = driver->fluid_grid->getElementVolume(e);
        density_fluxes(e) /= volume;
        momentum_fluxes(e) /= volume;
        energy_fluxes(e) /= volume;
        momentum_fluxes[e] += driver->fluid_body->rho(e) * driver->gravity;
    }

    //convert resulting fluxes to output
    Eigen::VectorXd result = Eigen::VectorXd(vector_size);
    convertStateSpaceToVector(job, driver, result, density_fluxes, momentum_fluxes, energy_fluxes);

    //if nu > 0, adjust residual
    if (nu_i > 0) {
        for (int j=0; j<vector_size; j++) {
            result(j) += nu*(u_0(j) - u(j));
        }
    }

    return result;
}

Eigen::VectorXd FVMSteadyStateSolver::DF(Job* job, FiniteVolumeDriver* driver, const Eigen::VectorXd& x_i, double h_i){
    //check that input is non-zero
    if (x_i.norm() < ABS_TOL){
        return Eigen::VectorXd::Zero(vector_size);
    }

    //initialize step in x-direction
    Eigen::VectorXd tmp_plus = u_n + h_i*x_i/x_i.norm();
    Eigen::VectorXd tmp_minus = u_n - h_i*x_i/x_i.norm();

    //calculate numerical derivative
    //tmp_plus = (F(job, driver, tmp_plus) - r_n) / h_i;
    tmp_plus = (F(job, driver, tmp_plus, nu) - F(job, driver, tmp_minus, nu)) / (2*h_i);

    //return result
    return tmp_plus;
}

void FVMSteadyStateSolver::convertVectorToStateSpace(Job* job, FiniteVolumeDriver* driver,
                                                     const Eigen::VectorXd& v,
                               Eigen::VectorXd& rho,
                               KinematicVectorArray& p,
                               Eigen::VectorXd& rhoE){
    //check vector for correct dimensions
    if (v.rows() != vector_size){
        std::cerr << "ERROR! Input to FVMSteadyStateSolver::convertVectorToStateSpace() has wrong dimensions." << std::endl;
        std::cerr << v.rows() << " != " << vector_size << "! Exiting." << std::endl;
        exit(0);
    }

    //fill density values
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        rho(e) = v(e);
    }

    //fill momentum values
    int offset = driver->fluid_grid->element_count;
    int GRID_DIM = driver->fluid_grid->GRID_DIM;
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        for (int i=0; i<GRID_DIM; i++){
            p(e,i) = v(offset + e*GRID_DIM + i);
        }
    }

    //if thermal simulation, fill energy values
    if (driver->TYPE == FiniteVolumeDriver::THERMAL){
        offset = driver->fluid_grid->element_count * (1 + GRID_DIM);
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            rhoE(e) = v(offset + e);
        }
    }

    return;
}


void FVMSteadyStateSolver::convertStateSpaceToVector(Job* job, FiniteVolumeDriver* driver,
                                                     Eigen::VectorXd& v,
                               const Eigen::VectorXd& rho,
                               const KinematicVectorArray& p,
                               const Eigen::VectorXd& rhoE){

    //check vector for correct dimensions
    if (v.rows() != vector_size){
        v = Eigen::VectorXd(vector_size);
    }

    //fill density values
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        v(e) = rho(e);
    }

    //fill momentum values
    int offset = driver->fluid_grid->element_count;
    int GRID_DIM = driver->fluid_grid->GRID_DIM;
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        for (int i=0; i<GRID_DIM; i++){
            v(offset + e*GRID_DIM + i) = p[e][i];
        }
    }

    //if thermal simulation, fill energy values
    if (driver->TYPE == FiniteVolumeDriver::THERMAL){
        offset = driver->fluid_grid->element_count * (1 + GRID_DIM);
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            v(offset + e) = rhoE(e);
        }
    }

    return;
}



void FVMSteadyStateSolver::step(Job* job, FiniteVolumeDriver* driver){
    //initial residual
    double step_norm = 1.0;
    double r0_norm = 1.0;
    double r_norm = 1.0;

    //initialize guess vector and BiCG-STAB vectors
    Eigen::VectorXd x_i = Eigen::VectorXd(vector_size);
    Eigen::VectorXd v_i = Eigen::VectorXd(vector_size);
    Eigen::VectorXd p_i = Eigen::VectorXd(vector_size);
    Eigen::VectorXd r_i = Eigen::VectorXd(vector_size);
    Eigen::VectorXd r_0 = Eigen::VectorXd(vector_size);
    Eigen::VectorXd s = Eigen::VectorXd(vector_size);
    Eigen::VectorXd t = Eigen::VectorXd(vector_size);

    //jacobi iteration vectors
    Eigen::VectorXd k_1 = Eigen::VectorXd(vector_size);
    Eigen::VectorXd k_2 = Eigen::VectorXd(vector_size);
    Eigen::VectorXd k_3 = Eigen::VectorXd(vector_size);
    Eigen::VectorXd k_4 = Eigen::VectorXd(vector_size);

    //continuation step
    int step_count = 0;
    do {
        //initialize state vector
        convertStateSpaceToVector(job, driver,
                                  u_0,
                                  driver->fluid_body->rho,
                                  driver->fluid_body->p,
                                  driver->fluid_body->rhoE);

        //if first iteration, save value
        //norm calculated without continuation
        r_norm = F(job, driver, u_0, 0.0).norm();
        if (step_count < 1) {
            step_norm = r_norm;
        }

        //check for convergence
        std::cout << "Continuation Iteration: " << step_count << std::endl;
        std::cout << "    nu:                 " << nu << std::endl;
        std::cout << "    Error:              " << r_norm << std::endl;
        std::cout << "    Relative Error:     " << r_norm / step_norm << std::endl;
        if (r_norm < ABS_TOL || r_norm < REL_TOL * step_norm) {
            std::cout << "Converged!" << std::endl << std::endl;
            break;
        }

        //newton iterator
        int n = 0;

        //do while loop for Newton-Raphson solver
        do {
            //initialize state vector and residual
            convertStateSpaceToVector(job, driver,
                                      u_n,
                                      driver->fluid_body->rho,
                                      driver->fluid_body->p,
                                      driver->fluid_body->rhoE);

            r_n = F(job, driver, u_n, nu);

            //if first iteration, save value
            r_norm = r_n.norm();
            if (n < 1) {
                r0_norm = r_norm;
            }

            //check for convergence
            std::cout << "    Newton-Raphson Iteration: " << n << std::endl;
            std::cout << "        Error:                " << r_norm << std::endl;
            std::cout << "        Relative Error:       " << r_norm / r0_norm << std::endl;
            if (r_norm < ABS_TOL || r_norm < REL_TOL * r0_norm) {
                std::cout << "    Converged!" << std::endl << std::endl;
                break;
            }

            //if you get here, we haven't converged yet
            if (INNER_SOLVER == BICGSTAB) {
                //bicgstab solver
                int i = 0;
                double rho_iminus = 1;
                double alpha = 1;
                double w_iminus = 1;
                double rho_i, w_i, beta;
                p_i.setZero();
                v_i.setZero();
                x_i.setZero();
                r_0 = r_n;
                r_i = r_n;

                //do while loop to solve BiCG-STAB problem
                do {
                    rho_i = r_0.dot(r_i);

                    //check that r_0 and r_i not orthogonal
                    if (std::sqrt(rho_i) < ABS_TOL || std::sqrt(rho_i) < REL_TOL * r_0.norm()) {
                        //if so, resample r_i and re-initialize
                        r_0 = r_i;
                        rho_i = r_0.dot(r_i);
                        rho_iminus = 1;
                        alpha = 1;
                        w_iminus = 1;
                        std::cout << "    Bump!" << std::endl;
                    }

                    beta = (rho_i / rho_iminus) * (alpha / w_iminus);
                    p_i = r_i + beta * (p_i - w_iminus * v_i);
                    v_i = DF(job, driver, p_i, h) * p_i.norm();
                    alpha = rho_i / (r_0.dot(v_i));

                    //check that r_0 and v_i not orthogonal
                    if (!std::isfinite(alpha)) {
                        //if so, resample r_i and re-initialize
                        r_0 = r_i;
                        rho_i = r_0.dot(r_i);
                        alpha = 1;
                        w_iminus = 1;
                        w_i = 1;
                        std::cout << "    Bump!" << std::endl;
                        s = r_i;
                    } else {
                        s = r_i - alpha * v_i;
                    }

                    //check that s is non-zero
                    if (s.norm() > ABS_TOL) {
                        t = DF(job, driver, s, h) * s.norm();
                        w_i = t.dot(s) / t.squaredNorm();
                        x_i = x_i + alpha * p_i + w_i * s;
                        r_i = s - w_i * t;
                    } else {
                        x_i = x_i + alpha * p_i;
                        r_i = s;
                    }

                    //check for convergence
                    std::cout << "        BiCGSTAB Iteration: " << i << "    " << std::endl;
                    std::cout << "            Error:          " << r_i.norm() << "    " << std::endl;
                    std::cout << "            Relative Error: " << r_i.norm() / r_norm << "    " << std::endl;
                    std::cout << "\e[A\e[A\e[A" << std::flush;
                    /*
                    std::cout << "x " << x_i.norm() << std::endl;
                    std::cout << "s " << s.norm() << " t " << t.norm() << std::endl;
                    std::cout << "p " << p_i.norm() << std::endl;
                    std::cout << "alpha " << alpha << std::endl;
                    std::cout << "omega " << w_i << std::endl;
                    std::cout << "beta " << beta << std::endl;
                     */

                    if (r_i.norm() < ABS_TOL || r_i.norm() < REL_TOL * r_norm) {
                        //converged
                        break;
                    } else if (r_i.norm() > BiCGSTAB_TOL * r_norm){
                        //probably not going to converge
                        break;
                    }

                    //advance residual
                    rho_iminus = rho_i;
                    w_iminus = w_i;

                    //iterate on i
                    i++;

                } while (i < max_i);

                std::cout << std::endl << std::endl << std::endl;
                //if i exceeds maximum, print statement
                if (i >= max_i) {
                    std::cout << "        Exceeded max iterations! Exiting BiCGSTAB loop." << std::endl;
                } else if (r_i.norm() > BiCGSTAB_TOL * r_norm){
                    std::cout << "        Residual not shrinking! Exiting BiCGSTAB loop." << std::endl;
                }

            } else if (INNER_SOLVER == DIRECT) {
                //direct solver
                //form sparse matrix
                Eigen::SparseMatrix<double> dFdU(vector_size, vector_size);

                //fill entries of dFdU using numerical derivative
                v_i.setZero();
                for (int ii = 0; ii < vector_size; ii++) {
                    //initialize ii component of vector as 1
                    v_i(ii) = 1.0;

                    //calculate resulting gradient
                    p_i = DF(job, driver, v_i, h);

                    //loop over elements of p_i and add to sparse matrix
                    for (int jj = 0; jj < vector_size; jj++) {
                        //only add non-zero entries
                        if (std::abs(p_i(jj)) > ABS_TOL) {
                            dFdU.insert(jj, ii) = p_i(jj);
                        }
                    }

                    //re-zero ii component of vector
                    v_i(ii) = 0.0;
                    std::cout << "    " << ii << " / " << vector_size << "\r" << std::flush;
                }
                std::cout << std::endl << "        dFdU filled!\r" << std::endl;

                //use Eigen to solve system of equations
                dFdU.makeCompressed();
                Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
                solver.compute(dFdU);

                if (solver.info() != Eigen::Success) {
                    // decomposition failed
                    std::cerr << "ERROR! Eigen::SparseLU failed to decompose dFdU. Exiting." << std::endl;
                    exit(0);
                }

                x_i = solver.solve(r_n);

                if (solver.info() != Eigen::Success) {
                    // solving failed
                    std::cerr << "ERROR! Eigen::SparseLU failed to solver dFdU*x = r. Exiting." << std::endl;
                    exit(0);
                }
            } else if (INNER_SOLVER == JACOBI){
                //essentially explicit time marching
                k_1 = F(job, driver, u_n, 0.0)/nu;
                k_2 = F(job, driver, (u_n + k_1/2.0), 0.0)/nu;
                k_3 = F(job, driver, (u_n + k_2/2.0), 0.0)/nu;
                k_4 = F(job, driver, (u_n + k_3), 0.0)/nu;

                //adjust u_n
                u_n = u_n + (k_1 + 2.0*k_2 + 2.0*k_3 + k_4)/6.0;

                //need to break out of newton loop, too
                break;
            }

            /*
            std::cout << "u " << u_n.norm() << std::endl;
            std::cout << "x " << x_i.norm() << std::endl;
             */

            //search along x for best fit
            //u_n = u_n - x_i;
            r_0 = r_n;
            r_n = F(job, driver, (u_n - x_i), nu);
            while (r_n.norm() > r_0.norm() || !std::isfinite(r_n.norm())){
                //drawback
                x_i = 0.5 * x_i;

                //check residual
                r_n = F(job, driver, (u_n - x_i), nu);

                std::cout << "        Line Search Residual: " << r_n.norm() << "\r";

            }
            std::cout << std::endl;

            u_n = u_n - x_i;

            /*
            //first, ensure initial guess is finite
            do {
                //check residual
                r_0 = F(job, driver, (u_n - x_i), nu);
                x_i = 0.5*x_i;

            } while (!std::isfinite(r_0.norm()));
            x_i = 2.0*x_i;

            //minimize along x_i
            double rho_0, rho, lambda;
            double drdl_max, drdl_min, drdl;
            double lambda_max = 1.0;
            double lambda_min = 0.0;

            drdl_max = (F(job, driver, (u_n - h*x_i),nu).norm() - F(job, driver, (u_n + h*x_i),nu).norm())/(2.0*h); //gradient at lambda = 0
            drdl_min = (F(job, driver, (u_n - (1.0 + h)*x_i),nu).norm() - F(job, driver, (u_n - (1.0 - h)*x_i),nu).norm())/(2.0*h); //gradient at lambda = 1

            int ii=0;
            do {
                lambda = 0.5*(lambda_max + lambda_min);
                drdl = (F(job, driver, (u_n - (lambda + h)*x_i),nu).norm() - F(job, driver, (u_n - (lambda - h)*x_i),nu).norm())/(2.0*h);

                //check sign
                if (drdl*drdl_min > 0){
                    //same sign
                    lambda_max = lambda;
                } else {
                    //dif sign
                    lambda_min = lambda;
                }

                std::cout << "    Line Search Iteration: " << ii << std::endl;
                std::cout << "        Error:             " << F(job, driver, (u_n - (lambda)*x_i),nu).norm() << std::endl;
                std::cout << "        Gradient:          " << drdl << std::endl;

                ii++;
            } while (drdl > ABS_TOL && ii < max_i);

            //lambda determined
            x_i = lambda*x_i;

            u_n = u_n - x_i;
            */


            convertVectorToStateSpace(job, driver,
                                      u_n,
                                      driver->fluid_body->rho,
                                      driver->fluid_body->p,
                                      driver->fluid_body->rhoE);

            //iterate on n
            n++;

            if (n >= max_n) {
                std::cout << "    Exceeded max iterations! Exiting Newton-Raphson solver." << std::endl;
            } else if (x_i.norm() < ABS_TOL || (r_0.norm() - r_n.norm())/(r_0.norm()) < LineSearch_TOL){
                //singular jacobian
                std::cout << "    Singular Jacobian! Exiting Newton-Raphson solver." << std::endl;
                n = max_n;
                break;
            }

        } while (n < max_n);

        if (n < max_n && INNER_SOLVER != JACOBI) {
            //converged
            //be more aggresive
            nu /= 10.0;

            //proceed to next step
            u_0 = u_n;

            //iterate step count
            step_count++;

        } else if (INNER_SOLVER == JACOBI){
            //do not change nu
            u_0 = u_n;

            //iterate step count
            step_count++;

        } else {
            //be more cautious
            nu *= 2.0;

            //try again
            //u_0 unchanged
        }

        /*
        else if (F(job, driver, u_n, 0.0).norm() < F(job, driver, u_0, 0.0).norm()){
            //did not converge, but still an ok guess
            nu *= 2.0;
            u_0 = u_n;

            //iterate step count
            step_count++;
        }*/

        convertVectorToStateSpace(job, driver,
                                  u_0,
                                  driver->fluid_body->rho,
                                  driver->fluid_body->p,
                                  driver->fluid_body->rhoE);

        if (step_count >= max_s) {
            std::cout << "Exceeded max iterations! Continuing Psuedo-Transient Continuation solver." << std::endl;

            //write frame to file
            driver->serializer->writeFrame(job,driver);

            //reset counter
            step_count = 1;
        }

    } while (step_count < max_s);

    //end step
    return;
}
//
// Created by aaron on 12/24/19.
// fvm_runge_kutta_solver.cpp
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
#include "fvm_solvers.hpp"

void FVMRungeKuttaSolver::init(Job* job, FiniteVolumeDriver* driver){
    if (int_props.size() < 1){
        //do nothing
    } else {
        order = int_props[0];
        if (order != 4){
            std::cout << "ERROR: FVMRungeKuttaSolver only has RK4 implemented." << std::endl;
            order = 4;
        }
    }
    std::cout << "FiniteVolumeSolver properties: (ORDER = " << order << ")" << std::endl;

    //size flux containers
    density_fluxes = Eigen::VectorXd(driver->fluid_grid->element_count);
    momentum_fluxes = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    energy_fluxes = Eigen::VectorXd(driver->fluid_grid->element_count);

    //size containers for RK method
    rho_k1 = density_fluxes;
    rho_k2 = density_fluxes;
    rho_k3 = density_fluxes;
    rho_k4 = density_fluxes;

    p_k1 = momentum_fluxes;
    p_k2 = momentum_fluxes;
    p_k3 = momentum_fluxes;
    p_k4 = momentum_fluxes;

    rhoE_k1 = energy_fluxes;
    rhoE_k2 = energy_fluxes;
    rhoE_k3 = energy_fluxes;
    rhoE_k4 = energy_fluxes;

    //done
    std::cout << "FiniteVolumeSolver initialized." << std::endl;
    return;
}

void FVMRungeKuttaSolver::step(Job* job, FiniteVolumeDriver* driver){

    if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
        /* STEP ONE */
        //call grid to reconstruct conserved fields
        driver->fluid_grid->constructDensityField(job, driver);
        driver->fluid_grid->constructMomentumField(job, driver);
        driver->fluid_grid->constructPorosityField(job, driver);

        //k1 = f(y_0)
        //fluid fluxes associated with each volume
        rho_k1 = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
        p_k1 = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2

        //update state to calculate k2 = f(y_0 + k1*dt/2)
        double volume;
        for (int e = 0; e < driver->fluid_grid->element_count; e++) {
            volume = driver->fluid_grid->getElementVolume(e);

            //update p_k1
            p_k1[e] += driver->fluid_body->rho(e) * driver->gravity * volume;

            driver->fluid_body->rho(e) += rho_k1(e) / volume * job->dt / 2.0;     //d(rho dv)/dt   = flux
            driver->fluid_body->p[e] += p_k1[e] / volume * job->dt / 2.0;    //d(rho u dv)/dt = flux
        }

        /*STEP TWO*/
        //k2 = f(y_0 + k1*dt/2)
        driver->fluid_grid->constructDensityField(job, driver);
        driver->fluid_grid->constructMomentumField(job, driver);
        driver->fluid_grid->constructPorosityField(job, driver);
        rho_k2 = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
        p_k2 = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2

        //update state to calculate k3 = f(y_0 + k2*dt/2)
        for (int e = 0; e < driver->fluid_grid->element_count; e++) {
            volume = driver->fluid_grid->getElementVolume(e);

            //update p_k2
            p_k2[e] += driver->fluid_body->rho(e) * driver->gravity * volume;

            driver->fluid_body->rho(e) += (rho_k2(e) - rho_k1(e)) / volume * job->dt / 2.0;     //d(rho dv)/dt   = flux
            driver->fluid_body->p[e] += (p_k2[e] - p_k1[e]) / volume * job->dt / 2.0;    //d(rho u dv)/dt = flux
        }

        /*STEP THREE*/
        //k3 = f(y_0 + k2*dt/2)
        driver->fluid_grid->constructDensityField(job, driver);
        driver->fluid_grid->constructMomentumField(job, driver);
        driver->fluid_grid->constructPorosityField(job, driver);
        rho_k3 = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
        p_k3 = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2

        //update state to calculate k4 = f(y_0 + k3*dt)
        for (int e = 0; e < driver->fluid_grid->element_count; e++) {
            volume = driver->fluid_grid->getElementVolume(e);

            //update p_k2
            p_k3[e] += driver->fluid_body->rho(e) * driver->gravity * volume;

            driver->fluid_body->rho(e) +=
                    (2.0 * rho_k3(e) - rho_k2(e)) / volume * job->dt / 2.0;     //d(rho dv)/dt   = flux
            driver->fluid_body->p[e] += (2.0 * p_k3[e] - p_k2[e]) / volume * job->dt / 2.0;    //d(rho u dv)/dt = flux
        }

        /*STEP FOUR*/
        //k4 = f(y_0 + k3*dt)
        driver->fluid_grid->constructDensityField(job, driver);
        driver->fluid_grid->constructMomentumField(job, driver);
        driver->fluid_grid->constructPorosityField(job, driver);
        rho_k4 = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
        p_k4 = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2

        //current state is y_0 + k_3*dt
        for (int e = 0; e < driver->fluid_grid->element_count; e++) {
            volume = driver->fluid_grid->getElementVolume(e);

            //update p_k2
            p_k4[e] += driver->fluid_body->rho(e) * driver->gravity * volume;

            driver->fluid_body->rho(e) +=
                    (rho_k1(e) / 6.0 + rho_k2(e) / 3.0 - 2.0 * rho_k3(e) / 3.0 + rho_k4(e) / 6.0) / volume *
                    job->dt;     //d(rho dv)/dt   = flux
            driver->fluid_body->p[e] += (p_k1[e] / 6.0 + p_k2[e] / 3.0 - 2.0 * p_k3[e] / 3.0 + p_k4[e] / 6.0) / volume *
                                        job->dt;    //d(rho u dv)/dt = flux
        }
    } else if (driver->TYPE == FiniteVolumeDriver::THERMAL) {
        /* STEP ONE */
        //call grid to reconstruct conserved fields
        driver->fluid_grid->constructDensityField(job, driver);
        driver->fluid_grid->constructMomentumField(job, driver);
        driver->fluid_grid->constructEnergyField(job, driver);
        driver->fluid_grid->constructPorosityField(job, driver);

        //k1 = f(y_0)
        //fluid fluxes associated with each volume
        rho_k1 = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
        p_k1 = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2
        rhoE_k1 = driver->fluid_grid->calculateElementEnergyFluxes(job, driver);

        //update state to calculate k2 = f(y_0 + k1*dt/2)
        double volume;
        for (int e = 0; e < driver->fluid_grid->element_count; e++) {
            volume = driver->fluid_grid->getElementVolume(e);

            //update p_k1
            p_k1[e] += driver->fluid_body->rho(e) * driver->gravity * volume;
            rhoE_k1(e) += driver->gravity.dot(driver->fluid_body->p[e]) * volume;

            driver->fluid_body->rho(e) += rho_k1(e) / volume * job->dt / 2.0;     //d(rho dv)/dt   = flux
            driver->fluid_body->p[e] += p_k1[e] / volume * job->dt / 2.0;    //d(rho u dv)/dt = flux
            driver->fluid_body->rhoE(e) += rhoE_k1(e) / volume *job->dt / 2.0;
        }

        /*STEP TWO*/
        //k2 = f(y_0 + k1*dt/2)
        driver->fluid_grid->constructDensityField(job, driver);
        driver->fluid_grid->constructMomentumField(job, driver);
        driver->fluid_grid->constructEnergyField(job, driver);
        driver->fluid_grid->constructPorosityField(job, driver);
        rho_k2 = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
        p_k2 = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2
        rhoE_k2 = driver->fluid_grid->calculateElementEnergyFluxes(job, driver);

        //update state to calculate k3 = f(y_0 + k2*dt/2)
        for (int e = 0; e < driver->fluid_grid->element_count; e++) {
            volume = driver->fluid_grid->getElementVolume(e);

            //update p_k2
            p_k2[e] += driver->fluid_body->rho(e) * driver->gravity * volume;
            rhoE_k2(e) += driver->gravity.dot(driver->fluid_body->p[e]) * volume;

            driver->fluid_body->rho(e) += (rho_k2(e) - rho_k1(e)) / volume * job->dt / 2.0;     //d(rho dv)/dt   = flux
            driver->fluid_body->p[e] += (p_k2[e] - p_k1[e]) / volume * job->dt / 2.0;    //d(rho u dv)/dt = flux
            driver->fluid_body->rhoE(e) += (rhoE_k2(e) - rhoE_k1(e)) / volume * job->dt / 2.0;
        }

        /*STEP THREE*/
        //k3 = f(y_0 + k2*dt/2)
        driver->fluid_grid->constructDensityField(job, driver);
        driver->fluid_grid->constructMomentumField(job, driver);
        driver->fluid_grid->constructEnergyField(job, driver);
        driver->fluid_grid->constructPorosityField(job, driver);
        rho_k3 = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
        p_k3 = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2
        rhoE_k3 = driver->fluid_grid->calculateElementEnergyFluxes(job, driver);

        //update state to calculate k4 = f(y_0 + k3*dt)
        for (int e = 0; e < driver->fluid_grid->element_count; e++) {
            volume = driver->fluid_grid->getElementVolume(e);

            //update p_k2
            p_k3[e] += driver->fluid_body->rho(e) * driver->gravity * volume;
            rhoE_k3(e) += driver->gravity.dot(driver->fluid_body->p[e]) * volume;

            driver->fluid_body->rho(e) +=
                    (2.0 * rho_k3(e) - rho_k2(e)) / volume * job->dt / 2.0;     //d(rho dv)/dt   = flux
            driver->fluid_body->p[e] += (2.0 * p_k3[e] - p_k2[e]) / volume * job->dt / 2.0;    //d(rho u dv)/dt = flux
            driver->fluid_body->rhoE(e) += (2.0 * rhoE_k3(e) - rhoE_k2(e)) / volume * job->dt / 2.0;
        }

        /*STEP FOUR*/
        //k4 = f(y_0 + k3*dt)
        driver->fluid_grid->constructDensityField(job, driver);
        driver->fluid_grid->constructMomentumField(job, driver);
        driver->fluid_grid->constructEnergyField(job, driver);
        driver->fluid_grid->constructPorosityField(job, driver);
        rho_k4 = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
        p_k4 = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2
        rhoE_k4 = driver->fluid_grid->calculateElementEnergyFluxes(job, driver);

        //current state is y_0 + k_3*dt
        for (int e = 0; e < driver->fluid_grid->element_count; e++) {
            volume = driver->fluid_grid->getElementVolume(e);

            //update p_k2
            p_k4[e] += driver->fluid_body->rho(e) * driver->gravity * volume;
            rhoE_k4(e) += driver->gravity.dot(driver->fluid_body->p[e]) * volume;

            driver->fluid_body->rho(e) +=
                    (rho_k1(e) / 6.0 + rho_k2(e) / 3.0 - 2.0 * rho_k3(e) / 3.0 + rho_k4(e) / 6.0) / volume *
                    job->dt;     //d(rho dv)/dt   = flux
            driver->fluid_body->p[e] += (p_k1[e] / 6.0 + p_k2[e] / 3.0 - 2.0 * p_k3[e] / 3.0 + p_k4[e] / 6.0) / volume *
                                        job->dt;    //d(rho u dv)/dt = flux

            driver->fluid_body->rhoE(e) += (rhoE_k1(e) / 6.0 + rhoE_k2(e) / 3.0 - 2.0 * rhoE_k3(e) / 3.0 + rhoE_k4(e) / 6.0) / volume *
                                            job->dt;
        }
    } else {
        //shouldn't have another flag...
        std::cerr << "ERROR! FVMDefaultSolver does not have simulation TYPE " << driver->TYPE << " implemented! Exiting." << std::endl;
        exit(0);
    }

    //end step
    return;
}
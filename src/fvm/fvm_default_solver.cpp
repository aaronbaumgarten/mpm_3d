//
// Created by aaron on 12/23/19.
// fvm_default_solver.cpp
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

void FVMDefaultSolver::init(Job* job, FiniteVolumeDriver* driver){
    //size flux containers
    density_fluxes = Eigen::VectorXd(driver->fluid_grid->element_count);
    momentum_fluxes = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    //done
    std::cout << "FiniteVolumeSolver initialized." << std::endl;
    return;
}

void FVMDefaultSolver::step(Job* job, FiniteVolumeDriver* driver){
    //call grid to reconstruct conserved fields
    driver->fluid_grid->constructDensityField(job, driver);
    driver->fluid_grid->constructMomentumField(job, driver);

    //fluid fluxes associated with each volume
    density_fluxes = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
    momentum_fluxes = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2

    //use fluxes to update element-wise average density and momentum (Forward Euler)
    double volume;
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        volume = driver->fluid_grid->getElementVolume(e);
        driver->fluid_body->rho(e) += density_fluxes(e) / volume * job->dt;     //d(rho dv)/dt   = flux
        driver->fluid_body->p[e]   += momentum_fluxes[e] / volume * job->dt;    //d(rho u dv)/dt = flux

        //add gravitational contribution to momentum update
        driver->fluid_body->p[e]   += driver->fluid_body->rho(e)*driver->gravity * job->dt;
    }

    //end step
    return;
}
//
// Created by aaron on 12/23/19.
// fvm_default_body.cpp
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
#include "fvm_bodies.hpp"


void FVMDefaultBody::init(Job *job, FiniteVolumeDriver *driver) {
    double rho_0, theta_0;
    if (fp64_props.size() < 2) {
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined (rho_0, theta_0).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        rho_0 = fp64_props[0];
        theta_0 = fp64_props[1];
        std::cout << "FiniteVolumeBody properties (rho_0 = " << rho_0 << ", theta_0 = " << theta_0 << ")." << std::endl;
    }

    //size vectors using grid (which should be initialized first)
    //momentum gradient
    p_x = KinematicTensorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    p_x.setZero();

    //momentum
    p = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    p.setZero();

    //density gradient
    rho_x = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    rho_x.setZero();

    //shear stress
    tau = MaterialTensorArray(driver->fluid_grid->element_count);
    tau.setZero();

    //density
    rho = Eigen::VectorXd(driver->fluid_grid->element_count);
    rho.setConstant(rho_0);

    //pressure
    P = Eigen::VectorXd(driver->fluid_grid->element_count);
    P.setZero();

    //temperature
    theta = Eigen::VectorXd(driver->fluid_grid->element_count);
    theta.setConstant(theta_0);

    std::cout << "FiniteVolumeBody initialized." << std::endl;
}
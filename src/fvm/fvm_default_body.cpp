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

    if (int_props.size() > 0) {
        if (int_props[0] == 1) {
            HYDROSTATIC_INITIALIZATION = true;
            std::cout << "Initializing FiniteVolumeBody as hydrostatic." << std::endl;
        }
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

    //porosity fields
    n = Eigen::VectorXd(job->grid->node_count);
    n.setConstant(1.0);

    n_e = Eigen::VectorXd(driver->fluid_grid->element_count);
    n_e.setConstant(1.0);
    true_density_x = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    true_density_x.setZero();

    //adjust density and porosity fields
    if (driver->fluid_material->calculatePorosity(job, driver) == 1){
        //if 1 is returned, then this is a mixture problem

        //integrate porosity field over element volumes
        n_e = driver->fluid_grid->M * driver->fluid_body->n;

        //adjust density accordingly
        for (int e=0; e<driver->fluid_grid->element_count; e++) {
            //average porosity over element volume
            n_e(e) /= driver->fluid_grid->getElementVolume(e);

            //\bar{\rho} = n*rho
            rho(e) *= n_e(e);
        }


        /*
        //use centroid value of porosity
        double tmp_n = 1;
        std::vector<int> nvec(0);
        std::vector<double> valvec(0);
        for (int e=0; e<driver->fluid_grid->element_count; e++) {
            nvec.resize(0);
            valvec.resize(0);

            KinematicVector tmp_x = driver->fluid_grid->getElementCentroid(job, e);
            job->grid->evaluateBasisFnValue(job,tmp_x,nvec,valvec);

            tmp_n = 0;
            for (int i=0; i<nvec.size(); i++){
                tmp_n += n(nvec[i])*valvec[i];
            }

            //\bar{\rho} = n*rho
            rho(e) *= tmp_n;
        }
         */
    }

    //temperature
    theta = Eigen::VectorXd(driver->fluid_grid->element_count);
    theta.setConstant(theta_0);

    //pressure
    P = Eigen::VectorXd(driver->fluid_grid->element_count);
    double P_0 = driver->fluid_material->getPressureFromDensityAndTemperature(job, driver, rho_0, theta_0, 1.0);
    P.setConstant(P_0);

    //create hydrostatic pressure distribution
    if (HYDROSTATIC_INITIALIZATION){
        //find approximate centroid of domain
        double domain_volume = 0;
        KinematicVector x_centroid = KinematicVector(job->JOB_TYPE);
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            domain_volume += driver->fluid_grid->getElementVolume(e);
            x_centroid += driver->fluid_grid->getElementCentroid(job, e) * driver->fluid_grid->getElementVolume(e);
        }
        x_centroid /= domain_volume;

        //define pressure gradient
        KinematicVector P_grad = rho_0*driver->gravity;

        //assign pressure based on location relative to centroid
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            P(e) = P_0 + P_grad.dot(driver->fluid_grid->getElementCentroid(job, e) - x_centroid);
        }

        //adjust densities
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            rho(e) = driver->fluid_material->getDensityFromPressureAndTemperature(job, driver, P(e), theta_0, n_e(e));
        }
    }

    //energy
    rhoE = Eigen::VectorXd(driver->fluid_grid->element_count);
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        rhoE(e) = driver->fluid_material->getInternalEnergyFromPressureAndTemperature(job, driver, P(e), theta_0, n_e(e));
    }

    //energy gradient
    rhoE_x = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    rhoE_x.setZero();

    //solid phase velocity
    v_s = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    v_s.setZero();

    std::cout << "FiniteVolumeBody initialized." << std::endl;
}
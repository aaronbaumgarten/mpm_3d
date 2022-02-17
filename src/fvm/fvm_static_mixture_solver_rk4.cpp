//
// Created by aaron on 5/18/20.
// fvm_static_mixture_solver_rk4.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <job.hpp>
#include <time.h>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "mpm_objects.hpp"
#include "fvm_objects.hpp"
#include "fvm_solvers.hpp"
#include "objects/bodies/bodies.hpp"

/*----------------------------------------------------------------------------*/
//one forward mixed time-stpe

void FVMStaticMixtureSolverRK4::step(Job* job, FiniteVolumeDriver* driver){

    /*----------------------*/
    /*  Begin FVM-MPM Step  */
    /*----------------------*/

    //no need to do anything to MPM (static)

    //map mixture properties to finite volumes (only on first step)
    if (job->t == 0) {
        mapMixturePropertiesToElements(job, driver);
        std::cout << "SUCCESS: Mapped mixture properties to grid quadrature points." << std::endl;
    }

    /*----------------------*/
    /*  Begin FVM-RK4 Step  */
    /*----------------------*/

    convertStateSpaceToVector(job,
                              driver,
                              u_0,
                              driver->fluid_body->rho,
                              driver->fluid_body->p,
                              driver->fluid_body->rhoE);

    //get drag and corrected coefficients from initial state
    //driver->fluid_grid->calculateSplitIntegralDragForces(job, driver, f_d, f_d_e);
    K_n = driver->fluid_grid->getCorrectedDragCoefficients(job, driver);

    /*----------------------*/
    /* Estimate k1 w/o drag */
    /*----------------------*/

    f_d.setZero();
    f_d_e.setZero();
    k1 = job->dt * F(job, driver, u_0);
    f_b_e = f_e/6.0;                        //store intermediate buoyancy term

    /*-----------------------*/
    /* Estimate Drag at t+dt */
    /*-----------------------*/

    //assign u_0 + k1 to fluid phase
    convertVectorToStateSpace(job,
                              driver,
                              (u_0 + k1),
                              driver->fluid_body->rho,
                              driver->fluid_body->p,
                              driver->fluid_body->rhoE);

    driver->fluid_grid->mapMixturePropertiesToQuadraturePoints(job, driver);
    driver->fluid_grid->constructDensityField(job, driver);
    driver->fluid_grid->constructMomentumField(job, driver);
    if (driver->TYPE == driver->THERMAL) {
        driver->fluid_grid->constructEnergyField(job, driver);
    }
    driver->fluid_grid->constructPorosityField(job, driver);
    driver->fluid_grid->constructVelocityField(job, driver);


    //get corrected estimate for drag at t+dt
    driver->fluid_grid->calculateSplitIntegralCorrectedDragForces(job, driver, f_d, f_d_e, K_n);

    //add f_d to k1
    convertVectorToStateSpace(job,
                              driver,
                              k1,
                              density_fluxes,
                              momentum_fluxes,
                              energy_fluxes);

    momentum_fluxes -= job->dt*f_d_e;

    convertStateSpaceToVector(job,
                              driver,
                              k1,
                              density_fluxes,
                              momentum_fluxes,
                              energy_fluxes);

    /*----------------------*/
    /*    Find k2 to k4     */
    /*----------------------*/

    k2 = job->dt * F(job, driver, (u_0 + k1/2.0));
    f_b_e += (f_e - f_d_e)/3.0;

    k3 = job->dt * F(job, driver, (u_0 + k2/2.0));
    f_b_e += (f_e - f_d_e)/3.0;

    k4 = job->dt * F(job, driver, (u_0 + k3));
    f_b_e += (f_e - f_d_e)/6.0;

    u_n = u_0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;

    convertVectorToStateSpace(job,
                              driver,
                              u_n,
                              driver->fluid_body->rho,
                              driver->fluid_body->p,
                              driver->fluid_body->rhoE);

    /*-----------------------*/
    /* Estimate Drag at t+dt */
    /*-----------------------*/

    //subtract predicted drag force from fluid to get u_plus
    for (int e=0; e<driver->fluid_grid->element_count; e++) {
        driver->fluid_body->p[e] += job->dt * f_d_e[e];
    }

    //reconstruct fluid fields
    driver->fluid_grid->constructDensityField(job, driver);
    driver->fluid_grid->constructMomentumField(job, driver);
    if (driver->TYPE == driver->THERMAL) {
        driver->fluid_grid->constructEnergyField(job, driver);
    }
    driver->fluid_grid->constructPorosityField(job, driver);
    driver->fluid_grid->constructVelocityField(job, driver);

    //get correction term
    driver->fluid_grid->calculateSplitIntegralCorrectedDragForces(job, driver, second_drag_correction, f_d_e, K_n);

    //add corrected drag to fluid momentum
    double volume;
    for (int e = 0; e < driver->fluid_grid->element_count; e++) {
        driver->fluid_body->p[e] -= job->dt * f_d_e[e];
    }

    //get energy contribution at end of step
    if (driver->TYPE == driver->THERMAL) {
        energy_fluxes = driver->fluid_grid->calculateInterphaseEnergyFlux(job, driver);
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            driver->fluid_body->rhoE(e) += job->dt * energy_fluxes(e) / driver->fluid_grid->getElementVolume(e);
        }
    }

    /*----------------------*/
    /*   End FVM-RK4 Step   */
    /*----------------------*/

    //do nothing

    /*----------------------*/
    /*     End MPM Step     */
    /*----------------------*/


    return;
}
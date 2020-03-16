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

void FVMMixtureSolverRK4::init(Job* job, FiniteVolumeDriver* driver){

    //call parent initialization f'n
    FVMMixtureSolver::init(job, driver);

    //use Runge-Kutta 4 to find dU/dt
    //U is concatenated element state vector (rho, rho*u, rho*E)

    //size of vectors
    if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL){
        //rho, rho*u, mx_t
        vector_size = driver->fluid_grid->element_count*(1 + driver->fluid_grid->GRID_DIM);
    } else if (driver->TYPE == FiniteVolumeDriver::THERMAL){
        //rho, rho*u, rho*E, mx_t
        vector_size = driver->fluid_grid->element_count*(2 + driver->fluid_grid->GRID_DIM);
    } else {
        std::cerr << "ERROR! FVMSteadyStateSolver not implemented for simulationg TYPE " << driver->TYPE << "! Exiting." << std::endl;
        exit(0);
    }

    //done
    std::cout << "FiniteVolumeSolver initialized." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//one forward mixed time-stpe

void FVMMixtureSolverRK4::step(Job* job, FiniteVolumeDriver* driver){

    /*----------------------*/
    /*  Begin FVM-MPM Step  */
    /*----------------------*/

    //create map
    createMappings(job);

    //map particles to grid
    mapPointsToNodes(job);

    //add arbitrary loading conditions
    generateLoads(job);
    applyLoads(job);

    //add contact forces
    generateContacts(job);
    addContacts(job);

    //enforce boundary conditions
    generateBoundaryConditions(job);
    addBoundaryConditions(job);

    //map mixture properties to finite volumes
    mapMixturePropertiesToElements(job, driver);

    /*----------------------*/
    /*  Begin FVM-RK4 Step  */
    /*----------------------*/

    convertStateSpaceToVector(job,
                              driver,
                              u_0,
                              driver->fluid_body->rho,
                              driver->fluid_body->p,
                              driver->fluid_body->rhoE);

    tmp_mx_t = job->bodies[solid_body_id]->nodes->mx_t;
    k1 = job->dt * F(job, driver, u_0, tmp_mx_t);
    f_i1 = f_i;                             //store intermediate drag calculations

    for (int i=0; i<tmp_mx_t.size(); i++){
        tmp_mx_t(i) = job->bodies[solid_body_id]->nodes->mx_t[i]
                      + job->dt*(job->bodies[solid_body_id]->nodes->f[i]
                                 + f_i[i]*job->grid->nodeVolume(job,i)) / 2.0;
    }
    k2 = job->dt * F(job, driver, (u_0 + k1/2.0), tmp_mx_t);
    f_i2 = f_i;

    for (int i=0; i<tmp_mx_t.size(); i++){
        tmp_mx_t(i) = job->bodies[solid_body_id]->nodes->mx_t[i]
                      + job->dt*(job->bodies[solid_body_id]->nodes->f[i]
                                 + f_i[i]*job->grid->nodeVolume(job,i)) / 2.0;
    }
    k3 = job->dt * F(job, driver, (u_0 + k2/2.0), tmp_mx_t);
    f_i3 = f_i;

    for (int i=0; i<tmp_mx_t.size(); i++){
        tmp_mx_t(i) = job->bodies[solid_body_id]->nodes->mx_t[i]
                      + job->dt*(job->bodies[solid_body_id]->nodes->f[i]
                                 + f_i[i]*job->grid->nodeVolume(job,i));
    }
    k4 = job->dt * F(job, driver, (u_0 + k3), tmp_mx_t);
    f_i4 = f_i;

    u_n = u_0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    f_i = (f_i1 + 2.0*f_i2 + 2.0*f_i3 + f_i4)/6.0;

    convertVectorToStateSpace(job,
                              driver,
                              u_n,
                              driver->fluid_body->rho,
                              driver->fluid_body->p,
                              driver->fluid_body->rhoE);

    //add force contribution to solid phase
    for (int i=0; i<job->grid->node_count; i++) {
        job->bodies[solid_body_id]->nodes->f[i] += f_i[i]*job->grid->nodeVolume(job,i);
    }

    /*----------------------*/
    /*   End FVM-RK4 Step   */
    /*----------------------*/

    //enforce boundary conditions AGAIN
    generateBoundaryConditions(job);
    addBoundaryConditions(job);

    //move grid
    moveGrid(job);

    //move particles
    movePoints(job);

    //calculate strainrate
    calculateStrainRate(job);

    //update density
    updateDensity(job, Material::UPDATE);

    //add body forces
    job->driver->generateGravity(job);
    job->driver->applyGravity(job);

    //update stress
    updateStress(job, Material::UPDATE);

    /*----------------------*/
    /*     End MPM Step     */
    /*----------------------*/


    return;
}

/*----------------------------------------------------------------------------*/
//fvm conversion f'ns
Eigen::VectorXd FVMMixtureSolverRK4::F(Job* job,
                                       FiniteVolumeDriver* driver,
                                       const Eigen::VectorXd& u,
                                       const KinematicVectorArray& mx_t){

    //convert mx_t to x_t
    double m;
    for (int i=0; i<mx_t.size(); i++){
        m = job->bodies[solid_body_id]->nodes->m(i);
        if (m > 0) {
            driver->fluid_body->v_s(i) = mx_t(i) / m;
        } else {
            driver->fluid_body->v_s(i).setZero();
        }
    }

    //map to fvm
    driver->fluid_grid->mapSolidVelocityToQuadraturePoints(job, driver);

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

    //fluid fluxes associated with each volume
    density_fluxes = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
    momentum_fluxes = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2
    energy_fluxes = driver->fluid_grid->calculateElementEnergyFluxes(job, driver);      //J/s?

    //get discretized inter-phase force
    f_i = driver->fluid_grid->calculateInterphaseForces(job, driver);

    //map interphase force to fvm elements
    f_e = driver->fluid_grid->M * f_i;

    //add gravity and scale by element volume
    double volume;
    for (int e = 0; e < driver->fluid_grid->element_count; e++) {
        volume = driver->fluid_grid->getElementVolume(e);
        density_fluxes(e) /= volume;
        momentum_fluxes(e) /= volume;
        energy_fluxes(e) /= volume;
        momentum_fluxes[e] += driver->fluid_body->rho(e) * driver->gravity;
        momentum_fluxes[e] -= f_e[e] / volume;
    }

    //convert resulting fluxes to output
    Eigen::VectorXd result = Eigen::VectorXd(vector_size);
    convertStateSpaceToVector(job, driver, result, density_fluxes, momentum_fluxes, energy_fluxes);

    return result;
}

void FVMMixtureSolverRK4::convertVectorToStateSpace(Job* job, FiniteVolumeDriver* driver,
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
    offset = driver->fluid_grid->element_count * (1 + GRID_DIM);
    if (driver->TYPE == FiniteVolumeDriver::THERMAL){
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            rhoE(e) = v(offset + e);
        }
        offset = driver->fluid_grid->element_count * (2 + GRID_DIM);
    }

    return;
}


void FVMMixtureSolverRK4::convertStateSpaceToVector(Job* job, FiniteVolumeDriver* driver,
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

    offset = driver->fluid_grid->element_count * (1 + GRID_DIM);

    //if thermal simulation, fill energy values
    if (driver->TYPE == FiniteVolumeDriver::THERMAL){
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            v(offset + e) = rhoE(e);
        }
        offset = driver->fluid_grid->element_count * (2 + GRID_DIM);
    }

    return;
}

/*----------------------------------------------------------------------------*/
//functions to avoid updating stress

void FVMMixtureSolverRK4::updateDensity(Job* job, int SPEC){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (int i=0;i<job->bodies[b]->points->v.rows();i++) {
            job->bodies[b]->points->v(i) *= std::exp(job->dt * job->bodies[b]->points->L(i).trace());
        }

        if (SPEC == Material::UPDATE) {
            //this is new, but maybe useful
            job->bodies[b]->points->updateIntegrators(job, job->bodies[b].get());
        }
    }
    return;
}

void FVMMixtureSolverRK4::updateStress(Job* job, int SPEC){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b]->activeMaterial == 0){
            continue;
        }
        job->bodies[b]->material->calculateStress(job, job->bodies[b].get(), SPEC);
    }
    return;
}
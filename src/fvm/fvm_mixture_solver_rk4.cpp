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
        //rho, rho*u
        vector_size = driver->fluid_grid->element_count*(1 + driver->fluid_grid->GRID_DIM);
    } else if (driver->TYPE == FiniteVolumeDriver::THERMAL){
        //rho, rho*u, rho*E
        vector_size = driver->fluid_grid->element_count*(2 + driver->fluid_grid->GRID_DIM);
    } else {
        std::cerr << "ERROR! FVMSteadyStateSolver not implemented for simulationg TYPE " << driver->TYPE << "! Exiting." << std::endl;
        exit(0);
    }

    //initialize vectors
    f_d = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    f_b = f_d;
    f_i = f_d;
    first_drag_correction = f_d;
    second_drag_correction = f_d;
    f_d_e = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    f_b_e = f_d_e;
    f_e = f_d_e;

    //ask MPM to map properties to grid for FVM initialization
    job->bodies[solid_body_id]->generateMap(job, cpdi_spec);
    job->bodies[solid_body_id]->nodes->m = job->bodies[solid_body_id]->S * job->bodies[solid_body_id]->points->m;

    //loop over str-props and assign relevant flags
    std::vector<std::string> options = {"USE_THERMAL_FLOOR"};
    for (int i=0; i<str_props.size(); i++){
        switch (Parser::findStringID(options, str_props[i])){
            case 0:
                //USE_THERMAL_FLOOR
                USE_THERMAL_FLOOR = true;
                T_floor = fp64_props[0];
                cv_floor = fp64_props[1];
                std::cout << "Using thermal floor. WARNING: energy WILL be inserted into simulation. T = " << T_floor;
                std::cout << " [K]. cv = " << cv_floor << " J/kg*K." << std::endl;
                break;
            default:
                //do nothing
                break;
        }
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

    //add body forces
    driver->generateGravity(job);
    driver->applyGravity(job);
    //add specific body force to points
    for (int i=0; i<job->bodies[solid_body_id]->points->b.size(); i++){
        job->bodies[solid_body_id]->points->b[i] += job->bodies[solid_body_id]->points->v(i)/job->bodies[solid_body_id]->points->m(i)
                                                   * driver->getSolidLoading(job, job->bodies[solid_body_id]->points->x(i));
    }

    //map particles to grid
    mapPointsToNodes(job);

    //enforce boundary conditions on solid body only
    f_bc = -job->bodies[solid_body_id]->nodes->f;
    job->bodies[solid_body_id]->boundary->generateRules(job,job->bodies[solid_body_id].get());
    job->bodies[solid_body_id]->boundary->applyRules(job,job->bodies[solid_body_id].get());
    //store difference in forces
    f_bc += job->bodies[solid_body_id]->nodes->f;

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

    //get drag and corrected coefficients from initial state
    //driver->fluid_grid->calculateSplitIntegralDragForces(job, driver, f_d, f_d_e);
    K_n = driver->fluid_grid->getCorrectedDragCoefficients(job, driver);

    /*----------------------*/
    /* Estimate k1 w/o drag */
    /*----------------------*/

    f_d.setZero();
    f_d_e.setZero();
    k1 = job->dt * F(job, driver, u_0);
    f_i1 = f_i;                             //store intermediate drag calculations
    f_b_e = f_e/6.0;                        //store intermediate buoyancy term

    /*-----------------------*/
    /* Estimate Drag at t+dt */
    /*-----------------------*/

    //assign u_0 + k1 to fluid phase
    //assign v_0 + f_i to solid phase
    convertVectorToStateSpace(job,
                              driver,
                              (u_0 + k1),
                              driver->fluid_body->rho,
                              driver->fluid_body->p,
                              driver->fluid_body->rhoE);

    for (int i = 0; i<job->bodies[solid_body_id]->nodes->mx_t.size(); i++){
        if (job->bodies[solid_body_id]->nodes->m(i) > 0) {
            job->bodies[solid_body_id]->nodes->x_t(i) = (job->bodies[solid_body_id]->nodes->mx_t(i) +
                                                         job->dt*(job->bodies[solid_body_id]->nodes->f(i)
                                                                  + f_i1(i)*job->grid->nodeVolume(job,i)))
                                                        / job->bodies[solid_body_id]->nodes->m(i);
        } else {
            job->bodies[solid_body_id]->nodes->x_t(i).setZero();
        }
    }

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

    //reset v_0
    for (int i = 0; i<job->bodies[solid_body_id]->nodes->mx_t.size(); i++){
        if (job->bodies[solid_body_id]->nodes->m(i) > 0) {
            job->bodies[solid_body_id]->nodes->x_t(i) = job->bodies[solid_body_id]->nodes->mx_t(i)
                                                        / job->bodies[solid_body_id]->nodes->m(i);
        } else {
            job->bodies[solid_body_id]->nodes->x_t(i).setZero();
        }
    }

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

    f_i1 += f_d;

    /*----------------------*/
    /*    Find k2 to k4     */
    /*----------------------*/

    k2 = job->dt * F(job, driver, (u_0 + k1/2.0));
    f_i2 = f_i;
    f_b_e += (f_e - f_d_e)/3.0;

    k3 = job->dt * F(job, driver, (u_0 + k2/2.0));
    f_i3 = f_i;
    f_b_e += (f_e - f_d_e)/3.0;

    k4 = job->dt * F(job, driver, (u_0 + k3));
    f_i4 = f_i;
    f_b_e += (f_e - f_d_e)/6.0;

    u_n = u_0 + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    f_i = (f_i1 + 2.0*f_i2 + 2.0*f_i3 + f_i4)/6.0;

    convertVectorToStateSpace(job,
                              driver,
                              u_n,
                              driver->fluid_body->rho,
                              driver->fluid_body->p,
                              driver->fluid_body->rhoE);

    //add buoyant force contribution to solid phase
    for (int i=0; i<job->grid->node_count; i++) {
        job->bodies[solid_body_id]->nodes->f[i] += (f_i[i] - f_d[i])*job->grid->nodeVolume(job,i);
    }

    /*-----------------------------*/
    /* Ensure Thermal Floor is Set */
    /*-----------------------------*/

    //if using Thermal Floor:
    double tmpT;
    if (driver->TYPE == driver->THERMAL && USE_THERMAL_FLOOR){
        // loop over elements
        for (int e=0; e<driver->fluid_grid->element_count; e++) {

            //check thermal floor
            tmpT = driver->fluid_material->getTemperature(job,
                                                          driver,
                                                          driver->fluid_body->rho(e),
                                                          driver->fluid_body->p(e),
                                                          driver->fluid_body->rhoE(e),
                                                          driver->fluid_body->n_e(e));
            if (tmpT <= T_floor || !std::isfinite(tmpT)) {

                //Highlight Error
                std::cout << std::endl;
                std::cout << "    Thermal Floor: [" << e << "]: " << driver->fluid_material->getTemperature(job,
                                                                                                            driver,
                                                                                                            driver->fluid_body->rho(e),
                                                                                                            driver->fluid_body->p(e),
                                                                                                            driver->fluid_body->rhoE(e),
                                                                                                            driver->fluid_body->n_e(e));
                std::cout << " --> " << T_floor << std::endl;

                //Assign New Temperature
                driver->fluid_body->rhoE(e) = driver->fluid_body->rho(e) * cv_floor * T_floor
                                              + 0.5 * driver->fluid_body->p[e].dot(driver->fluid_body->p[e]) / driver->fluid_body->rho(e);

            }
        }
    }

    /*-----------------------*/
    /* Estimate Drag at t+dt */
    /*-----------------------*/

    //assign v_plus grid velocity
    for (int i = 0; i<job->bodies[solid_body_id]->nodes->mx_t.size(); i++){
        if (job->bodies[solid_body_id]->nodes->m(i) > 0) {
            job->bodies[solid_body_id]->nodes->x_t(i) = (job->bodies[solid_body_id]->nodes->mx_t(i) +
                                                         job->dt*(job->bodies[solid_body_id]->nodes->f(i)))
                                                        / job->bodies[solid_body_id]->nodes->m(i);
        } else {
            job->bodies[solid_body_id]->nodes->x_t(i).setZero();
        }
    }

    //subtract predicted drag force from fluid to get u_plus
    for (int e=0; e<driver->fluid_grid->element_count; e++) {
        driver->fluid_body->p[e] += job->dt * f_d_e[e];
    }

    //tell fluid grid to update mixture properties
    driver->fluid_grid->mapMixturePropertiesToQuadraturePoints(job, driver);
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

    //add correction to solid momentum and reset velocity
    for (int i = 0; i<job->bodies[solid_body_id]->nodes->mx_t.size(); i++){
        if (job->bodies[solid_body_id]->nodes->m(i) > 0) {
            job->bodies[solid_body_id]->nodes->x_t(i) = job->bodies[solid_body_id]->nodes->mx_t(i)
                                                         / job->bodies[solid_body_id]->nodes->m(i);

            //only non-zero here
            job->bodies[solid_body_id]->nodes->f[i] += second_drag_correction[i]*job->grid->nodeVolume(job,i);
        }
    }

    //remove original BC forces from solid body
    job->bodies[solid_body_id]->nodes->f -= f_bc;

    /*----------------------*/
    /*   End FVM-RK4 Step   */
    /*----------------------*/

    //add arbitrary loading conditions
    generateLoads(job);
    applyLoads(job);

    //add contact forces
    generateContacts(job);
    addContacts(job);

    //enforce boundary conditions
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

    //update stress
    updateStress(job, Material::UPDATE);

    /*----------------------*/
    /*     End MPM Step     */
    /*----------------------*/


    return;
}

/*----------------------------------------------------------------------------*/
//fvm conversion f'ns
Eigen::VectorXd FVMMixtureSolverRK4::F(Job* job, FiniteVolumeDriver* driver, const Eigen::VectorXd& u){
    //convert input vector into fluid_body state space
    convertVectorToStateSpace(job, driver,
                              u,
                              driver->fluid_body->rho,
                              driver->fluid_body->p,
                              driver->fluid_body->rhoE);

    //if using Thermal Floor:
    double tmpT;
    if (driver->TYPE == driver->THERMAL && USE_THERMAL_FLOOR){
        // loop over elements
        for (int e=0; e<driver->fluid_grid->element_count; e++) {

            //check thermal floor
            tmpT = driver->fluid_material->getTemperature(job,
                                                       driver,
                                                       driver->fluid_body->rho(e),
                                                       driver->fluid_body->p(e),
                                                       driver->fluid_body->rhoE(e),
                                                       driver->fluid_body->n_e(e));
            if (tmpT <= T_floor || !std::isfinite(tmpT)) {

                //Highlight Error
                std::cout << std::endl;
                std::cout << "    Thermal Floor: [" << e << "]: " << driver->fluid_material->getTemperature(job,
                                                                                                            driver,
                                                                                                            driver->fluid_body->rho(e),
                                                                                                            driver->fluid_body->p(e),
                                                                                                            driver->fluid_body->rhoE(e),
                                                                                                            driver->fluid_body->n_e(e));
                std::cout << " --> " << T_floor << std::endl;

                //Assign New Temperature
                driver->fluid_body->rhoE(e) = driver->fluid_body->rho(e) * cv_floor * T_floor
                                                    + 0.5 * driver->fluid_body->p[e].dot(driver->fluid_body->p[e]) / driver->fluid_body->rho(e);

            }
        }
    }

    //calculate fluxes
    driver->fluid_grid->constructDensityField(job, driver);
    driver->fluid_grid->constructMomentumField(job, driver);
    if (driver->TYPE == driver->THERMAL) {
        driver->fluid_grid->constructEnergyField(job, driver);
    }
    driver->fluid_grid->constructPorosityField(job, driver);
    driver->fluid_grid->constructVelocityField(job, driver);

    //fluid fluxes associated with each volume
    density_fluxes = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
    momentum_fluxes = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2
    if (driver->TYPE == driver->THERMAL) {
        energy_fluxes = driver->fluid_grid->calculateElementEnergyFluxes(job, driver);      //J/s?
    }

    //get discretized inter-phase force
    driver->fluid_grid->calculateSplitIntegralBuoyantForces(job, driver, f_b, f_e);
    f_i = f_b + f_d;
    f_e += f_d_e;

    //add gravity and scale by element volume
    double volume;
    KinematicVector load;
    for (int e = 0; e < driver->fluid_grid->element_count; e++) {
        volume = driver->fluid_grid->getElementVolume(e);
        load = driver->getFluidLoading(job, driver->fluid_grid->getElementCentroid(job, e));
        density_fluxes(e) /= volume;
        momentum_fluxes[e] /= volume;
        if (driver->TYPE == driver->THERMAL) {
            energy_fluxes(e) /= volume;
        }
        momentum_fluxes[e] += driver->fluid_body->rho(e) * driver->gravity;
        momentum_fluxes[e] += load;
        momentum_fluxes[e] -= f_e[e];// / volume;
        if (driver->TYPE == driver->THERMAL) {
            energy_fluxes(e) += driver->gravity.dot(driver->fluid_body->p[e]);
            energy_fluxes(e) += load.dot(driver->fluid_body->p[e])/driver->fluid_body->rho(e);
        }
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

/*------------------------------------------------------------------------------*/
void FVMMixtureSolverRK4::writeFrame(Job *job, FiniteVolumeDriver *driver) {
    driver->serializer->writeVectorArray(f_d_e, "drag_force");
    driver->serializer->writeVectorArray(f_b_e, "buoyant_force");
    return;
}
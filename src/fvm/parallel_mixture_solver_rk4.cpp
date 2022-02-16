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
#include "objects/solvers/solvers.hpp"

void ParallelMixtureSolverRK4::init(Job* job, FiniteVolumeDriver* driver){

    //copy properties to parents to avoid ambiguity?
    ThreadPoolExplicitUSL::int_props = FVMMixtureSolverRK4::int_props;
    ThreadPoolExplicitUSL::fp64_props = FVMMixtureSolverRK4::fp64_props;
    ThreadPoolExplicitUSL::str_props = FVMMixtureSolverRK4::str_props;

    //call parent initialization f'ns
    FVMMixtureSolverRK4::init(job, driver);
    ThreadPoolExplicitUSL::init(job);

    return;
}

/*----------------------------------------------------------------------------*/
//one forward mixed time-stpe

void ParallelMixtureSolverRK4::step(Job* job, FiniteVolumeDriver* driver){

    //tell mixture solver base class to take step
    //FVMMixtureSolverRK4::step(job, driver);


    /*----------------------*/
    /*  Begin FVM-MPM Step  */
    /*----------------------*/

    //create map
    FVMMixtureSolverRK4::createMappings(job);

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
    FVMMixtureSolverRK4::generateLoads(job);
    FVMMixtureSolverRK4::applyLoads(job);

    //add contact forces
    FVMMixtureSolverRK4::generateContacts(job);
    FVMMixtureSolverRK4::addContacts(job);

    //enforce boundary conditions
    FVMMixtureSolverRK4::generateBoundaryConditions(job);
    FVMMixtureSolverRK4::addBoundaryConditions(job);

    //move grid
    moveGrid(job);

    //move particles
    movePoints(job);

    //calculate strainrate
    calculateStrainRate(job);

    //update density
    FVMMixtureSolverRK4::updateDensity(job, Material::UPDATE);

    //update stress
    FVMMixtureSolverRK4::updateStress(job, Material::UPDATE);

    /*----------------------*/
    /*     End MPM Step     */
    /*----------------------*/

    return;
}

/*---------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ParallelMixtureSolverRK4::mapPointsToNodes(Job* job){
    //check that threadpool exists and not serial
    if (job->thread_count <= 1){
        ExplicitUSL::mapPointsToNodes(job);
        return;
    }

    Body *body;
    Points *points;
    Nodes *nodes;
    Eigen::VectorXd pval;
    Eigen::VectorXd nval;
    KinematicVectorArray pvec;
    MaterialVectorArray nvec;
    MaterialTensorArray tmpMat;
    double tmpVAL;

    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = job->bodies[b].get();
        points = job->bodies[b]->points.get();
        nodes = job->bodies[b]->nodes.get();

        //map mass
        parallelMultiply(body->S, points->m, nodes->m, MPMScalarSparseMatrix::NORMAL, true, b);

        //map momentum
        parallelMultiply(points->x_t, points->m, 1.0, points->mx_t, true);

        parallelMultiply(body->S, points->mx_t, nodes->mx_t, MPMScalarSparseMatrix::NORMAL, true, b);

        parallelDivide(nodes->mx_t, nodes->m, 1.0, nodes->x_t, true);

        //map body force
        pvec = KinematicVectorArray(points->b.size(), points->b.VECTOR_TYPE);
        parallelMultiply(points->b, points->m, 1.0, pvec, true);
        parallelMultiply(body->S, pvec, nodes->f, MPMScalarSparseMatrix::NORMAL, true, b);

        //map divergence of stress
        tmpMat = MaterialTensorArray(points->T.size());
        parallelMultiply(points->T, points->v, 1.0, tmpMat, true);

        nvec = MaterialVectorArray(nodes->x.size());
        parallelMultiply(body->gradS, tmpMat, nvec, MPMSparseMatrixBase::NORMAL, true, b);
        for (int i = 0; i < nvec.size(); i++) {
            nodes->f[i] -= KinematicVector(nvec[i], nodes->f.VECTOR_TYPE);
        }

        if (job->JOB_TYPE == job->JOB_AXISYM){
            pval = Eigen::VectorXd(points->x.size());
            nval = Eigen::VectorXd(nodes->x.size());

            //scale s_tt by r and add contribution to f_r
            for (int i=0;i<pval.rows();i++){
                pval(i) = tmpMat(i,2,2) / points->x(i,0);
            }
            parallelMultiply(body->S, pval, nval, MPMScalarSparseMatrix::NORMAL, true, b);
            for (int i=0; i<nval.rows(); i++){
                nodes->f(i,0) -= nval(i);
            }

            //scale s_rt by r and add contribution to f_t
            for (int i=0; i<pval.rows(); i++){
                pval(i) = tmpMat(i,0,2) / points->x(i,0);
            }
            parallelMultiply(body->S, pval, nval, MPMScalarSparseMatrix::NORMAL, true, b);
            for (int i=0; i<nval.rows(); i++){
                nodes->f(i,2) += nval(i);
            }
        }
        //std::cout << points->m.sum() << " ?= " << nodes->m.sum() << std::endl;
    }
    return;
}

void ParallelMixtureSolverRK4::moveGrid(Job* job){
    //check that threadpool exists and not serial
    if (job->thread_count <= 1){
        ExplicitUSL::moveGrid(job);
        return;
    }

    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //update momentum
        job->bodies[b]->nodes->mx_t += job->dt * job->bodies[b]->nodes->f;

        //calculate velocity
        parallelDivide(job->bodies[b]->nodes->mx_t, job->bodies[b]->nodes->m, 1.0, job->bodies[b]->nodes->x_t, true);

        //set displacement
        job->bodies[b]->nodes->u = job->dt * job->bodies[b]->nodes->x_t;

        //calculate difference in velocity
        parallelDivide(job->bodies[b]->nodes->f, job->bodies[b]->nodes->m, job->dt, job->bodies[b]->nodes->diff_x_t, true);
    }
    return;
}

void ParallelMixtureSolverRK4::movePoints(Job* job){
    //check that threadpool exists and not serial
    if (job->thread_count <= 1){
        ExplicitUSL::movePoints(job);
        return;
    }

    Body* body;
    Points* points;
    Nodes* nodes;
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = job->bodies[b].get();
        points = job->bodies[b]->points.get();
        nodes = job->bodies[b]->nodes.get();

        //map nodal displacement to point positions
        parallelMultiply(body->S, nodes->u, points->x, MPMScalarSparseMatrix::TRANSPOSED, false, b);
        parallelMultiply(body->S, nodes->u, points->u, MPMScalarSparseMatrix::TRANSPOSED, false, b);

        //fix position for out of plane dimension
        if (job->grid->GRID_DIM < job->DIM){
            for (int i=0; i<points->x.size(); i++){
                for (int pos=job->grid->GRID_DIM; pos<job->DIM; pos++){
                    points->x(i,pos) = 0;
                }
            }
        }

        //map nodal velocity diff to points
        parallelMultiply(body->S, nodes->diff_x_t, points->x_t, MPMSparseMatrixBase::TRANSPOSED, false, b);

        //calculate momentum
        parallelMultiply(points->x_t, points->m, 1.0, points->mx_t, true);
    }
    return;
}

void ParallelMixtureSolverRK4::calculateStrainRate(Job* job){
    //check that threadpool exists and not serial
    if (job->thread_count <= 1){
        ExplicitUSL::calculateStrainRate(job);
        return;
    }

    Body* body;
    Points* points;
    Nodes* nodes;
    KinematicVectorArray pvec;

    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = job->bodies[b].get();
        points = job->bodies[b]->points.get();
        nodes = job->bodies[b]->nodes.get();

        //calculate gradient of nodal velocity at points
        parallelMultiply(body->gradS, nodes->x_t, points->L, MPMSparseMatrixBase::TRANSPOSED, true, b);

        //correct L if axisymmetric
        if (job->JOB_TYPE == job->JOB_AXISYM){
            pvec = KinematicVectorArray(points->L.size(), points->L.TENSOR_TYPE);
            parallelMultiply(body->S, nodes->x_t, pvec, MPMSparseMatrixBase::TRANSPOSED, true, b);
            for (int i=0; i<points->L.size(); i++){
                points->L(i,0,2) = -pvec(i,2) / points->x(i,0);
                points->L(i,2,2) = pvec(i,0) / points->x(i,0);
            }
        }
    }
    return;
}

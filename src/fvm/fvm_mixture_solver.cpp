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
#include "objects/bodies/bodies.hpp"

void FVMMixtureSolver::init(Job* job, FiniteVolumeDriver* driver){
    //check that contact properties are set
    if (int_props.size() == 0){
        //do nothing
        cpdi_spec = DefaultBody::CPDI_ON;
        contact_spec = Contact::IMPLICIT;
    } if (int_props.size() == 1){
        //cpdi_spec given as argument
        cpdi_spec = int_props[0];
        contact_spec = Contact::IMPLICIT;
    } if (int_props.size() >= 2){
        //cpdi_spec and contact_spec given
        cpdi_spec = int_props[0];
        contact_spec = int_props[1];
    }

    if (str_props.size() < 1){
        std::cout << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 1 property defined (solid_body).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //set body ids by name
        if (str_props.size() == 1){
            for (int b = 0; b < job->bodies.size(); b++) {
                if (str_props[0].compare(job->bodies[b]->name) == 0){
                    solid_body_id = b;
                    break;
                }
            }
        }

        // or set body ids by int
        if (solid_body_id < 0){
            std::cout << str_props.size() << "\n";
            fprintf(stderr,
                    "%s:%s: Need at least 1 property defined (solid_body).\n",
                    __FILE__, __func__);
            exit(0);
        }
    }

    printf("Solver properties (cpdi_spec = %i, contact_spec = %i, solid_body_id = %i).\n", cpdi_spec, contact_spec, solid_body_id);

    //size flux containers
    density_fluxes = Eigen::VectorXd(driver->fluid_grid->element_count);
    momentum_fluxes = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
    energy_fluxes = Eigen::VectorXd(driver->fluid_grid->element_count);

    f_i = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
    f_e = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);

    //done
    std::cout << "FiniteVolumeSolver initialized." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//one forward mixed time-stpe

void FVMMixtureSolver::step(Job* job, FiniteVolumeDriver* driver){

    /*----------------------*/
    /*  Begin FVM-MPM Step  */
    /*----------------------*/

    //create map
    createMappings(job);

    //map particles to grid
    mapPointsToNodes(job);

    //construct finite volume gradients
    calculateElementGradients(job, driver);

    //map mixture properties to finite volumes
    mapMixturePropertiesToElements(job, driver);

    //generate fluxes (don't add yet)
    generateFluxes(job, driver);

    //generate mixture forces (FVM only)
    //generateMixtureForces(job, driver);

    //add fluxes and forces (FVM and MPM)
    applyFluxes(job, driver);
    //applyMixtureForces(job, driver);

    return;

    /*----------------------*/
    /*     End FVM Step     */
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
    updateDensity(job);

    //add body forces
    job->driver->generateGravity(job);
    job->driver->applyGravity(job);

    //update stress
    updateStress(job);

    /*----------------------*/
    /*     End MPM Step     */
    /*----------------------*/

    return;
}

/*----------------------------------------------------------------------------*/
//
void FVMMixtureSolver::calculateElementGradients(Job* job, FiniteVolumeDriver* driver){

    if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {

        //call grid to reconstruct conserved fields
        driver->fluid_grid->constructDensityField(job, driver);
        driver->fluid_grid->constructMomentumField(job, driver);
        driver->fluid_grid->constructPorosityField(job, driver);

    } else if (driver->TYPE == FiniteVolumeDriver::THERMAL) {

        //will need to implement this later
        std::cerr << "ERROR! Missing implementation of mixture energy fluxes. Exiting." << std::endl;
        exit(0);

    } else {

        //shouldn't have another flag...
        std::cerr << "ERROR! FVMMixtureSolver does not have simulation TYPE " << driver->TYPE << " implemented! Exiting." << std::endl;
        exit(0);
    }

    return;
}

void FVMMixtureSolver::mapMixturePropertiesToElements(Job* job, FiniteVolumeDriver* driver){
    //call grid to map mixture properties from MPM to FVM
    driver->fluid_grid->mapMixturePropertiesToQuadraturePoints(job, driver);
    return;
}

void FVMMixtureSolver::generateFluxes(Job* job, FiniteVolumeDriver* driver){

    if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
        //fluid fluxes associated with each volume
        density_fluxes = driver->fluid_grid->calculateElementMassFluxes(job, driver);       //kg/s
        momentum_fluxes = driver->fluid_grid->calculateElementMomentumFluxes(job, driver);  //kg m/s^2

    } else if (driver->TYPE == FiniteVolumeDriver::THERMAL) {
        //will need to implement this later
        std::cerr << "ERROR! Missing implementation of mixture energy fluxes. Exiting." << std::endl;
        exit(0);

    } else {
        //shouldn't have another flag...
        std::cerr << "ERROR! FVMMixtureSolver does not have simulation TYPE " << driver->TYPE << " implemented! Exiting." << std::endl;
        exit(0);
    }

    return;
}

void FVMMixtureSolver::applyFluxes(Job* job, FiniteVolumeDriver* driver){
    if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
        //use fluxes to update element-wise average density and momentum (Forward Euler)
        double volume;
        for (int e = 0; e < driver->fluid_grid->element_count; e++) {
            volume = driver->fluid_grid->getElementVolume(e);
            driver->fluid_body->rho(e) += density_fluxes(e) / volume * job->dt;     //d(rho dv)/dt   = flux
            driver->fluid_body->p[e] += momentum_fluxes[e] / volume * job->dt;    //d(rho u dv)/dt = flux

            //add gravitational contribution to momentum update
            driver->fluid_body->p[e] += driver->fluid_body->rho(e) * driver->gravity * job->dt;
        }

    } else if (driver->TYPE == FiniteVolumeDriver::THERMAL) {
        //will need to implement this later
        std::cerr << "ERROR! Missing implementation of mixture energy fluxes. Exiting." << std::endl;
        exit(0);

    } else {
        //shouldn't have another flag...
        std::cerr << "ERROR! FVMMixtureSolver does not have simulation TYPE " << driver->TYPE << " implemented! Exiting." << std::endl;
        exit(0);
    }

    return;
}

void FVMMixtureSolver::generateMixtureForces(Job* job, FiniteVolumeDriver* driver){
    //get discretized inter-phase force
    f_i = driver->fluid_grid->calculateInterphaseForces(job, driver);
    return;
}

void FVMMixtureSolver::applyMixtureForces(Job* job, FiniteVolumeDriver* driver){
    //subtract force contribution from fluid phase
    f_e = driver->fluid_grid->M * f_i;
    for (int e=0; e<driver->fluid_grid->element_count; e++){
        driver->fluid_body->p[e] -= f_e[e] / driver->fluid_grid->getElementVolume(e) * job->dt;
    }

    //add force contribution to solid phase
    for (int i=0; i<job->grid->node_count; i++) {
        job->bodies[solid_body_id]->nodes->f[i] += f_i[i]*job->grid->nodeVolume(job,i);
    }
}


/*----------------------------------------------------------------------------*/
//standard MPM time-stepping functions
void FVMMixtureSolver::createMappings(Job *job){
    for (int b=0;b<job->bodies.size();b++){
        //job->bodies[b]->generateMap(job, DefaultBody::CPDI_ON); //use_cpdi by default
        job->bodies[b]->generateMap(job, cpdi_spec);
    }
    return;
}

void FVMMixtureSolver::mapPointsToNodes(Job* job){
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
        nodes->m = body->S * points->m; //m_i = S_ip * m_p

        //map momentum
        for (int i = 0; i < points->mx_t.size(); i++) {
            points->mx_t(i) = points->m(i) * points->x_t(i);
        }
        nodes->mx_t = body->S * points->mx_t;

        //calculate velocity
        for (int i = 0; i < nodes->x_t.size(); i++) {
            if (nodes->m(i) > 0) {
                nodes->x_t(i) = nodes->mx_t(i) / nodes->m(i);
            } else {
                nodes->x_t(i).setZero();
            }
        }

        //map body force
        pvec = KinematicVectorArray(points->b.size(), points->b.VECTOR_TYPE);
        for (int i = 0; i < points->b.size(); i++) {
            pvec(i) = points->m(i) * points->b(i);
        }
        nodes->f = body->S * pvec;

        //map divergence of stress
        tmpMat = points->T;
        for (int i = 0; i < tmpMat.size(); i++) {
            tmpMat[i] *= points->v[i];
        }
        nvec = body->gradS.left_multiply_by_tensor(tmpMat); //f_i = v_p T_p * gradS_ip
        for (int i = 0; i < nvec.size(); i++) {
            nodes->f[i] -= KinematicVector(nvec[i], nodes->f.VECTOR_TYPE);
        }

        if (job->JOB_TYPE == job->JOB_AXISYM){
            pval = Eigen::VectorXd(points->x.size());

            //scale s_tt by r and add contribution to f_r
            for (int i=0;i<pval.rows();i++){
                pval(i) = tmpMat(i,2,2) / points->x(i,0);
            }
            nval = body->S * pval;
            for (int i=0; i<nval.rows(); i++){
                nodes->f(i,0) -= nval(i);
            }

            //scale s_rt by r and add contribution to f_t
            for (int i=0; i<pval.rows(); i++){
                pval(i) = tmpMat(i,0,2) / points->x(i,0);
            }
            nval = body->S * pval;
            for (int i=0; i<nval.rows(); i++){
                nodes->f(i,2) += nval(i);
            }
        }
        //std::cout << points->m.sum() << " ?= " << nodes->m.sum() << std::endl;
    }
    return;
}

void FVMMixtureSolver::generateContacts(Job* job){
    for (int c=0;c<job->contacts.size();c++){
        if (job->activeContacts[c] == 0){
            continue;
        }
        job->contacts[c]->generateRules(job);
    }
    return;
}

void FVMMixtureSolver::addContacts(Job* job){
    for (int c=0;c<job->contacts.size();c++){
        if (job->activeContacts[c] == 0){
            continue;
        }
        //job->contacts[c]->applyRules(job,Contact::EXPLICIT);
        job->contacts[c]->applyRules(job, contact_spec);
    }
    return;
}

void FVMMixtureSolver::generateBoundaryConditions(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b]->activeBoundary == 0){
            continue;
        }
        job->bodies[b]->boundary->generateRules(job,job->bodies[b].get());
    }
    return;
}

void FVMMixtureSolver::addBoundaryConditions(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b]->activeBoundary == 0){
            continue;
        }
        job->bodies[b]->boundary->applyRules(job,job->bodies[b].get());
    }
    return;
}

void FVMMixtureSolver::moveGrid(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //update momentum
        job->bodies[b]->nodes->mx_t += job->dt * job->bodies[b]->nodes->f;

        //calculate velocity
        for (int i=0;i<job->bodies[b]->nodes->x_t.size();i++){
            if (job->bodies[b]->nodes->m(i) > 0) {
                job->bodies[b]->nodes->x_t(i) = job->bodies[b]->nodes->mx_t(i) / job->bodies[b]->nodes->m(i);
            } else {
                job->bodies[b]->nodes->x_t(i).setZero();
            }
        }

        //set displacement
        job->bodies[b]->nodes->u = job->dt * job->bodies[b]->nodes->x_t;

        //calculate difference in velocity
        for (int i=0;i<job->bodies[b]->nodes->diff_x_t.size();i++){
            if (job->bodies[b]->nodes->m(i) > 0) {
                job->bodies[b]->nodes->diff_x_t(i) = job->dt * job->bodies[b]->nodes->f(i) / job->bodies[b]->nodes->m(i);
            } else {
                job->bodies[b]->nodes->diff_x_t(i).setZero();
            }
        }
    }
    return;
}

void FVMMixtureSolver::movePoints(Job* job){
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
        points->x += body->S.operate(nodes->u, MPMSparseMatrixBase::TRANSPOSED);
        points->u += body->S.operate(nodes->u, MPMSparseMatrixBase::TRANSPOSED);

        //fix position for out of plane dimension
        if (job->grid->GRID_DIM < job->DIM){
            for (int i=0; i<points->x.size(); i++){
                for (int pos=job->grid->GRID_DIM; pos<job->DIM; pos++){
                    points->x(i,pos) = 0;
                }
            }
        }

        //map nodal velocity diff to points
        points->x_t += body->S.operate(nodes->diff_x_t, MPMSparseMatrixBase::TRANSPOSED);

        //calculate momentum
        for (int i=0;i<points->mx_t.size();i++){
            points->mx_t(i) = points->m(i) * points->x_t(i);
        }
    }
    return;
}

void FVMMixtureSolver::calculateStrainRate(Job* job){
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
        //body->bodyCalcPointGradient(job,points->L,nodes->x_t,Body::SET);
        points->L = body->gradS.tensor_product_transpose(nodes->x_t, MPMSparseMatrixBase::TRANSPOSED);

        //correct L if axisymmetric
        if (job->JOB_TYPE == job->JOB_AXISYM){
            pvec = body->S.operate(nodes->x_t, MPMSparseMatrixBase::TRANSPOSED);
            for (int i=0; i<points->L.size(); i++){
                points->L(i,0,2) = -pvec(i,2) / points->x(i,0);
                points->L(i,2,2) = pvec(i,0) / points->x(i,0);
            }
        }
    }
    return;
}

void FVMMixtureSolver::updateDensity(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (int i=0;i<job->bodies[b]->points->v.rows();i++) {
            job->bodies[b]->points->v(i) *= std::exp(job->dt * job->bodies[b]->points->L(i).trace());
        }

        //this is new, but maybe useful
        job->bodies[b]->points->updateIntegrators(job,job->bodies[b].get());
    }
    return;
}

void FVMMixtureSolver::updateStress(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b]->activeMaterial == 0){
            continue;
        }
        job->bodies[b]->material->calculateStress(job, job->bodies[b].get(), Material::UPDATE);
    }
    return;
}

void FVMMixtureSolver::generateLoads(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        job->bodies[b]->generateLoads(job);
    }
    return;
}

void FVMMixtureSolver::applyLoads(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        job->bodies[b]->applyLoads(job);
    }
    return;
}
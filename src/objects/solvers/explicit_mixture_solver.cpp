//
// Created by aaron on 5/16/18.
// explicit_usl.cpp
//

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "solvers.hpp"
#include "objects/bodies/bodies.hpp"
#include "objects/materials/materials.hpp"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

/*----------------------------------------------------------------------------*/
//
void ExplicitMixtureSolver::init(Job* job){

    //initialize parent class
    ExplicitUSL::init(job);

    //identify bodies
    if (str_props.size() < 3){
        std::cerr << str_props.size() << " != 3\n";
        fprintf(stderr, "%s:%s:", __FILE__, __func__);

        std::cerr << "ExplicitMixtureSolver needs at least 3 bodies defined by name.\n";
        std::cerr << "    granular_body:\n";
        std::cerr << "    fluid_body:\n";
        std::cerr << "    solid_body:\n";
        exit(0);
    } else {

        //initialization boolean
        bool initilization_fail = false;

        //set bodyIDs by name
        std::vector<int> bodyIDs = {-1, -1, -1};
        for (int i = 0; i < bodyIDs.size(); i++) {
            for (int b = 0; b < job->bodies.size(); b++) {
                if (str_props[i].compare(job->bodies[b]->name) == 0) {
                    bodyIDs[i] = b;
                    break;
                }
            }
            if (bodyIDs[i] <= 0){
                std::cerr << "Uh oh! [" << i << "]" << std::endl;
                initilization_fail = true;
                break;
            }
        }

        if (initilization_fail){
            std::cerr << "ExplicitMixtureSolver needs at least 3 bodies defined by name.\n";
            std::cerr << "    granular_body: " << str_props[0] << "\n";
            std::cerr << "    fluid_body: " << str_props[1] << "\n";
            std::cerr << "    solid_body: " << str_props[2] << "\n";
            exit(0);
        }

        //assign pointers
        granular_body = job->bodies[bodyIDs[0]].get();
        fluid_body    = job->bodies[bodyIDs[1]].get();
        solid_body    = job->bodies[bodyIDs[2]].get();

        //tell console about findings!
        std::cout << "ExplicitMixtureSolver has identified 3 bodies!\n";
        std::cout << "    granular_body: " << granular_body->name << "\n";
        std::cout << "    fluid_body: "    << fluid_body->name << "\n";
        std::cout << "    solid_body: "    << solid_body->name << "\n";

        //attempt to identify granular material model
        if (granular_body->material->object_name.compare("CompressibleBreakageMechanicsSand") == 0){
            //success!

            //attempt high risk of pointer to granular material model
            compressible_breakage_mechanics_sand_model = dynamic_cast<CompressibleBreakageMechanicsSand*>(granular_body->material.get());

            //assume success and deal with the consequences
            std::cout << "ExplicitMixtureSolver has identified granular material model!\n";
            std::cout << "   " << compressible_breakage_mechanics_sand_model->object_name << " =? CompressibleBreakageMechanicsSand\n";
        } else {
            std::cerr << "ExplicitMixtureSolver has failed to identify granular material model!\n";
            std::cout << "   " << granular_body->material->object_name << " =? CompressibleBreakageMechanicsSand\n";
        }

        //attempt to identify fluid material model
        if (fluid_body->material->object_name.compare("BarotropicViscousFluid") == 0){
            //success!
            fluid_model = 0;

            //attempt high risk of pointer to granular material model
            barotropic_viscous_fluid_model = dynamic_cast<BarotropicViscousFluid*>(fluid_body->material.get());

            //assume success and deal with the consequences
            std::cout << "ExplicitMixtureSolver has identified fluid material model!\n";
            std::cout << "   " << barotropic_viscous_fluid_model->object_name << " =? BarotropicViscousFluid\n";
        } else {
            std::cerr << "ExplicitMixtureSolver has failed to identify granular material model!\n";
            std::cout << "   " << granular_body->material->object_name << " =? CompressibleBreakageMechanicsSand\n";
        }

    }

    return;
}

/*----------------------------------------------------------------------------*/
//
void ExplicitMixtureSolver::step(Job* job){
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

    return;
}

/*----------------------------------------------------------------------------*/
//
std::string ExplicitMixtureSolver::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

/*----------------------------------------------------------------------------*/
//
int ExplicitMixtureSolver::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//
void ExplicitMixtureSolver::createMappings(Job *job){
    for (int b=0;b<job->bodies.size();b++){
        //job->bodies[b]->generateMap(job, DefaultBody::CPDI_ON); //use_cpdi by default
        job->bodies[b]->generateMap(job, cpdi_spec);
    }
    return;
}

void ExplicitMixtureSolver::mapPointsToNodes(Job* job){
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

void ExplicitMixtureSolver::generateContacts(Job* job){
    for (int c=0;c<job->contacts.size();c++){
        if (job->activeContacts[c] == 0){
            continue;
        }
        job->contacts[c]->generateRules(job);
    }
    return;
}

void ExplicitMixtureSolver::addContacts(Job* job){
    for (int c=0;c<job->contacts.size();c++){
        if (job->activeContacts[c] == 0){
            continue;
        }
        //job->contacts[c]->applyRules(job,Contact::EXPLICIT);
        job->contacts[c]->applyRules(job, contact_spec);
    }
    return;
}

void ExplicitMixtureSolver::generateBoundaryConditions(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b]->activeBoundary == 0){
            continue;
        }
        job->bodies[b]->boundary->generateRules(job,job->bodies[b].get());
    }
    return;
}

void ExplicitMixtureSolver::addBoundaryConditions(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b]->activeBoundary == 0){
            continue;
        }
        job->bodies[b]->boundary->applyRules(job,job->bodies[b].get());
    }
    return;
}

void ExplicitMixtureSolver::moveGrid(Job* job){
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

void ExplicitMixtureSolver::movePoints(Job* job){
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

void ExplicitMixtureSolver::calculateStrainRate(Job* job){
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

void ExplicitMixtureSolver::updateDensity(Job* job){
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

void ExplicitMixtureSolver::updateStress(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b]->activeMaterial == 0){
            continue;
        }
        job->bodies[b]->material->calculateStress(job, job->bodies[b].get(), Material::UPDATE);
    }
    return;
}
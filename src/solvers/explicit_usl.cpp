//
// Created by aaron on 5/25/17.
// explicit_usl.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <driver.hpp>

#include "job.hpp"
#include "serializer.hpp"

#include "solver.hpp"

#include "grid.hpp"
#include "contact.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "points.hpp"

#include "material.hpp"
#include "boundary.hpp"

extern "C" void solverInit(Job* job); //initialize solver
extern "C" void solverStep(Job* job); //one forward mpm step

extern "C" std::string solverSaveState(Job* job, Serializer* serializer, std::string filepath); //save solver state to returned filename in serializer folder
extern "C" int solverLoadState(Job* job, Serializer* serializer, std::string fullpath); //load state from given full path

/*----------------------------------------------------------------------------*/

void createMappings(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        job->bodies[b].bodyGenerateMap(job, Body::CPDI_ON); //use_cpdi by default
    }
    return;
}

void mapPointsToNodes(Job* job){
    Body *body;
    Points *points;
    Nodes *nodes;
    Eigen::Matrix pvec;
    Eigen::Matrix nvec;
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = &(job->bodies[b]);
        points = &(job->bodies[b].points);
        nodes = &(job->bodies[b].nodes);

        //map mass
        body->bodyCalcNodalValues(job,nodes->m,points->m,Body::SET);

        //map momentum
        for (size_t i=0;i<points->mx_t.cols();i++){
            points->mx_t.col(i) = points->m.array() * points->x_t.col(i).array();
        }
        body->bodyCalcNodalValues(job,nodes->mx_t,points->mx_t,Body::SET);

        //calculate velocity
        for (size_t i=0;i<nodes->x_t.rows();i++){
            if (nodes->m(i) > 0){
                nodes->x_t.row(i) = nodes->mx_t.row(i) / nodes->m(i);
            } else {
                nodes->x_t.row(i).setZero();
            }
        }

        //map body force
        pvec = job->jobVectorArray(points->b.rows());
        for (size_t i=0;i<points->b.cols();i++){
            pvec.col(i) = points->m.array() * points->b.col(i).array();
        }
        body->bodyCalcNodalValues(job,nodes->f,pvec,Body::SET);

        //map divergence of stress
        body->bodyCalcNodalDivergence(job,nodes->f,points->T,Body::ADD);
        //nvec = job->jobVectorArray(nodes->f.rows());
        //body->bodyCalcNodalDivergence(job,nvec,points->T);
        //nodes->f += nvec;
    }
    return;
}

void generateContacts(Job* job){
    for (size_t c=0;c<job->contacts.size();c++){
        if (job->activeContacts[c] == 0){
            continue;
        }
        job->contacts[c].contactGenerateRules(job);
    }
    return;
}

void addContacts(Job* job){
    for (size_t c=0;c<job->contacts.size();c++){
        if (job->activeContacts[c] == 0){
            continue;
        }
        job->contacts[c].contactApplyRules(job);
    }
    return;
}

void generateBoundaryConditions(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b].activeBoundary == 0){
            continue;
        }
        job->bodies[b].boundary.boundaryGenerateRules(job,&(job->bodies[b]));
    }
    return;
}

void addBoundaryConditions(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b].activeBoundary == 0){
            continue;
        }
        job->bodies[b].boundary.boundaryApplyRules(job,&(job->bodies[b]));
    }
    return;
}

void moveGrid(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //update momentum
        job->bodies[b].nodes.mx_t += job->dt * job->bodies[b].nodes.f;

        //calculate velocity
        for (size_t i=0;i<job->bodies[b].nodes.x_t.rows();i++){
            if (job->bodies[b].nodes.m(i) > 0) {
                job->bodies[b].nodes.x_t.row(i) = job->bodies[b].nodes.mx_t.row(i) / job->bodies[b].nodes.m(i);
            } else {
                job->bodies[b].nodes.x_t.setZero();
            }
        }

        //set displacement
        job->bodies[b].nodes.u = job->dt * job->bodies[b].nodes.x_t;

        //calculate difference in velocity
        for (size_t i=0;i<job->bodies[b].nodes.diff_x_t.rows();i++){
            if (job->bodies[b].nodes.m(i) > 0) {
                job->bodies[b].nodes.diff_x_t.row(i) = job->bodies[b].nodes.f.row(i) / job->bodies[b].nodes.m(i);
            } else {
                job->bodies[b].nodes.diff_x_t.setZero();
            }
        }
    }
    return;
}

void movePoints(Job* job){
    Body* body;
    Points* points;
    Nodes* nodes;
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = &(job->bodies[b]);
        points = &(job->bodies[b].points);
        nodes = &(job->bodies[b].nodes);

        //map nodal displacement to point positions
        body->bodyCalcPointValues(job,points->x,nodes->u,Body::ADD);
        body->bodyCalcPointValues(job,points->u,nodes->u,Body::ADD);

        //map nodal velocity diff to points
        body->bodyCalcPointValues(job,points->x_t,nodes->diff_x_t,Body::ADD);

        //calculate momentum
        for (size_t i=0;i<points->mx_t.cols();i++){
            points->mx_t.col(i) = points->m.array() * points->x_t.col(i).array();
        }
    }
    return;
}

void calculateStrainRate(Job* job){
    Body* body;
    Points* points;
    Nodes* nodes;
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = &(job->bodies[b]);
        points = &(job->bodies[b].points);
        nodes = &(job->bodies[b].nodes);

        //calculate gradient of nodal velocity at points
        body->bodyCalcPointGradient(job,points->L,nodes->x_t,Body::SET);
    }
    return;
}

void updateDensity(Job* job){
    Eigen::MatrixXd L;
    Eigen::VectorXd tmpVec;

    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (size_t i=0;i<job->bodies[b].points.v.rows();i++) {
            tmpVec << job->bodies[b].points.L.row(i).transpose();
            L = job->jobTensor(tmpVec.data());
            job->bodies[b].points.v(i) *= std::exp(job->dt * L.trace());
        }
    }
    return;
}

void updateStress(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b].activeMaterial == 0){
            continue;
        }
        job->bodies[b].material.materialCalculateStress(job,&(job->bodies[b]),Material::UPDATE);
    }
    return;
}

/*----------------------------------------------------------------------------*/

void solverInit(Job* job){
    std::cout << "Solver Initialized." << std::endl;
    return;
}

void solverStep(Job* job){
    //create map
    createMappings(job);

    //map particles to grid
    mapPointsToNodes(job);

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
    job->driver.driverGenerateGravity(job);
    job->driver.driverApplyGravity(job);

    //update stress
    updateStress(job);

    return;
}

std::string solverSaveState(Job* job, Serializer* serializer, std::string filepath){
    //nothing to save
    return "";
}
int solverLoadState(Job* job, Serializer* serializer, std::string fullpath){
    //nothing to load
    return 1;
}

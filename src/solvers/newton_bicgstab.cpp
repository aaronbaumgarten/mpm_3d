//
// Created by aaron on 6/13/17.
// newton_bicgstab.cpp
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

double TOL = 1e-14;
double h = 1e-5;
int maxIter = 100;

std::vector<Body> tmp_bodies(0);
Eigen::VectorXd mv_k(0,1); //initial momentum state
Eigen::VectorXd v_k(0,1); //initial velocity state
Eigen::VectorXd mv_n(0,1); //guess momentum state
Eigen::VectorXd v_n(0,1); //guess velocity state
Eigen::VectorXd f_n(0,1); //guess force state
Eigen::VectorXd m(0,1); //mass state
Eigen::VectorXd x(0,1); //search direction

Eigen::VectorXd Fv(0,1); //residual of velocity guess
Eigen::VectorXd Fv_plus(0,1); //residual of velocity plus search step

Eigen::VectorXd p(0,1);
Eigen::VectorXd r(0,1);
Eigen::VectorXd rhat(0,1);
Eigen::VectorXd v(0,1);
Eigen::VectorXd s(0,1);
Eigen::VectorXd t(0,1);

Eigen::VectorXd tmpVec(0,1);
double alpha = 0;
double beta = 0;
double rho = 0;
double omega = 0;
double rr = 0;

double rnorm = 0;
double bnorm = 0;

Eigen::VectorXd x_min(0,1); //best search direction
Eigen::VectorXd mv_n_trial(0,1); //trial guess
double mv_n_res = 1;
double rnorm_min = 1;
double lambda = 1;

extern "C" void solverInit(Job* job); //initialize solver
extern "C" void solverStep(Job* job); //one forward mpm step

extern "C" std::string solverSaveState(Job* job, Serializer* serializer, std::string filepath); //save solver state to returned filename in serializer folder
extern "C" int solverLoadState(Job* job, Serializer* serializer, std::string fullpath); //load state from given full path

/*----------------------------------------------------------------------------*/

void solverInit(Job* job){
    if (job->solver.fp64_props.size() < 2 && job->solver.int_props.size() > 1){
        std::cout << job->solver.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 3 properties defined ({TOL, linear_step}, {maxIter}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        TOL = job->solver.fp64_props[0];
        h = job->solver.fp64_props[1];
        maxIter = job->solver.int_props[0];
        printf("Solver properties (TOL = %g, linear_step = %g, maxIter = %i).\n",
               TOL, h, maxIter);
    }

    //create tmp bodies vector
    tmp_bodies.reserve(job->bodies.size());
    for (size_t b=0;b<job->bodies.size();b++){
        tmp_bodies.push_back(Body());
    }

    std::cout << "Solver Initialized." << std::endl;
    return;
}

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
    Eigen::MatrixXd pvec;
    Eigen::MatrixXd nvec;
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
        pvec = job->jobVectorArray<double>(points->b.rows());
        for (size_t i=0;i<points->b.cols();i++){
            pvec.col(i) = points->m.array() * points->b.col(i).array();
        }
        body->bodyCalcNodalValues(job,nodes->f,pvec,Body::SET);

        //map divergence of stress
        body->bodyCalcNodalDivergence(job,nodes->f,points->T,Body::ADD);
        //nvec = job->jobVectorArray<double>(nodes->f.rows());
        //body->bodyCalcNodalDivergence(job,nvec,points->T);
        //nodes->f += nvec;
    }
    return;
}

void mapPointForcesToNodes(Job* job){
    Body *body;
    Points *points;
    Nodes *nodes;
    Eigen::MatrixXd pvec;
    Eigen::MatrixXd nvec;
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = &(job->bodies[b]);
        points = &(job->bodies[b].points);
        nodes = &(job->bodies[b].nodes);

        //map body force
        pvec = job->jobVectorArray<double>(points->b.rows());
        for (size_t i=0;i<points->b.cols();i++){
            pvec.col(i) = points->m.array() * points->b.col(i).array();
        }
        body->bodyCalcNodalValues(job,nodes->f,pvec,Body::SET);

        //map divergence of stress
        body->bodyCalcNodalDivergence(job,nodes->f,points->T,Body::ADD);
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

void addContacts(Job* job, int SPEC){
    for (size_t c=0;c<job->contacts.size();c++){
        if (job->activeContacts[c] == 0){
            continue;
        }
        job->contacts[c].contactApplyRules(job, SPEC);
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
                job->bodies[b].nodes.x_t.row(i).setZero();
            }
        }

        //set displacement
        //job->bodies[b].nodes.u = job->dt * (v_k + job->bodies[b].nodes.x_t)/2.0;
        size_t body_start_pos, rel_start_pos;
        body_start_pos = b * job->DIM * job->grid.node_count;
        for (size_t pos = 0; pos < job->DIM; pos++) {
            rel_start_pos = pos * job->grid.node_count;
            job->bodies[b].nodes.u.col(pos) = job->dt * (v_k.block(body_start_pos + rel_start_pos, 0, job->grid.node_count, 1) + job->bodies[b].nodes.x_t.col(pos))/2.0;
        }

        //calculate difference in velocity
        for (size_t i=0;i<job->bodies[b].nodes.diff_x_t.rows();i++){
            if (job->bodies[b].nodes.m(i) > 0) {
                job->bodies[b].nodes.diff_x_t.row(i) = job->dt * job->bodies[b].nodes.f.row(i) / job->bodies[b].nodes.m(i);
            } else {
                job->bodies[b].nodes.diff_x_t.row(i).setZero();
            }
        }
    }
    return;
}

void displaceGrid(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //set displacement
        //job->bodies[b].nodes.u = job->dt * (v_k + job->bodies[b].nodes.x_t)/2.0;
        size_t body_start_pos, rel_start_pos;
        body_start_pos = b * job->DIM * job->grid.node_count;
        for (size_t pos = 0; pos < job->DIM; pos++) {
            rel_start_pos = pos * job->grid.node_count;
            job->bodies[b].nodes.u.col(pos) = job->dt * (v_k.block(body_start_pos + rel_start_pos, 0, job->grid.node_count, 1) + job->bodies[b].nodes.x_t.col(pos))/2.0;
        }

        //calculate difference in velocity
        for (size_t i=0;i<job->bodies[b].nodes.diff_x_t.rows();i++){
            if (job->bodies[b].nodes.m(i) > 0) {
                job->bodies[b].nodes.diff_x_t.row(i) = job->dt * job->bodies[b].nodes.f.row(i) / job->bodies[b].nodes.m(i);
            } else {
                job->bodies[b].nodes.diff_x_t.row(i).setZero();
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
            tmpVec = job->bodies[b].points.L.row(i).transpose();
            L = job->jobTensor<double>(tmpVec.data());
            job->bodies[b].points.v(i) *= std::exp(job->dt * L.trace());
        }
    }
    return;
}

void updateStress(Job* job, int SPEC){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b].activeMaterial == 0){
            continue;
        }
        job->bodies[b].material.materialCalculateStress(job,&(job->bodies[b]), SPEC);
    }
    return;
}

/*----------------------------------------------------------------------------*/

void sizeVectors(Job* job){
    size_t len = job->bodies.size() * job->DIM * job->grid.node_count;
    mv_k.resize(len);
    v_k.resize(len);
    mv_n.resize(len);
    v_n.resize(len);
    f_n.resize(len);
    m.resize(len);
    x.resize(len);

    Fv.resize(len);
    Fv_plus.resize(len);

    p.resize(len);
    r.resize(len);
    rhat.resize(len);
    v.resize(len);
    s.resize(len);
    t.resize(len);

    tmpVec.resize(len);

    x_min.resize(len);
    mv_n_trial.resize(len);
    return;
}

void fillMomentumVector(Job* job, Eigen::VectorXd& mvOut){
    if (mvOut.rows() < job->bodies.size() * job->DIM * job->grid.node_count){
        std::cerr << "ERROR: fillMomentumVector() passed vector of size: " << mvOut.rows() << ". Need: " << job->bodies.size() * job->DIM * job->grid.node_count << "." << std::endl;
        return;
    }

    size_t body_start_pos, rel_start_pos;
    for (size_t b=0;b<job->bodies.size();b++){
        body_start_pos = b * job->DIM * job->grid.node_count;
        if(job->activeBodies[b] == 0){
            mvOut.block(body_start_pos, 0, job->DIM * job->grid.node_count, 1).setZero();
        } else {
            for (size_t pos = 0; pos < job->DIM; pos++) {
                rel_start_pos = pos * job->grid.node_count;
                mvOut.block(body_start_pos + rel_start_pos, 0, job->grid.node_count, 1) = job->bodies[b].nodes.mx_t.col(
                        pos);
            }
        }
    }

    return;
}

void fillForceVector(Job* job, Eigen::VectorXd& fOut){
    if (fOut.rows() < job->bodies.size() * job->DIM * job->grid.node_count){
        std::cerr << "ERROR: fillForceVector() passed vector of size: " << fOut.rows() << ". Need: " << job->bodies.size() * job->DIM * job->grid.node_count << "." << std::endl;
        return;
    }

    size_t body_start_pos, rel_start_pos;
    for (size_t b=0;b<job->bodies.size();b++){
        body_start_pos = b * job->DIM * job->grid.node_count;
        if(job->activeBodies[b] == 0){
            fOut.block(body_start_pos, 0, job->DIM * job->grid.node_count, 1).setZero();
        } else {
            for (size_t pos = 0; pos < job->DIM; pos++) {
                rel_start_pos = pos * job->grid.node_count;
                fOut.block(body_start_pos + rel_start_pos, 0, job->grid.node_count, 1) = job->bodies[b].nodes.f.col(
                        pos);
            }
        }
    }

    return;
}

void fillMassVector(Job* job, Eigen::VectorXd& mOut){
    if (mOut.rows() < job->bodies.size() * job->DIM * job->grid.node_count){
        std::cerr << "ERROR: fillMassVector() passed vector of size: " << mOut.rows() << ". Need: " << job->bodies.size() * job->DIM * job->grid.node_count << "." << std::endl;
        return;
    }

    size_t body_start_pos, rel_start_pos;
    for (size_t b=0;b<job->bodies.size();b++){
        body_start_pos = b * job->DIM * job->grid.node_count;
        if(job->activeBodies[b] == 0){
            mOut.block(body_start_pos, 0, job->DIM * job->grid.node_count, 1).setZero();
        } else {
            for (size_t pos = 0; pos < job->DIM; pos++) {
                rel_start_pos = pos * job->grid.node_count;
                mOut.block(body_start_pos + rel_start_pos, 0, job->grid.node_count, 1) = job->bodies[b].nodes.m;
            }
        }
    }

    return;
}

void pushMomentumVector(Job* job, Eigen::VectorXd& mvIn){
    if (mvIn.rows() < job->bodies.size() * job->DIM * job->grid.node_count){
        std::cerr << "ERROR: pushMomentumVector() passed vector of size: " << mvIn.rows() << ". Need: " << job->bodies.size() * job->DIM * job->grid.node_count << "." << std::endl;
        return;
    }

    size_t body_start_pos, rel_start_pos;
    for (size_t b=0;b<job->bodies.size();b++){
        body_start_pos = b * job->DIM * job->grid.node_count;
        if(job->activeBodies[b] == 0){
            continue;
        } else {
            for (size_t pos = 0; pos < job->DIM; pos++) {
                rel_start_pos = pos * job->grid.node_count;
                job->bodies[b].nodes.mx_t.col(pos) = mvIn.block(body_start_pos + rel_start_pos, 0, job->grid.node_count, 1);
                for (size_t i = 0; i < job->grid.node_count; i++) {
                    if (job->bodies[b].nodes.m(i) > 0) {
                        job->bodies[b].nodes.x_t(i,pos) = job->bodies[b].nodes.mx_t(i, pos) / job->bodies[b].nodes.m(i);
                    } else {
                        job->bodies[b].nodes.x_t(i,pos) = 0;
                    }
                }
            }
        }
    }
}

void pushForceVector(Job* job, Eigen::VectorXd& fIn){
    if (fIn.rows() < job->bodies.size() * job->DIM * job->grid.node_count){
        std::cerr << "ERROR: pushForceVector() passed vector of size: " << fIn.rows() << ". Need: " << job->bodies.size() * job->DIM * job->grid.node_count << "." << std::endl;
        return;
    }

    size_t body_start_pos, rel_start_pos;
    for (size_t b=0;b<job->bodies.size();b++){
        body_start_pos = b * job->DIM * job->grid.node_count;
        if(job->activeBodies[b] == 0){
            continue;
        } else {
            for (size_t pos = 0; pos < job->DIM; pos++) {
                rel_start_pos = pos * job->grid.node_count;
                job->bodies[b].nodes.f.col(pos) = fIn.block(body_start_pos + rel_start_pos, 0, job->grid.node_count, 1);
            }
        }
    }
}

void storeBodies(Job* job){
    for (size_t b=0;b<tmp_bodies.size();b++){
        //points
        tmp_bodies[b].points.v = job->bodies[b].points.v;
        tmp_bodies[b].points.T = job->bodies[b].points.T;
        //tmp_bodies[b].points.L = job->bodies[b].points.L;

        //nodes
        tmp_bodies[b].nodes.mx_t = job->bodies[b].nodes.mx_t;
        tmp_bodies[b].nodes.x_t = job->bodies[b].nodes.x_t;
        tmp_bodies[b].nodes.f = job->bodies[b].nodes.f;
    }
    return;
}

void revertBodies(Job* job){
    //CAREFUL WITH THIS!
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        //points
        job->bodies[b].points.v = tmp_bodies[b].points.v;
        job->bodies[b].points.T = tmp_bodies[b].points.T;
        //job->bodies[b].points.L = tmp_bodies[b].points.L;

        //nodes
        job->bodies[b].nodes.mx_t = tmp_bodies[b].nodes.mx_t;
        job->bodies[b].nodes.x_t = tmp_bodies[b].nodes.x_t;
        job->bodies[b].nodes.f = tmp_bodies[b].nodes.f;
    }
    return;
}

/*----------------------------------------------------------------------------*/

void trialStep(Job* job, Eigen::VectorXd& forceOut, Eigen::VectorXd& momentumIn){
    //revert bodies
    revertBodies(job);

    //push momentum to nodes --------------------------
    pushMomentumVector(job, momentumIn);
    addBoundaryConditions(job);

    //calculate resulting force -----------------------
    calculateStrainRate(job);

    updateDensity(job);

    job->driver.driverApplyGravity(job);

    updateStress(job, Material::TRIAL);

    mapPointForcesToNodes(job);

    addContacts(job, Contact::IMPLICIT);

    addBoundaryConditions(job);

    fillForceVector(job, forceOut);
    return;
}

void applyA(Job* job, Eigen::VectorXd& Js_vec, Eigen::VectorXd& s_vec){
    //calculate v = Ap --------------------------------
    if (v_n.norm() > 0) {
        tmpVec = v_n + h * s_vec * v_n.norm() / s_vec.norm(); //velocity
    } else {
        tmpVec = h * s_vec / s_vec.norm(); //velocity
    }

    tmpVec = tmpVec.array() * m.array(); //momentum

    trialStep(job,f_n,tmpVec);

    Fv_plus = mv_k + job->dt * f_n - tmpVec; //momentum

    //approximate Ap  ---------------------------------
    if (v_n.norm() > 0) {
        Js_vec = s_vec.norm() / (h * v_n.norm()) * (Fv_plus - Fv); //momentum
    } else {
        Js_vec = s_vec.norm() / h * (Fv_plus - Fv);
    }

    //adjust v = Ap -----------------------------------
    for (size_t i = 0; i < v.rows(); i++) {
        if (m(i) > 0) {
            Js_vec(i) /= m(i); //velocity
        } else {
            Js_vec(i) = 0;
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void solverStep(Job* job){
    /*------------------------------------------------------------------------*/
    //find initial state of nodes

    //create map
    createMappings(job);

    //map particles to grid
    mapPointsToNodes(job);

    generateBoundaryConditions(job); //only call once
    generateContacts(job); //only call once
    job->driver.driverGenerateGravity(job); //only call once

    //enforce boundary conditions (before any forces are calculated)
    addBoundaryConditions(job);

    //create temporary bodies (to store initial state of mv, T, v)
    tmp_bodies.reserve(job->bodies.size());
    storeBodies(job);

    //create vectors (mv_n, f_n, etc)
    sizeVectors(job);

    //assign initial state
    fillMomentumVector(job, mv_k);
    fillMassVector(job, m);
    for (size_t i=0;i<mv_k.rows();i++){
        if (m(i) > 0){
            v_k(i) = mv_k(i) / m(i);
        } else {
            v_k(i) = 0;
        }
    }

    /*------------------------------------------------------------------------*/
    //take trial step for initial guess

    //add contact forces
    addContacts(job, Contact::EXPLICIT);

    //enforce boundary conditions
    addBoundaryConditions(job);

    //move grid
    moveGrid(job);

    //store nodal velocity
    fillMomentumVector(job, mv_n);


    /*------------------------------------------------------------------------*/
    //find initial residual

    //trial step
    trialStep(job,f_n,mv_n);

    Fv = mv_k + job->dt*f_n - mv_n;
    mv_n_res = Fv.norm(); //momentum residual

    for (size_t i=0;i<mv_n.rows();i++){
        if (m(i) > 0){
            v_n(i) = mv_n(i)/m(i); //velocity state
        } else {
            v_n(i) = 0;
        }
    }

    /*------------------------------------------------------------------------*/
    //correct final state

    size_t n = 0;
    size_t k = 0;
    size_t j = 0;
    std::cout << std::endl; //move cursor down one line
    while (mv_n_res / mv_k.norm() > TOL && n < maxIter) {
        //printf("\33[2K"); //clear line
        std::cout << "Newton Step: [" << n << " <? " << maxIter << "]. Residual: [" << mv_n_res / mv_k.norm() << " <? " << TOL << "]." << std::endl;
        //set residual vector
        for (size_t i = 0; i < m.rows(); i++) {
            if (m(i) > 0) {
                r(i) = -Fv(i) / m(i); //velocity
            } else {
                r(i) = 0;
            }
        }

        x.setZero();
        rhat = r; //velocity
        rnorm = r.norm();
        bnorm = rnorm;
        rnorm_min = rnorm; //save minimum residual
        x_min = x; //save best search direction

        rho = 1; alpha = 1; omega = 1;
        v.setZero(); p.setZero();

        //cg step *************************************************************
        k = 0;
        while (rnorm / bnorm > TOL && k < maxIter) {
            //printf("\33[2K"); //clear line
            std::cout << "BICGSTAB Step: [" << k << " <? " << maxIter << "]. Residual: [" << rnorm / bnorm << " <? " << TOL << "].\r";

            //bicgstab step------------------------------------
            k += 1;
            rr = rhat.dot(r);
            beta = (rr/rho) * (alpha/omega);
            rho = rr;

            p = r + beta*(p - omega*v);

            //v = Ap
            applyA(job,v,p);

            alpha = rho / rhat.dot(v);
            s = r - alpha*v;

            //t = As
            applyA(job,t,s);

            omega = t.dot(s) / t.dot(t);
            x = x + alpha*p + omega*s;
            r = s - omega*t;
            rnorm = r.norm();

            //check for new minumum ---------------------------
            if (rnorm <= rnorm_min) {
                x_min = x;
            }
        }
        //printf("\33[2K"); //clear line
        std::cout << "BICGSTAB Step: [" << k << " <? " << maxIter << "]. Residual: [" << rnorm / bnorm << " <? " << TOL << "].\r";
        printf("\033[A"); //move cursor up one line

        //newton step *********************************************************
        lambda = 1.0;
        j = 0;
        do {
            //update momentum guess------------------------
            mv_n_trial = mv_n.array() + lambda * x_min.array() * m.array(); //momentum

            trialStep(job, f_n, mv_n_trial);

            Fv = mv_k + job->dt * f_n - mv_n_trial; //momentum
            lambda /= 2.0;
            j += 1;

            //check for minimun----------------------------
        } while (Fv.norm() > mv_n_res);

        mv_n_res = Fv.norm();
        mv_n = mv_n_trial;

        for (size_t i=0;i<mv_n.rows();i++){
            if (m(i) > 0){
                v_n(i) = mv_n(i)/m(i); //velocity state
            } else {
                v_n(i) = 0;
            }
        }

        n += 1;
    }
    //printf("\33[2K"); //clear line
    std::cout << "Newton Step: [" << n << " <? " << maxIter << "]. Residual: [" << mv_n_res / mv_k.norm() << " <? " << TOL << "].\r";
    printf("\033[A"); //move cursor up one line

    /*------------------------------------------------------------------------*/
    //propogate final state

    //final revert
    revertBodies(job);

    //pushMomentumVector(job,mv_n);
    pushMomentumVector(job,mv_n);
    pushForceVector(job,f_n);

    addBoundaryConditions(job);

    displaceGrid(job);

    movePoints(job);

    calculateStrainRate(job);

    updateDensity(job);

    updateStress(job, Material::UPDATE);

    job->driver.driverApplyGravity(job);

    return;
}

std::string solverSaveState(Job* job, Serializer* serializer, std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename/directory name
    std::ostringstream s;
    s << "mpm_v2.solver." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    //create filename
    filename = s.str();

    std::ofstream ffile((filepath+filename), std::ios::trunc);

    //write data
    if (ffile.is_open()) {
        size_t len = tmp_bodies.size();

        ffile << "# mpm_v2 Solver\n";
        ffile << len << "\n"; //length of body vector
        ffile << TOL << "\n";
        ffile << h << "\n";
        ffile << maxIter << "\n";
        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filename << "\" !\n";
        return "ERR";
    }

    //assume that original bodies were mangled. revert bodies
    revertBodies(job);

    return filename;
}
int solverLoadState(Job* job, Serializer* serializer, std::string fullpath){
    std::string line; //read line

    std::ifstream fin(fullpath); //file to load from

    if (fin.is_open()) {
        //if open, read lines
        std::getline(fin,line); //header
        std::getline(fin,line); //len
        tmp_bodies.reserve(std::stoi(line)); //reserve length of body vector
        for (size_t b=0;b<std::stoi(line); b++){
            tmp_bodies.push_back(Body());
        }
        std::getline(fin,line); //TOL
        TOL = std::stod(line);
        std::getline(fin,line); //h
        h = std::stod(line);
        std::getline(fin,line); //maxIter
        maxIter = std::stoi(line);
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Solver Loaded." << std::endl;
    return 1;
}

//
// Created by aaron on 10/10/22.
//

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

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

/*----------------------------------------------------------------------------*/
//
void GeneralizedVortexErrorSolver::init(Job* job){
    //must be 2D
    if (job->JOB_TYPE != 2){
        std::cerr << "ERROR! GeneralizedVortexErrorSolver is only designed for 2D problems! Exiting." << std::endl;
        exit(0);
    }

    if (fp64_props.size() < 3 || str_props.size() < 1){
        std::cerr << "Need 4 properties defined: {A, B, density} and {filename}." << std::endl;
        exit(0);
    } else {
        A = fp64_props[0];
        B = fp64_props[1];
        density = fp64_props[2];
        output_filename = str_props[0];
    }

    //open and clear file
    std::ofstream file (output_filename,std::ios::trunc);
    if (file.is_open()){
        //success!
        //write file header
        file << "t, ||u^* - u^Q||_L2, ||v^* - v^Q||_L2, ||a^* - a^Q||_L2\n";

        file.close();
    } else {
        std::cerr << "ERROR! Cannot open " << output_filename << "! Exiting." << std::endl;
        exit(0);
    }

    std::cout << "Solver properties (A = " << A << ", B = " << B << ", filename = " << output_filename << ")." << std::endl;
    std::cout << "Solver Initialized." << std::endl;
    return;
}

double GeneralizedVortexErrorSolver::alpha(double r, double t){
    double d = r - 1.0;
    return A * std::sin(B * M_PI * t) * (1.0 - 32.0 * d * d + 256 * d * d * d * d);
}

double GeneralizedVortexErrorSolver::dalpha_dr(double r, double t){
    double d = r - 1.0;
    return A * std::sin(B * M_PI * t) * (-64.0 * d + 1024 * d * d * d);
}

double GeneralizedVortexErrorSolver::d2alpha_dr2(double r, double t) {
    double d = r - 1.0;
    return A * std::sin(B * M_PI * t) * (-64.0 + 3072 * d * d);
}

double GeneralizedVortexErrorSolver::dalpha_dt(double r, double t){
    double d = r - 1.0;
    return A * B * M_PI * std::cos(B * M_PI * t) * (1.0 - 32.0 * d * d + 256 * d * d * d * d);
}

double GeneralizedVortexErrorSolver::d2alpha_dt2(double r, double t){
    double d = r - 1.0;
    return -A * B * B * M_PI * M_PI * std::sin(B * M_PI * t) * (1.0 - 32.0 * d * d + 256 * d * d * d * d);
}

KinematicVector GeneralizedVortexErrorSolver::getDisplacement(Job* job, KinematicVector const &x){
    //radial position
    double r = std::sqrt((x[0] - x0)*(x[0] - x0) + (x[1] - y0)*(x[1] - y0));

    //theta position
    double t = std::atan((x[1] - y0)/(x[0] - x0));
    if (std::abs(x[0] - x0) < 1e-10){
        if (x[1] - y0 >= 0){
            t = M_PI / 2.0;
        } else {
            t = -M_PI / 2.0;
        }
    } else if (x[0] - x0 < 0){
        t += M_PI;
    }

    //radial displacement
    double dr = 0.0;

    //theta displacement
    double dt = alpha(r, job->t);

    KinematicVector result = KinematicVector(job->JOB_TYPE);
    if (r > a && r < b) {
        result[0] = x[0] - x0 - r * std::cos(t - dt);
        result[1] = x[1] - y0 - r * std::sin(t - dt);
    } else {
        result.setZero();
    }
    return result;
}

KinematicVector GeneralizedVortexErrorSolver::getVelocity(Job* job, KinematicVector const &x){
    //radial position
    double r = std::sqrt((x[0] - x0)*(x[0] - x0) + (x[1] - y0)*(x[1] - y0));

    //radial displacement
    double vr = 0.0;

    //theta displacement
    double vt = r * dalpha_dt(r, job->t);

    KinematicVector result = KinematicVector(job->JOB_TYPE);
    if (r > a && r < b) {
        result[0] = vr * (x[0] - x0) / r - vt * (x[1] - y0) / r;
        result[1] = vr * (x[1] - y0) / r + vt * (x[0] - x0) / r;
    } else {
        result.setZero();
    }
    return result;
}

KinematicVector GeneralizedVortexErrorSolver::getAcceleration(Job* job, KinematicVector const &x){
    //radial position
    double r = std::sqrt((x[0] - x0)*(x[0] - x0) + (x[1] - y0)*(x[1] - y0));

    //radial displacement
    double ar = -r * dalpha_dt(r, job->t) * dalpha_dt(r, job->t);

    //theta displacement
    double at = r * d2alpha_dt2(r, job->t);

    KinematicVector result = KinematicVector(job->JOB_TYPE);
    if (r > a && r < b) {
        result[0] = ar * (x[0] - x0) / r - at * (x[1] - y0) / r;
        result[1] = ar * (x[1] - y0) / r + at * (x[0] - x0) / r;
    } else {
        result.setZero();
    }
    return result;
}

KinematicVector GeneralizedVortexErrorSolver::getBodyForce(Job* job, KinematicVector const &x){
    //radial position
    double r = std::sqrt((x[0] - x0)*(x[0] - x0) + (x[1] - y0)*(x[1] - y0));

    //radial displacement
    double br = -density * r * dalpha_dt(r, job->t) * dalpha_dt(r, job->t) +
                G * r * dalpha_dr(r, job->t) * dalpha_dr(r, job->t);

    //theta displacement
    double bt = density * r * d2alpha_dt2(r, job->t) - G * r * d2alpha_dr2(r, job->t) -
                3.0 * G * dalpha_dr(r, job->t);

    KinematicVector result = KinematicVector(job->JOB_TYPE);
    if (r > a && r < b) {
        result[0] = br * (x[0] - x0) / r - bt * (x[1] - y0) / r;
        result[1] = br * (x[1] - y0) / r + bt * (x[0] - x0) / r;
    } else {
        result.setZero();
    }
    return result;
}

MaterialTensor GeneralizedVortexErrorSolver::getStress(Job* job, KinematicVector const &x){
    //radial position
    double r = std::sqrt((x[0] - x0)*(x[0] - x0) + (x[1] - y0)*(x[1] - y0));

    //radial stress
    double srr = 0.0;

    //shear stress
    double srt = G * r * dalpha_dr(r, job->t);

    //theta stress
    double stt = G * r * r * dalpha_dr(r, job->t) * dalpha_dr(r, job->t);

    MaterialTensor result = MaterialTensor();
    MaterialTensor Q = MaterialTensor();
    if (r > a && r < b) {
        result(0,0) = srr;
        result(0, 1) = srt;
        result(1, 0) = srt;
        result(1, 1) = stt;

        Q(0, 0) = (x[0] - x0) / r;     Q(0, 1) = (x[1] - y0) / r;
        Q(1, 0) = -(x[1] - y0) / r;    Q(1, 1) = (x[0] - x0) / r;

        //transform from rotated coordinate system back to x,y
        result = Q * result * Q.transpose();
    } else {
        result.setZero();
    }
    return result;
}

void GeneralizedVortexErrorSolver::setGridLengths(Job* job){
    //get grid dimensions
    KinematicVector x_min = job->grid->nodeIDToPosition(job, 0);
    KinematicVector x_max = x_min;
    KinematicVector tmp_x;
    for (int i=1; i<job->grid->node_count; i++){
        tmp_x = job->grid->nodeIDToPosition(job, i);
        for (int pos=0; pos<job->grid->GRID_DIM; pos++){
            if (tmp_x[pos] > x_max[pos]){
                x_max[pos] = tmp_x[pos];
            } else if (tmp_x[pos] < x_min[pos]){
                x_min[pos] = tmp_x[pos];
            }
        }
    }
    Lx = x_max - x_min;
    return;
}

void GeneralizedVortexErrorSolver::setInitialVelocity(Job* job){

    //assign momentum and velocity to points
    KinematicVector tmpX;
    for (int p = 0; p<job->bodies[0]->points->x.size(); p++){
        tmpX = job->bodies[0]->points->x[p];
        job->bodies[0]->points->x_t[p] = getVelocity(job, tmpX);
        job->bodies[0]->points->mx_t[p] = job->bodies[0]->points->m(p) *
                job->bodies[0]->points->x_t[p];
    }

    return;
}

/*----------------------------------------------------------------------------*/
//
void GeneralizedVortexErrorSolver::step(Job* job){
    //get grid dimensions on first step
    if (job->t < job->dt){
        setGridLengths(job);
        setInitialVelocity(job);
    }

    //create map
    createMappings(job);

    //map points to grid
    mapPointsToNodes(job);

    //add boundary conditions
    addBoundaryConditions(job);

    //update grid quantities
    moveGrid(job);

    //calculate error and write to file
    writeErrorInfo(job);

    //map grid quantities to points
    movePoints(job);

    //calculate deformation and stress
    calculateStrainRate(job);
    updateDensity(job);
    updateStress(job);

    //update body force
    addLoads(job);

    return;
}

/*----------------------------------------------------------------------------*/
//
std::string GeneralizedVortexErrorSolver::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

/*----------------------------------------------------------------------------*/
//
int GeneralizedVortexErrorSolver::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//
void GeneralizedVortexErrorSolver::createMappings(Job *job){
    for (int b=0;b<job->bodies.size();b++){
        //job->bodies[b]->generateMap(job, DefaultBody::CPDI_ON); //use_cpdi by default
        job->bodies[b]->generateMap(job, 0);
    }
    return;
}

void GeneralizedVortexErrorSolver::mapPointsToNodes(Job* job){
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
        //std::cout << points->m.sum() << " ?= " << nodes->m.sum() << std::endl;
    }
    return;
}

void GeneralizedVortexErrorSolver::addBoundaryConditions(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b]->activeBoundary == 0){
            continue;
        }
        job->bodies[b]->boundary->applyRules(job,job->bodies[b].get());
    }
    return;
}

void GeneralizedVortexErrorSolver::moveGrid(Job* job){
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

void GeneralizedVortexErrorSolver::movePoints(Job* job){
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

void GeneralizedVortexErrorSolver::calculateStrainRate(Job* job){
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
    }
    return;
}

void GeneralizedVortexErrorSolver::updateDensity(Job* job){
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

void GeneralizedVortexErrorSolver::updateStress(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b]->activeMaterial == 0){
            continue;
        }
        job->bodies[b]->material->calculateStress(job, job->bodies[b].get(), Material::UPDATE);
    }
    return;
}

void GeneralizedVortexErrorSolver::addLoads(Job* job){
    KinematicVector tmpX;
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //add body force to points
        for (int p=0; p<job->bodies[b]->points->x.size(); p++){
            if (job->bodies[b]->points->active(p)) {
                tmpX = job->bodies[b]->points->x[p];
                job->bodies[b]->points->b[p] = job->bodies[b]->points->v(p) *
                        getBodyForce(job, tmpX) / job->bodies[b]->points->m(p);
            }
        }
    }
    return;
}


void GeneralizedVortexErrorSolver::writeErrorInfo(Job* job){
    //calculate and write error measures to output file for first body in simulation
    //BE CAREFUL!!!

    //initialize arrays
    Eigen::VectorXd V_i = Eigen::VectorXd(job->grid->node_count);
    Eigen::VectorXd v_i = Eigen::VectorXd(job->grid->node_count);

    //initialize figures of merit\
    //t, ||d^* - d^Q||_L2, ||v^* - v^Q||_L2, ||a^* - a^Q||_L2
    double u_L2 = 0;
    double v_L2 = 0;
    double a_L2 = 0;

    //get exact node volumes form grid
    for (int i=0; i<job->grid->node_count;i++){
        V_i(i) = job->grid->nodeVolume(job,i);
    }

    //get integrated node volume from points
    v_i = job->bodies[0]->S * job->bodies[0]->points->v;

    //calculate error measures
    KinematicVector tmpDis = KinematicVector(job->JOB_TYPE);
    KinematicVector tmpVel = KinematicVector(job->JOB_TYPE);
    KinematicVector tmpAcc = KinematicVector(job->JOB_TYPE);
    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    for (int p=0; p<job->bodies[0]->points->x.size(); p++){
        if (job->bodies[0]->points->active(p)){
            //||d^* - d^Q||_L2
            tmpDis = getDisplacement(job, job->bodies[0]->points->x[p]);
            tmpVec = (tmpDis - job->bodies[0]->points->u[p]);
            u_L2 += tmpVec.dot(tmpVec) * job->bodies[0]->points->v(p);
        }
    }
    for (int i=0; i<V_i.rows(); i++){
        if (job->bodies[0]->nodes->m(i) > 0) {
            //||v^* - v^Q||_L2
            tmpVel = getVelocity(job, job->bodies[0]->nodes->x[i]);
            tmpVec = (tmpVel - job->bodies[0]->nodes->mx_t[i] / job->bodies[0]->nodes->m(i));
            v_L2 += tmpVec.dot(tmpVec) * V_i(i);

            //||a^* - a^Q||_L2
            tmpAcc = getAcceleration(job, job->bodies[0]->nodes->x[i]);
            tmpVec = (tmpAcc - job->bodies[0]->nodes->f[i] / job->bodies[0]->nodes->m(i));
            a_L2 += tmpVec.dot(tmpVec) * V_i(i);
        }
    }
    //sqrt of ||a_err||_L2^2
    u_L2 = std::sqrt(u_L2);
    v_L2 = std::sqrt(v_L2);
    a_L2 = std::sqrt(a_L2);

    //open and write to file
    std::ofstream file (output_filename,std::ios::app);
    if (file.is_open()){
        //success!
        //write to file
        file << job->t << ", ";
        file << u_L2 << ", ";
        file << v_L2 << ", ";
        file << a_L2 << "\n";

        file.close();
    } else {
        std::cerr << "ERROR! Cannot open " << output_filename << "!" << std::endl;;
    }

}

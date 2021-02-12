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
void TGVErrorSolver::init(Job* job){
    //must be 2D
    if (job->JOB_TYPE != 2){
        std::cerr << "ERROR! TGVErrorSolver is only designed for 2D problems! Exiting." << std::endl;
        exit(0);
    }

    if (fp64_props.size() < 3 || str_props.size() < 1){
        std::cerr << "Need 4 properties defined: {a, u_max, density} and {filename}." << std::endl;
        exit(0);
    } else {
        a = 2.0*M_PI*fp64_props[0];
        u_max = fp64_props[1];
        density = fp64_props[2];
        output_filename = str_props[0];
    }

    //open and clear file
    std::ofstream file (output_filename,std::ios::trunc);
    if (file.is_open()){
        //success!
        //write file header
        file << "t, ||H||_2^2, ||e||_\\infty, ||a^* - a^Q||_L2, ||g||_\\infty, ||G_rel||_\\infty \n";

        file.close();
    } else {
        std::cerr << "ERROR! Cannot open " << output_filename << "! Exiting." << std::endl;
        exit(0);
    }

    std::cout << "Solver properties (a = " << a << ", u_max = " << u_max << ", filename = " << output_filename << ")." << std::endl;
    std::cout << "Solver Initialized." << std::endl;
    return;
}

KinematicVector TGVErrorSolver::getVelocity(Job* job, KinematicVector const &x){
    KinematicVector result = KinematicVector(job->JOB_TYPE);
    result[0] = u_max * std::cos(a*x[0]) * std::sin(a*x[1]);
    result[1] = -u_max * std::sin(a*x[0]) * std::cos(a*x[1]);
    return result;
}

KinematicVector TGVErrorSolver::getAcceleration(Job* job, KinematicVector const &x){
    KinematicVector result = KinematicVector(job->JOB_TYPE);
    result[0] = -u_max * u_max * a * std::sin(a*x[0]) * std::cos(a*x[0]);
    result[1] = -u_max * u_max * a * std::sin(a*x[1]) * std::cos(a*x[1]);
    return result;
}


double TGVErrorSolver::getPressure(Job* job, KinematicVector const &x){
    return -density * u_max * u_max * 0.25 * (std::cos(2.0*a*x[0]) + std::cos(2.0*a*x[1]));
}

void TGVErrorSolver::setGridLengths(Job* job){
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

/*----------------------------------------------------------------------------*/
//
void TGVErrorSolver::step(Job* job){
    //get grid dimensions on first step
    if (job->t < job->dt){
        setGridLengths(job);
    }

    //create map
    createMappings(job);

    //assign point pressures
    assignPressure(job);

    //calculate acceleration
    calculateAcceleration(job);

    //calculate error and write to file
    writeErrorInfo(job);

    //assign grid velocity
    assignVelocity(job);

    //move points
    movePoints(job);

    //calculate material strain
    calculateStrainRate(job);

    //update density and integrators
    updateDensity(job);

    return;
}

/*----------------------------------------------------------------------------*/
//
std::string TGVErrorSolver::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

/*----------------------------------------------------------------------------*/
//
int TGVErrorSolver::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//
void TGVErrorSolver::createMappings(Job *job){
    for (int b=0;b<job->bodies.size();b++){
        //job->bodies[b]->generateMap(job, DefaultBody::CPDI_ON); //use_cpdi by default
        job->bodies[b]->generateMap(job, 0);
    }
    return;
}

void TGVErrorSolver::assignPressure(Job* job){
    //temporary stress tensor
    MaterialTensor tmpT = MaterialTensor::Identity();
    double pressure = 0;

    for (int b=0; b<job->bodies.size(); b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        //loop over points and assign stresses
        for (int p=0; p<job->bodies[b]->points->x.size(); p++){
            //get pressure
            pressure = getPressure(job, job->bodies[b]->points->x[p]);

            //assign pressure
            tmpT = -pressure * MaterialTensor::Identity();
            job->bodies[b]->points->T[p] = tmpT;
        }
    }

    return;
}

void TGVErrorSolver::calculateAcceleration(Job* job){
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

        //map divergence of stress
        nodes->f.setZero();
        tmpMat = points->T;
        for (int i = 0; i < tmpMat.size(); i++) {
            tmpMat[i] *= points->v[i];
        }
        nvec = body->gradS.left_multiply_by_tensor(tmpMat); //f_i = v_p T_p * gradS_ip
        for (int i = 0; i < nvec.size(); i++) {
            nodes->f[i] -= KinematicVector(nvec[i], nodes->f.VECTOR_TYPE);
        }
    }
    return;
}

void TGVErrorSolver::assignVelocity(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //set velocity
        for (int i=0;i<job->bodies[b]->nodes->x_t.size();i++) {
            job->bodies[b]->nodes->x_t[i] = getVelocity(job, job->bodies[b]->nodes->x[i]);
        }

        //set displacement
        job->bodies[b]->nodes->u = job->dt * job->bodies[b]->nodes->x_t;
    }
    return;
}

void TGVErrorSolver::movePoints(Job* job){
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
        for (int p=0; p<points->x.size(); p++) {
            points->x_t[p] = getVelocity(job, points->x[p]);
        }

        //calculate momentum
        for (int i=0;i<points->mx_t.size();i++){
            points->mx_t(i) = points->m(i) * points->x_t(i);
        }
    }
    return;
}

void TGVErrorSolver::calculateStrainRate(Job* job){
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

void TGVErrorSolver::updateDensity(Job* job){
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


void TGVErrorSolver::writeErrorInfo(Job* job){
    //calculate and write error measures to output file for first body in simulation
    //BE CAREFUL!!!

    //initialize arrays
    Eigen::VectorXd V_i = Eigen::VectorXd(job->grid->node_count);
    Eigen::VectorXd v_i = Eigen::VectorXd(job->grid->node_count);
    Eigen::VectorXd e = Eigen::VectorXd(job->grid->node_count);
    Eigen::VectorXd H = Eigen::VectorXd(job->grid->node_count);
    Eigen::VectorXd G_rel = Eigen::VectorXd(job->grid->node_count);
    KinematicVectorArray g = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);

    //initialize figures of merit
    double H_norm = 0;
    double e_norm = 0;
    double a_L2 = 0;
    double g_max = 0;
    double G_norm = 0;

    //get exact node volumes form grid
    for (int i=0; i<job->grid->node_count;i++){
        V_i(i) = job->grid->nodeVolume(job,i);
    }

    //get integrated node volume from points
    v_i = job->bodies[0]->S * job->bodies[0]->points->v;

    //get integrated gradient from points
    g = job->bodies[0]->gradS * job->bodies[0]->points->v;

    //get node-wise integral of (x-x_i).grad(N_i(x))
    G_rel.setZero();
    int p = 0;
    int n = 0;
    KinematicVector dx = KinematicVector(job->JOB_TYPE);
    KinematicVector gradN = KinematicVector(job->JOB_TYPE);
    for (int j = 0; j<job->bodies[0]->gradS.size(); j++){
        n = job->bodies[0]->gradS.i_vec[j];
        p = job->bodies[0]->gradS.j_vec[j];
        gradN = job->bodies[0]->gradS.buffer[j];
        dx = (job->bodies[0]->points->x[p] - job->bodies[0]->nodes->x[n]);

        //check that dx is valid (for periodic domains)
        for (int pos=0; pos<job->grid->GRID_DIM; pos++){
            if (dx[pos] > Lx[pos]/2.0){
                dx[pos] -= Lx[pos];
            } else if (dx[pos] < -Lx[pos]/2.0){
                dx[pos] += Lx[pos];
            }
        }

        G_rel(n) += job->bodies[0]->points->v(p)*dx.dot(gradN);
    }
    for (int i=0; i<job->grid->node_count; i++){
        G_rel(i) += job->grid->GRID_DIM*V_i(i);
        G_rel(i) /= job->grid->GRID_DIM*V_i(i);
        G_rel(i) = std::min(G_rel(i), 0.0); //G_rel >=0 for all i
    }

    //calculate arrays
    double tmpNum;
    for (int i = 0; i < V_i.rows(); i++) {
        tmpNum = (v_i(i) - V_i(i));
        H(i) = std::max(0.0, tmpNum);
        e(i) = H(i) / V_i(i);
    }

    //calculate error measures
    KinematicVector tmpAcc = KinematicVector(job->JOB_TYPE);
    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    for (int i=0; i<V_i.rows(); i++){
        //||H||_2^2
        H_norm += H(i)*H(i);

        //||e||_\infty
        if (e(i) > e_norm){
            e_norm = e(i);
        }

        //||a^* - a^Q||_L2
        tmpAcc = getAcceleration(job, job->bodies[0]->nodes->x[i]);
        if (job->bodies[0]->nodes->m(i) > 0) {
            tmpVec = (tmpAcc - job->bodies[0]->nodes->f[i] / job->bodies[0]->nodes->m(i));
            a_L2 += tmpVec.dot(tmpVec) * V_i(i);
        }

        //||g||_\infty
        if (g[i].norm() > g_max){
            g_max = g[i].norm();
        }

        //grid dimension is 2
        if (G_rel(i) < G_norm){
            G_norm = G_rel(i);
        }
    }
    //sqrt of ||a_err||_L2^2
    a_L2 = std::sqrt(a_L2);

    //open and write to file
    std::ofstream file (output_filename,std::ios::app);
    if (file.is_open()){
        //success!
        //write to file
        file << job->t << ", ";
        file << H_norm << ", ";
        file << e_norm << ", ";
        file << a_L2 << ", ";
        file << g_max << ", ";
        file << -G_norm << "\n";

        file.close();
    } else {
        std::cerr << "ERROR! Cannot open " << output_filename << "!" << std::endl;;
    }

}

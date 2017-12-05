//
// Created by aaron on 11/21/17.
// explicit_usl_apic.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <driver.hpp>
#include <eigen3/Eigen/Dense>

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

int use_cpdi = 0;
int use_flip = 0;
double m_tol = 0;

class AffineParticle {
public:
    bool initialized;
    int id;
    Eigen::MatrixXd B;
    Eigen::MatrixXd D;
    Eigen::MatrixXd Dinv;
    Eigen::MatrixXd tmpNodeTensor; //will be used frequently, preallocating
    Eigen::MatrixXd tmpNodeVector;
    Eigen::MatrixXd tmpPointTensor;
    Eigen::MatrixXd tmpPointVector;
    AffineParticle();
    void affineparticleInit(Job *job, int idIN);
    void affineparticleCalculateB(Job* job, int use_flip_IN);
    void affineparticleCalculateD(Job* job);
};

std::vector<AffineParticle> affineparticles(0);

extern "C" void solverInit(Job* job); //initialize solver
extern "C" void solverStep(Job* job); //one forward mpm step

extern "C" std::string solverSaveState(Job* job, Serializer* serializer, std::string filepath); //save solver state to returned filename in serializer folder
extern "C" int solverLoadState(Job* job, Serializer* serializer, std::string fullpath); //load state from given full path

/*----------------------------------------------------------------------------*/

AffineParticle::AffineParticle():
        //set matrices to zeros
        B(0,0),
        D(0,0),
        Dinv(0,0),
        tmpNodeTensor(0,0),
        tmpNodeVector(0,0),
        tmpPointTensor(0,0),
        tmpPointVector(0,0)
{
    //body is not initialized
    initialized = 0;
    id = -1;
}

void AffineParticle::affineparticleCalculateB(Job *job, int use_flip_IN) {
    Eigen::MatrixXd tmpMat = job->jobTensor<double>();
    Eigen::VectorXd tmpVec(tmpMat.size());

    std::vector<int> n_id(0);
    std::vector<double> n_val(0);

    for (size_t p=0; p<B.rows(); p++){
        n_id.clear();
        n_val.clear();
        tmpMat.setZero();
        job->grid.gridEvaluateShapeFnValue(job, job->bodies[id].points.x.row(p), n_id, n_val);
        for (size_t i = 0; i<n_id.size(); i++) {
            for (size_t j = 0; j < tmpMat.rows(); j++) {
                for (size_t k = 0; k < tmpMat.cols(); k++) {
                    if (use_flip_IN == 1) {
                        tmpMat(j, k) +=
                                n_val[i] * (job->bodies[id].nodes.diff_x_t(n_id[i], j)) *
                                (job->bodies[id].nodes.x(n_id[i], k) - job->bodies[id].points.x(p, k));
                    } else {
                        tmpMat(j, k) +=
                                n_val[i] * (job->bodies[id].nodes.x_t(n_id[i], j)) *
                                (job->bodies[id].nodes.x(n_id[i], k) - job->bodies[id].points.x(p, k));
                    }
                }
            }
        }
        for (size_t pos=0;pos<tmpVec.size();pos++){
            tmpVec(pos) = tmpMat(pos);
        }
        if (use_flip_IN==1){
            B.row(p) += tmpVec.transpose();
        } else {
            B.row(p) = tmpVec.transpose();
        }
    }

    /*
    for (size_t i=0; i < tmpNodeTensor.rows(); i++){
        //self code tensor product
        for (size_t j=0; j<tmpMat.rows(); j++){
            for (size_t k=0; k<tmpMat.cols(); k++){
                if (use_flip_IN == 0) {
                    //set B from velocity directly
                    tmpMat(j, k) = job->bodies[id].nodes.x_t(i,j) * (job->bodies[id].nodes.x(i,k));
                } else if (use_flip_IN == 1){
                    //update B from difference in velocity
                    tmpMat(j,k) = job->bodies[id].nodes.diff_x_t(i,j) * (job->bodies[id].nodes.x(i,k));
                }
            }
        }
        for (size_t pos=0;pos<tmpVec.size();pos++){
            tmpVec(pos) = tmpMat(pos);
        }
        tmpNodeTensor.row(i) = tmpVec.transpose();
    }

    if (use_flip_IN == 0) {
        job->bodies[id].bodyCalcPointValues(job, B, tmpNodeTensor, Body::SET);
        job->bodies[id].bodyCalcPointValues(job, tmpPointVector, job->bodies[id].nodes.x_t, Body::SET);
    } else if (use_flip_IN == 1) {
        job->bodies[id].bodyCalcPointValues(job, B, tmpNodeTensor, Body::ADD);
        job->bodies[id].bodyCalcPointValues(job, tmpPointVector, job->bodies[id].nodes.diff_x_t, Body::SET);
    }

    for (size_t p=0; p<B.rows(); p++){
        //self code tensor product
        for (size_t j=0; j<tmpMat.rows(); j++){
            for (size_t k=0; k<tmpMat.cols(); k++){
                //set B from velocity directly
                tmpMat(j, k) = tmpPointVector(p,j) * job->bodies[id].points.x(p,k);
            }
        }
        for (size_t pos=0;pos<tmpVec.size();pos++){
            tmpVec(pos) = tmpMat(pos);
        }
        B.row(p) = B.row(p) - tmpVec.transpose();
    }
     */

    return;
}

void AffineParticle::affineparticleCalculateD(Job *job) {
    Eigen::MatrixXd tmpMat = job->jobTensor<double>();
    Eigen::MatrixXd tmpInv = job->jobTensor<double>();
    Eigen::VectorXd tmpVec(tmpMat.size());

    std::vector<int> n_id(0);
    std::vector<double> n_val(0);

    for (size_t p=0; p<D.rows(); p++){
        n_id.clear();
        n_val.clear();
        tmpMat.setZero();
        job->grid.gridEvaluateShapeFnValue(job, job->bodies[id].points.x.row(p), n_id, n_val);
        for (size_t i = 0; i<n_id.size(); i++) {
            for (size_t j = 0; j < tmpMat.rows(); j++) {
                for (size_t k = 0; k < tmpMat.cols(); k++) {
                    tmpMat(j, k) += n_val[i]* (job->bodies[id].nodes.x(n_id[i],j) - job->bodies[id].points.x(p, j)) * (job->bodies[id].nodes.x(n_id[i],k) - job->bodies[id].points.x(p, k));
                }
            }
        }
        for (size_t pos=0;pos<tmpVec.size();pos++){
            tmpVec(pos) = tmpMat(pos);
        }
        D.row(p) = tmpVec.transpose();
        //std::cout << "a: " << D.row(p) << std::endl << std::endl;

        //invert D
        tmpInv = tmpMat.colPivHouseholderQr().solve(job->jobTensor<double>(Job::IDENTITY));
        for (size_t pos=0;pos<tmpVec.size();pos++){
            tmpVec(pos) = tmpInv(pos);
        }
        Dinv.row(p) = tmpVec.transpose();
    }

    /*
    for (size_t i=0; i < tmpNodeTensor.rows(); i++){
        //self code tensor product
        for (size_t j=0; j<tmpMat.rows(); j++){
            for (size_t k=0; k<tmpMat.cols(); k++){
                tmpMat(j,k) = (job->bodies[id].nodes.x(i,j)) * (job->bodies[id].nodes.x(i,k));
            }
        }
        for (size_t pos=0;pos<tmpVec.size();pos++){
            tmpVec(pos) = tmpMat(pos);
        }
        tmpNodeTensor.row(i) = tmpVec.transpose();
    }


    job->bodies[id].bodyCalcPointValues(job, D, tmpNodeTensor, Body::SET);

    for (size_t p=0; p<D.rows(); p++){
        //self code tensor product
        for (size_t j=0; j<tmpMat.rows(); j++){
            for (size_t k=0; k<tmpMat.cols(); k++){
                //set B from velocity directly
                tmpMat(j, k) = job->bodies[id].points.x(p,j) * job->bodies[id].points.x(p,k);
            }
        }

        for (size_t pos=0;pos<tmpVec.size();pos++){
            tmpVec(pos) = tmpMat(pos);
        }
        std::cout << "b: " << D.row(p) << std::endl;
        D.row(p) = D.row(p) - tmpVec.transpose();
        std::cout << "a: " << D.row(p) << std::endl << std::endl;

        //invert D
        tmpVec = D.row(p).transpose();
        tmpMat = job->jobTensor<double>(tmpVec.data());
        tmpInv = tmpMat.colPivHouseholderQr().solve(job->jobTensor<double>(Job::IDENTITY));
        for (size_t pos=0;pos<tmpVec.size();pos++){
            tmpVec(pos) = tmpInv(pos);
        }
        Dinv.row(p) = tmpVec.transpose();
    }
     */

    return;
}

void AffineParticle::affineparticleInit(Job *job, int idIN){
    if (initialized==0) {
        //states
        initialized = 1;
        id = idIN;
        //history dependent
        B = job->jobTensorArray<double>(job->bodies[id].points.x.rows());
        B.setZero();
        //time dependent
        D = job->jobTensorArray<double>(job->bodies[id].points.x.rows());
        D.setZero();
        Dinv = job->jobTensorArray<double>(job->bodies[id].points.x.rows());
        Dinv.setZero();
        //node vector
        tmpNodeTensor = job->jobTensorArray<double>(job->bodies[idIN].nodes.x.rows());
        tmpNodeVector = job->jobVectorArray<double>(job->bodies[idIN].nodes.x.rows());
        tmpPointTensor = job->jobTensorArray<double>(job->bodies[id].points.x.rows());
        tmpPointVector = job->jobVectorArray<double>(job->bodies[id].points.x.rows());

        //set initial B, D
        if (use_cpdi == 0) {
            job->bodies[id].bodyGenerateMap(job,Body::CPDI_OFF);
        } else {
            job->bodies[id].bodyGenerateMap(job,Body::CPDI_ON);
        }
        //map mass
        job->bodies[id].bodyCalcNodalValues(job,job->bodies[id].nodes.m,job->bodies[id].points.m,Body::SET);
        //map momentum
        for (size_t i=0;i<job->bodies[id].points.mx_t.cols();i++){
            job->bodies[id].points.mx_t.col(i) = job->bodies[id].points.m.array() * job->bodies[id].points.x_t.col(i).array();
        }
        job->bodies[id].bodyCalcNodalValues(job,job->bodies[id].nodes.mx_t,job->bodies[id].points.mx_t,Body::SET);
        //calculate velocity
        for (size_t i=0;i<job->bodies[id].nodes.x_t.rows();i++){
            if (job->bodies[id].nodes.m(i) > 0){
                job->bodies[id].nodes.x_t.row(i) = job->bodies[id].nodes.mx_t.row(i) / job->bodies[id].nodes.m(i);
            } else {
                job->bodies[id].nodes.x_t.row(i).setZero();
            }
        }
        //simple filtered version
        affineparticleCalculateB(job,0);
        affineparticleCalculateD(job);
    }
    return;
}

/*----------------------------------------------------------------------------*/

void solverInit(Job* job){
    if (job->solver.int_props.size() < 2){
        std::cout << job->solver.int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties (use_cpdi, use_flip) defined.\n",
                __FILE__, __func__);
        exit(0);
    } else if (job->solver.int_props.size() >= 2){
        //set body ids by name
        use_cpdi = job->solver.int_props[0];
        use_flip = job->solver.int_props[1];
        for (size_t b=0;b<job->bodies.size();b++){
            affineparticles.push_back(AffineParticle());
        }
        if (job->solver.fp64_props.size() >= 1){
            m_tol = job->solver.fp64_props[0];
            std::cout << "DETECTED MASS TOLERANCE INPUT TO SOLVER: m_tol = " << m_tol << "!" << std::endl;
        }
    }
    std::cout << "Solver properties: (use_cpdi = " << use_cpdi << ", use_flip = " << use_flip << ")." << std::endl;
    std::cout << "Solver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void checkAffinePoints(Job* job){
    for (size_t b=0; b<job->bodies.size(); b++){
        if (affineparticles[b].initialized == 0) {
            affineparticles[b].affineparticleInit(job,b);
        }
    }
}

void createMappings(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (use_cpdi == 1) {
            job->bodies[b].bodyGenerateMap(job, Body::CPDI_ON);
        } else {
            job->bodies[b].bodyGenerateMap(job, Body::CPDI_OFF);
        }
    }
    return;
}

void mapPointsToNodes(Job* job){
    Body *body;
    Points *points;
    Nodes *nodes;
    Eigen::MatrixXd pvec;
    Eigen::MatrixXd nvec;
    Eigen::MatrixXd tmpB = job->jobTensor<double>();
    Eigen::MatrixXd tmpDinv = job->jobTensor<double>();
    Eigen::MatrixXd tmpMat = job->jobTensor<double>();
    Eigen::VectorXd tmpVec;
    Eigen::VectorXd tmpVec2 = job->jobVector<double>();

    std::vector<int> n_id(0);
    std::vector<double> n_val(0);

    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = &(job->bodies[b]);
        points = &(job->bodies[b].points);
        nodes = &(job->bodies[b].nodes);

        //map mass
        body->bodyCalcNodalValues(job,nodes->m,points->m,Body::SET);

        //setup B*Dinv*(x_i - x_p)
        /*
        affineparticles[b].affineparticleCalculateD(job);
        for (size_t p=0;p<points->x.rows();p++){
            tmpVec = affineparticles[b].B.row(p).transpose();
            tmpB = job->jobTensor<double>(tmpVec.data());

            tmpVec = affineparticles[b].Dinv.row(p).transpose();
            tmpDinv = job->jobTensor<double>(tmpVec.data());

            std::cout << "B" << std::endl;
            std::cout << tmpB << std::endl << std::endl;
            std::cout << "Dinv" << std::endl;
            std::cout << tmpDinv << std::endl << std::endl;

            tmpMat = points->m[p]*tmpB*tmpDinv;
            for (size_t pos=0;pos<tmpVec.size();pos++){
                tmpVec(pos) = tmpMat(pos);
            }
            affineparticles[b].tmpPointTensor.row(p) = tmpVec.transpose();

            tmpVec2 = -points->m[p]*tmpB*tmpDinv*points->x.row(p).transpose();
            affineparticles[b].tmpPointVector.row(p) = tmpVec2.transpose();
        }

        //map momentum
        for (size_t i=0;i<points->mx_t.cols();i++){
            points->mx_t.col(i) = points->m.array() * points->x_t.col(i).array();
        }
        body->bodyCalcNodalValues(job,nodes->mx_t,points->mx_t,Body::SET);
        body->bodyCalcNodalValues(job,nodes->mx_t,affineparticles[b].tmpPointVector,Body::ADD);
        body->bodyCalcNodalValues(job,affineparticles[b].tmpNodeTensor,affineparticles[b].tmpPointTensor,Body::SET);
        for (size_t i=0;i<nodes->x.rows();i++){
            tmpVec = affineparticles[b].tmpNodeTensor.row(i).transpose();
            tmpMat = job->jobTensor<double>(tmpVec.data());
            tmpVec2 = tmpMat*nodes->x.row(i).transpose();
            nodes->mx_t.row(i) = nodes->mx_t.row(i) + tmpVec2.transpose();
        }
        */

        //map momentum
        for (size_t i=0;i<points->mx_t.cols();i++){
            points->mx_t.col(i) = points->m.array() * points->x_t.col(i).array();
        }
        body->bodyCalcNodalValues(job,nodes->mx_t,points->mx_t,Body::SET);

        //map affine velocity field
        /*
        affineparticles[b].affineparticleCalculateD(job);
        for (size_t p=0;p<points->x.rows();p++){
            n_id.clear();
            n_val.clear();
            tmpMat.setZero();
            job->grid.gridEvaluateShapeFnValue(job, points->x.row(p), n_id, n_val);

            tmpVec = affineparticles[b].B.row(p).transpose();
            tmpB = job->jobTensor<double>(tmpVec.data());

            tmpVec = affineparticles[b].Dinv.row(p).transpose();
            tmpDinv = job->jobTensor<double>(tmpVec.data());

            //std::cout << "B" << std::endl;
            //std::cout << tmpB << std::endl << std::endl;
            //std::cout << "Dinv" << std::endl;
            //std::cout << tmpDinv << std::endl << std::endl;

            for (size_t i=0; i<n_id.size(); i++) {
                nodes->mx_t.row(n_id[i]) += n_val[i] * points->m[p] * (tmpB * tmpDinv * (nodes->x.row(n_id[i]) - points->x.row(p)).transpose()).transpose();
            }
        }
        */
        for (size_t p=0;p<points->x.rows();p++){
            n_id.clear();
            n_val.clear();
            tmpMat.setZero();
            job->grid.gridEvaluateShapeFnValue(job, points->x.row(p), n_id, n_val);

            tmpVec = job->bodies[b].points.L.row(p).transpose();
            tmpMat = job->jobTensor<double>(tmpVec.data());

            for (size_t i=0; i<n_id.size(); i++) {
                if (nodes->m(n_id[i]) > m_tol) {
                    nodes->mx_t.row(n_id[i]) += n_val[i] * points->m[p] * (tmpMat * (nodes->x.row(n_id[i]) -
                                                                                     points->x.row(
                                                                                             p)).transpose()).transpose();
                }
            }
        }

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
        job->contacts[c].contactApplyRules(job,Contact::EXPLICIT);
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
        job->bodies[b].nodes.u = job->dt * job->bodies[b].nodes.x_t;

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

        //affineparticles[b].affineparticleCalculateB(job,0); //use_flip);

        //map nodal displacement to point positions
        body->bodyCalcPointValues(job,points->x,nodes->u,Body::ADD);
        body->bodyCalcPointValues(job,points->u,nodes->u,Body::ADD);

        if (use_flip == 1) {
            //map nodal velocity diff to points
            body->bodyCalcPointValues(job, points->x_t, nodes->diff_x_t, Body::ADD);
        } else {
            body->bodyCalcPointValues(job, points->x_t, nodes->x_t, Body::SET);
        }

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

void solverStep(Job* job){
    //ensure affine particle fields are intialized
    //checkAffinePoints(job);

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
    //nothing to save yet
    //need to rewrite this if it works
    std::cout << "WARNING! explicit_usl_apic.so IS NOT SETUP TO LOAD OR SAVE STATE! UNDEFINED BEHAVIOR EXPECTED!" << std::endl;
    return "";
}
int solverLoadState(Job* job, Serializer* serializer, std::string fullpath){
    //just initialize it from scratch...
    std::cout << "WARNING! explicit_usl_apic.so IS NOT SETUP TO LOAD OR SAVE STATE! UNDEFINED BEHAVIOR EXPECTED!" << std::endl;
    if (job->solver.int_props.size() < 2){
        std::cout << job->solver.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties (use_cpdi, use_flip) defined.\n",
                __FILE__, __func__);
        exit(0);
    } else if (job->solver.int_props.size() >= 2){
        //set body ids by name
        use_cpdi = job->solver.int_props[0];
        use_flip = job->solver.int_props[1];
        for (size_t b=0;b<job->bodies.size();b++){
            affineparticles.push_back(AffineParticle());
        }
    }
    std::cout << "Solver properties: [" << use_cpdi << ", " << use_flip << "]." << std::endl;
    std::cout << "Solver Loaded." << std::endl;
    return 1;
}

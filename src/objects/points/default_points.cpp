//
// Created by aaron on 5/14/18.
// default_points.cpp
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>
#include <sstream>
#include <Eigen/Core>
#include <ctime>

#include "mpm_objects.hpp"
#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "job.hpp"

#include "points.hpp"
#include "objects/bodies/bodies.hpp"

/*----------------------------------------------------------------------------*/
//initialize point state (assumes that readFromFile has been called)
//no safety check on this, so be careful please
void DefaultPoints::init(Job* job, Body* body){
    //v0 initialization
    v0 = v;

    //extent initialization (for grid integration)
    if(job->grid->GRID_DIM == 1){
        for (int i=0;i<v.rows();i++){
            extent[i] = 0.5 * v[i];
        }
    } else if (job->grid->GRID_DIM == 2){
        for (int i=0;i<v.rows();i++){
            extent[i] = 0.5 * std::sqrt(v[i]);
        }
    } else if (job->grid->GRID_DIM == 3){
        for (int i = 0; i < v.rows(); i++) {
            extent[i] = 0.5 * std::cbrt(v[i]);
        }
    }

    //A matrix initialization
    int cp = 1;
    for (int i=0;i<job->grid->GRID_DIM;i++){
        cp *= 2; //square or cube
    }
    //initialize A matrix
    //for mapping position in cube to id
    //0 -> -1,-1,-1
    //1 -> +1,-1,-1
    //...
    //8 -> +1,+1,+1
    Eigen::VectorXi onoff = -1*Eigen::VectorXi::Ones(job->grid->GRID_DIM);
    A = Eigen::MatrixXi(cp, job->grid->GRID_DIM);
    for (int c=0; c<cp;c++){
        for (int i=0;i<onoff.rows();i++){
            A(c,i) = onoff(i);
        }
        for (int i=0;i<onoff.rows();i++) {
            if (onoff(i) == -1){
                onoff(i) = 1;
                break;
            } else {
                onoff(i) = -1;
            }
        }
    }

    std::cout << "Points Initialized: [" << file << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//read point data from file
void DefaultPoints::readFromFile(Job *job, Body *body, std::string fileIN) {
    //first line lists number of points in file
    //subsequent lines contain m, v, x, x_t, active
    file = fileIN;

    std::string line;
    std::ifstream fin(file);
    std::stringstream ss;
    std::vector<std::string> s_vec;

    if (fin.is_open()) {
        std::getline(fin, line);
        int len = std::stoi(line);

        //size KinematicVectors
        x = KinematicVectorArray(len, job->JOB_TYPE);
        u = KinematicVectorArray(len, job->JOB_TYPE);
        x_t = KinematicVectorArray(len, job->JOB_TYPE);
        mx_t = KinematicVectorArray(len, job->JOB_TYPE);
        b = KinematicVectorArray(len, job->JOB_TYPE);

        //size scalar vectors
        m.resize(len);
        v.resize(len);
        v0.resize(len);
        active.resize(len);
        extent.resize(len);

        //size tensor arrays
        T = MaterialTensorArray(len);
        L = KinematicTensorArray(len, job->JOB_TYPE);

        //zero out all entries to start
        x.setZero();
        u.setZero();
        x_t.setZero();
        m.setZero();
        v.setZero();
        v0.setZero();
        mx_t.setZero();
        b.setZero();
        T.setZero();
        L.setZero();
        active.setZero();
        extent.setZero();

        for (int i = 0; i < len; i++) {
            std::getline(fin, line);
            s_vec = Parser::splitString(line, ' ');
            if (s_vec.size() == (1 + 1 + job->DIM + job->DIM + 1)) {
                m[i] = std::stod(s_vec[0]);                 //first column is mass
                v[i] = std::stod(s_vec[1]);                 //second column is volume
                for (int pos = 0; pos < job->DIM; pos++) {
                    x(i, pos) = std::stod(s_vec[2 + pos]);    //following cols are position
                }
                for (int pos = 0; pos < job->DIM; pos++) {
                    x_t(i, pos) = std::stod(s_vec[2 + job->DIM + pos]);     //following cols are velocity
                }
                active[i] = std::stod(s_vec[2 + 2 * job->DIM]);

            } else if (s_vec.size() == (1 + 1 + job->grid->GRID_DIM + job->grid->GRID_DIM + 1)) {
                m[i] = std::stod(s_vec[0]);                 //first column is mass
                v[i] = std::stod(s_vec[1]);                 //second column is volume
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++) {
                    x(i, pos) = std::stod(s_vec[2 + pos]);    //following cols are position
                }
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++) {
                    x_t(i, pos) = std::stod(s_vec[2 + job->grid->GRID_DIM + pos]);     //following cols are velocity
                }
                active[i] = std::stod(s_vec[2 + 2 * job->grid->GRID_DIM]);

            } else {
                std::cerr << "ERROR: Unable to read line: " << file << std::endl;
                return;
            }
        }

    } else {
        std::cerr << "ERROR: Unable to open file: " << file << std::endl;
        return;
    }

    return;
}


void DefaultPoints::generateMap(Job* job, Body* body, int SPEC) {
    //the default body will defer to points to generate the mapping

    body->S.clear();
    body->gradS.clear();

    //calculate phi and grad phi
    std::vector<int> nvec(0);
    std::vector<double> valvec(0);
    KinematicVectorArray gradvec(0,job->JOB_TYPE);
    KinematicVector tmpGrad(job->JOB_TYPE);
    KinematicVector tmpVec(job->JOB_TYPE);

    int ith_cpdi;

    for (int i = 0; i < x.size(); i++) {
        KinematicVector::Map x_i = x[i];
        ith_cpdi = SPEC;
        //check whether point is active:
        if (active(i) == 0) {
            continue;
        } else if (!job->grid->inDomain(job, x_i)) {
            active(i) = 0;
            continue;
        } else if (ith_cpdi == DefaultBody::CPDI_ON) {
            //check that corners are in domain
            for (int c = 0; c < A.rows(); c++) {
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                    tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos);
                }
                for (int pos = job->grid->GRID_DIM; pos < tmpVec.DIM; pos++){
                    tmpVec(pos) = 0;
                }
                if (!job->grid->inDomain(job, tmpVec)) {
                    //corner out of domain
                    ith_cpdi = DefaultBody::CPDI_OFF;
                    break;
                }
            }
        }

        //calculate map for ith point
        if (ith_cpdi == DefaultBody::CPDI_OFF) {
            nvec.resize(0);
            valvec.resize(0);
            job->grid->evaluateBasisFnValue(job, x_i, nvec, valvec);
            for (int j = 0; j < nvec.size(); j++) {
                body->S.push_back(nvec[j], i, valvec[j]); //node, point, value
            }

            nvec.resize(0);
            gradvec.resize(0);
            job->grid->evaluateBasisFnGradient(job, x_i, nvec, gradvec);
            for (int j = 0; j < nvec.size(); j++) {
                body->gradS.push_back(nvec[j], i, gradvec[j]); //node, point, value
            }
        } else if (ith_cpdi == DefaultBody::CPDI_ON) {
            nvec.resize(0);
            valvec.resize(0);
            for (int c = 0; c < A.rows(); c++) {
                //spread influence
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                    tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos);
                }
                for (int pos = job->grid->GRID_DIM; pos < tmpVec.DIM; pos++){
                    tmpVec(pos) = 0;
                }
                job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
            }
            for (int j = 0; j < nvec.size(); j++) {
                body->S.push_back(nvec[j], i, valvec[j] / A.rows()); //node, point, value
            }

            //average gradients along sides of extent
            nvec.resize(0);
            valvec.resize(0);
            gradvec.resize(0);
            for (int c = 0; c < A.rows(); c++) {
                //find shape function value at offset point location
                //add node ids to nodevec
                //add values to valvec
                valvec.resize(0);

                for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                    tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos);
                }
                for (int pos = job->grid->GRID_DIM; pos < tmpVec.DIM; pos++){
                    tmpVec(pos) = 0;
                }

                job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);

                for (int v = 0; v < valvec.size(); v++) {
                    //gradient contribution from corner
                    //G(x) = (S(x+a) - S(x-a))/(2a)
                    for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                        tmpGrad[pos] = valvec[v] / (A.rows() * 0.5 * extent(i)) * A(c,pos);
                    }
                    for (int pos = job->grid->GRID_DIM; pos < tmpGrad.DIM; pos++){
                        tmpGrad[pos] = 0;
                    }
                    gradvec.push_back(tmpGrad);
                }
            }
            for (int j = 0; j < nvec.size(); j++) {
                body->gradS.push_back(nvec[j], i, gradvec[j]); //node, point, value
            }
        } else {
            std::cerr << "Unrecognized Argument for use_cpdi in DefaultBody::generateMap(): " << SPEC << std::endl;
        }
    }
    return;
}


/*----------------------------------------------------------------------------*/
void DefaultPoints::writeHeader(Job *job, Body *body, Serializer *serializer, std::ofstream &pfile, int SPEC) {
    //use default header
    serializer->writeDefaultPointHeader(job, body, pfile, SPEC);
    return;
}

//write relavent point data to file
void DefaultPoints::writeFrame(Job* job, Body* body, Serializer* serializer){
    //serializer will use x-position to create format for file
    //serializer.serializerWriteVectorArray(&x, "position")
    serializer->writeVectorArray(u,"displacement");
    serializer->writeVectorArray(x_t,"velocity");
    serializer->writeScalarArray(m,"mass");
    serializer->writeScalarArray(v,"volume");
    serializer->writeVectorArray(mx_t,"momentum");
    serializer->writeVectorArray(b,"body_force");
    serializer->writeTensorArray(T,"cauchy_stress");
    serializer->writeTensorArray(L,"velocity_gradient");
    //need to make double
    Eigen::VectorXd tmpVec = active.cast<double>();
    serializer->writeScalarArray(tmpVec,"active");
    serializer->writeScalarArray(extent,"extent");

    //pressure
    for(int i=0;i<T.size();i++){
        tmpVec(i) = -1.0/3.0 * T[i].trace();
    }
    serializer->writeScalarArray(tmpVec,"pressure");

    //density
    tmpVec = m.array() / v.array();
    serializer->writeScalarArray(tmpVec,"density");

    return;
}

void DefaultPoints::updateIntegrators(Job* job, Body* body){
    //do nothing
    return;
}


/*----------------------------------------------------------------------------*/
//not implemented yet
std::string DefaultPoints::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int DefaultPoints::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}

/*----------------------------------------------------------------------------*/
//
void DefaultPoints::generateLoads(Job* job, Body* body){
    //do nothing
    return;
}

void DefaultPoints::applyLoads(Job* job, Body* body){
    //do nothing
    return;
}
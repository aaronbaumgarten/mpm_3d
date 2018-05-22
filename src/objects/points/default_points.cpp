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

/*----------------------------------------------------------------------------*/
//initialize point state (assumes that readFromFile has been called)
//no safety check on this, so be careful please
void DefaultPoints::init(Job* job, Body* body){
    //v0 initialization
    v0 = v;

    //extent initialization
    if(job->DIM == 1){
        for (size_t i=0;i<v.rows();i++){
            extent[i] = 0.5 * v[i];
        }
    } else if (job->DIM == 2){
        for (size_t i=0;i<v.rows();i++){
            extent[i] = 0.5 * std::sqrt(v[i]);
        }
    } else if (job->DIM == 3){
        for (size_t i = 0; i < v.rows(); i++) {
            extent[i] = 0.5 * std::cbrt(v[i]);
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
        size_t len = std::stoi(line);

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

        for (size_t i = 0; i < len; i++) {
            std::getline(fin, line);
            s_vec = Parser::splitString(line, ' ');
            if (s_vec.size() < (1 + 1 + job->DIM + job->DIM + 1)){
                std::cerr << "ERROR: Unable to read line: " << file << std::endl;
                return;
            }

            m[i] = std::stod(s_vec[0]);                 //first column is mass
            v[i] = std::stod(s_vec[1]);                 //second column is volume
            for (int pos = 0; pos < job->DIM; pos++){
                x(i, pos) = std::stod(s_vec[2 + pos]);    //following cols are position
            }
            for (int pos = 0; pos < job->DIM; pos++){
                x_t(i, pos) = std::stod(s_vec[2 + job->DIM + pos]);     //following cols are velocity
            }
            active[i] = std::stod(s_vec[2 + 2*job->DIM]);
        }

    } else {
        std::cerr << "ERROR: Unable to open file: " << file << std::endl;
        return;
    }

    return;
}

/*----------------------------------------------------------------------------*/
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


/*----------------------------------------------------------------------------*/
//not implemented yet
std::string DefaultPoints::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int DefaultPoints::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}
//
// Created by aaron on 12/21/18.
// delta_points.cpp
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

void DeltaPoints::init(Job* job, Body* body){
    //v0 initialization
    v0 = v;

    //extent initialization (for grid integration)
    if(job->grid->GRID_DIM == 1){
        for (int i=0;i<v.rows();i++){
            extent[i] = 0.5 * v[i];
        }
    } else if (job->JOB_TYPE == job->JOB_AXISYM){
        for (int i=0;i<v.rows();i++){
            extent[i] = 0.5 * std::sqrt(v[i]/x(i,0));
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

    if (fp64_props.size() > 3){
        h = fp64_props[3];
        alpha = fp64_props[4];std::cout << "[" << h << ", " << alpha << "]." << std::endl;
    } else {
        h = 0;
        alpha = 0;
        std::cout << std::endl;
    }

    V_i.resize(job->grid->node_count,1);
    v_i.resize(job->grid->node_count,1);
    e.resize(job->grid->node_count,1);
    grad_e = KinematicVectorArray(body->points->x.size(),job->JOB_TYPE);

    for (int i=0; i<job->grid->node_count;i++){
        V_i(i) = job->grid->nodeVolume(job,i);
    }

    std::cout << "Points Initialized: [" << file << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//
void DeltaPoints::generateLoads(Job* job, Body* body){
    //do nothing
    return;
}

void DeltaPoints::applyLoads(Job* job, Body* body){
    //do nothing
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
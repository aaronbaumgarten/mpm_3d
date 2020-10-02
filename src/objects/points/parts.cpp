//
// Created by aaron on 11/30/18.
// parts.cpp
//

#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <Eigen/Core>

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

void Ball::init(Job *job){
    //initialize r, o
    if (fp64_props.size() < (job->grid->GRID_DIM+1)){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need %i values defined.\n",
                __FILE__, __func__, (job->grid->GRID_DIM+1));
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        o = KinematicVector(job->JOB_TYPE);
        for (int i=0; i<job->grid->GRID_DIM; i++){
            o[i] = fp64_props[i];
        }
        r = fp64_props[job->grid->GRID_DIM];

        std::cout << "Part properties (o = " << EIGEN_MAP_OF_KINEMATIC_VECTOR(o).transpose() << ", r = " <<  r  << ")." << std::endl;
    }
    return;
    return;
}

bool Ball::encompasses(KinematicVector &xIN) {
    if ((xIN - o).norm() <= r){
        return true;
    } else {
        return false;
    }
}

/*----------------------------------------------------------------------------*/

void Box::init(Job *job) {
    //initialize x_min, x_max
    if (fp64_props.size() < 2*job->grid->GRID_DIM){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need %i values defined.\n",
                __FILE__, __func__, job->grid->GRID_DIM);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        x_min = KinematicVector(job->JOB_TYPE);
        x_max = KinematicVector(job->JOB_TYPE);
        for (int i=0; i<job->grid->GRID_DIM; i++){
            x_min[i] = fp64_props[i];
            x_max[i] = fp64_props[i+job->grid->GRID_DIM];
        }

        std::cout << "Part properties (x_min = " << EIGEN_MAP_OF_KINEMATIC_VECTOR(x_min).transpose() << ", x_max = " <<  EIGEN_MAP_OF_KINEMATIC_VECTOR(x_max).transpose()  << ")." << std::endl;
    }
    return;
}

bool Box::encompasses(KinematicVector &xIN) {
    for (int i=0; i<xIN.DIM; i++){
        if (xIN(i) > x_max(i) || xIN(i) < x_min(i)){
            return false;
        }
    }
    return true;
}

/*----------------------------------------------------------------------------*/

void SineWave::init(Job* job){
    //initialize amplitude, wavelength, offset
    if (fp64_props.size() < 4){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need 4 values defined.\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store values
        amplitude = fp64_props[0];
        wavelength = fp64_props[1];
        x_0 = fp64_props[2];
        y_0 = fp64_props[3];

        std::cout << "Part properties (amplitude = " << amplitude
                  << ", wavelength = " << wavelength
                  << ", x_0 = " << x_0
                  << ", y_0 = " << y_0 << ")." << std::endl;
    }
    return;
}

bool SineWave::encompasses(KinematicVector &xIN) {
    assert(xIN.DIM >= 2 && "SineWave not implented in 1D.");
    //check whether y-dimension lies underneath function of x-dimension
    if ((xIN[1]-y_0) <= amplitude*std::sin(2.0*M_PI*(xIN[0] - x_0)/wavelength)){
        return true;
    } else {
        return false;
    }
}

/*----------------------------------------------------------------------------*/

void SandPile::init(Job* job){
    //initialize amplitude, wavelength, offset
    if (fp64_props.size() < job->grid->GRID_DIM + 1){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need %i values defined.\n",
                __FILE__, __func__, job->grid->GRID_DIM + 1);
        exit(0);
    } else {
        //store values
        peak = KinematicVector(job->JOB_TYPE);
        for (int i=0; i<job->grid->GRID_DIM; i++){
            peak[i] = fp64_props[i];
        }

        mu = fp64_props[job->grid->GRID_DIM];
        H_index = job->grid->GRID_DIM-1;

        std::cout << "Part properties (peak = ";
        for (int i=0; i<job->grid->GRID_DIM; i++){
            std::cout << peak[i] << ", ";
        }
        std::cout << "mu = " << mu << ")." << std::endl;
    }
    return;
}

bool SandPile::encompasses(KinematicVector &xIN) {
    //find relative position to peak
    KinematicVector xdif = xIN - peak;

    //find relative height
    double hdif = xdif[H_index];
    xdif[H_index] = 0;

    //find radial distance
    double rdif = xdif.norm();

    //check that point lies underneath mu slope
    return (hdif < -mu*rdif);
}
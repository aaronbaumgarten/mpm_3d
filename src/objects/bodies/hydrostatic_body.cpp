//
// Created by aaron on 6/18/18.
// hydrostatic_body.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "bodies.hpp"

/*----------------------------------------------------------------------------*/
//
void HydrostaticBody::init(Job* job){
    points->init(job,this);
    nodes->init(job,this);

    if (activeMaterial != 0) {
        material->init(job, this);
    }
    if (activeBoundary != 0) {
        boundary->init(job, this);
    }

    if (fp64_props.size() < 2){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined ({g, effective_density}).\n",
                __FILE__, __func__);
        exit(0);
    } else if (fp64_props.size() == 1 + job->grid->GRID_DIM) {
        //gravity vector given
        KinematicVector gvec = KinematicVector(job->JOB_TYPE);
        for (int ii = 0; ii<job->grid->GRID_DIM; ii++){
            gvec[ii] = fp64_props[ii];
        }
        g = gvec.norm();
        effective_density = fp64_props[job->grid->GRID_DIM];

        //set hydrostatic pressure
        double height = 0;
        double pressure;
        for (int i=0; i < points->x.size(); i++){
            if (-points->x(i).dot(gvec)/g > height || i == 0){
                height = -points->x(i).dot(gvec)/g;
            }
        }

        for (int i = 0; i < points->x.size(); i++) {
            pressure = (height*g + points->x(i).dot(gvec)) * effective_density;
            material->assignPressure(job, this, pressure, i, Material::UPDATE);
        }
    } else {
        g = fp64_props[0];
        effective_density = fp64_props[1];

        //set hydrostatic pressure
        double height = 0;
        double pressure;
        for (int i=0; i < points->x.size(); i++){
            if (points->x(i,job->grid->GRID_DIM - 1) > height){
                height = points->x(i,job->grid->GRID_DIM - 1);
            }
        }

        for (int i = 0; i < points->x.size(); i++) {
            pressure = (height - points->x(i, job->grid->GRID_DIM - 1)) * effective_density * g;
            material->assignPressure(job, this, pressure, i, Material::UPDATE);
        }
    }

    //assign vector type
    S = MPMScalarSparseMatrix(nodes->x.size(), points->x.size());
    gradS = KinematicVectorSparseMatrix(nodes->x.size(), points->x.size(), job->JOB_TYPE);


    printf("Body properties (g = %g, effective_density = %g)\n" ,g, effective_density);
    std::cout << "Body Initialized: [" << name << "]." << std::endl;
    return;
}
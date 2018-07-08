//
// Created by aaron on 7/5/18.
// launched_body.cpp
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
void LaunchedBody::init(Job* job){
    points->init(job,this);
    nodes->init(job,this);
    material->init(job,this);
    boundary->init(job,this);

    if (fp64_props.size() < 2){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined ({launch_speed, launch_velocity}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        launch_speed = fp64_props[0];
        launch_angle = fp64_props[1];

        //set hydrostatic pressure
        double height = 0;
        double pressure;
        for (int i=0; i < points->x.size(); i++){
            if (points->x(i,job->DIM - 1) > height){
                height = points->x(i,job->DIM - 1);
            }
        }

        double pi = std::acos(-1.0);
        //set all point velocities
        //velocity direction given by angle 'up' from x-axis
        for (int i = 0; i < points->x.size(); i++) {
            points->x_t(i,0) = launch_speed*std::cos(launch_angle*pi/180.0);
            points->x_t(i,job->DIM-1) = launch_speed*std::sin(launch_angle*pi/180.0);
            points->mx_t(i,0) = points->m(i) * points->x_t(i,0);
            points->mx_t(i,job->DIM-1) = points->m(i) * points->x_t(i,job->DIM-1);
        }
    }

    //assign vector type
    S = MPMScalarSparseMatrix(nodes->x.size(), points->x.size());
    gradS = KinematicVectorSparseMatrix(nodes->x.size(), points->x.size(), job->JOB_TYPE);


    printf("Body properties (launch_speed = %g, launch_angle = %g)\n" , launch_speed, launch_angle);
    std::cout << "Body Initialized: [" << name << "]." << std::endl;
    return;
}
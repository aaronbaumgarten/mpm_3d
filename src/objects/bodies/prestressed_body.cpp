//
// Created by aaron on 11/20/18.
// prestressed_body.cpp
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
void PrestressedBody::init(Job* job){
    points->init(job,this);
    nodes->init(job,this);

    if (activeMaterial != 0) {
        material->init(job, this);
    }
    if (activeBoundary != 0) {
        boundary->init(job, this);
    }

    if (fp64_props.size() < 6){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 6 properties defined ({s_xx, s_xy, s_xz, s_yy, s_yz, s_zz}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        T = MaterialTensor();
        T(0,0) = fp64_props[0];
        T(0,1) = fp64_props[1];
        T(0,2) = fp64_props[2];
        T(1,1) = fp64_props[3];
        T(1,2) = fp64_props[4];
        T(2,2) = fp64_props[5];

        T(1,0) = T(0,1);
        T(2,0) = T(0,2);
        T(2,1) = T(1,2);

        for (int i = 0; i < points->x.size(); i++) {
            material->assignStress(job, this, T, i, Material::UPDATE);
        }
    }

    //assign vector type
    S = MPMScalarSparseMatrix(nodes->x.size(), points->x.size());
    gradS = KinematicVectorSparseMatrix(nodes->x.size(), points->x.size(), job->JOB_TYPE);


    printf("Body properties (s_xx = %g, s_xy = %g, s_xz = %g, s_yy = %g, s_yz = %g, s_zz = %g)\n" , T(0,0), T(0,1), T(0,2), T(1,1), T(1,2), T(2,2));
    std::cout << "Body Initialized: [" << name << "]." << std::endl;
    return;
}
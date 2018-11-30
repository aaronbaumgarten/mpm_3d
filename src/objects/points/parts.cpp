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

bool Ball::encompasses(KinematicVector &xIN) {
    if ((xIN - o).norm() <= r){
        return true;
    } else {
        return false;
    }
}

/*----------------------------------------------------------------------------*/

bool Box::encompasses(KinematicVector &xIN) {
    for (int i=0; i<xIN.DIM; i++){
        if (xIN(i) > x_max(i) || xIN(i) < x_min(i)){
            return false;
        }
    }
    return true;
}
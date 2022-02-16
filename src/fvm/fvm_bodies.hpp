//
// Created by aaron on 12/23/19.
// fvm_bodies.hpp
//

#ifndef MPM_V3_FVM_BODIES_HPP
#define MPM_V3_FVM_BODIES_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <fstream>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "mpm_objects.hpp"
#include "fvm_objects.hpp"

/*
//finite volume body class
//represents a single continuum body filling fluid domain
class FiniteVolumeBody : public MPMObject{
public:
    //initialize from job and driver
    virtual void init(Job*, FiniteVolumeDriver*) = 0;

    //container for fluid fields
    KinematicTensorArray p_x;      //momentum gradient
    KinematicVectorArray p, rho_x; //momentum and density gradient
    Eigen::VectorXd rho, P, theta; //density, pressure, and temperature
    MaterialTensorArray tau;       //shear stress

    //container for solid phase fields in mixture
    KinematicVectorArray v_s;   //solid phase velocity field (defined on MPM grid)
    Eigen::VectorXd n;          //mixture porosity field (defined on MPM grid)
};
 */

class FVMDefaultBody : public FiniteVolumeBody{
public:
    FVMDefaultBody(){
        object_name = "FVMDefaultBody"; //set object name here
    }

    bool HYDROSTATIC_INITIALIZATION = false;

    //initialize from job and driver
    virtual void init(Job* job, FiniteVolumeDriver* driver);
};

class FVMRocketBody : public FVMDefaultBody{
public:
    FVMRocketBody(){
        object_name = "FVMRocketBody";
    }

    double h,r,rho_1,theta_1,v_1;

    virtual void init(Job* job, FiniteVolumeDriver* driver);
};

#endif //MPM_V3_FVM_BODIES_HPP

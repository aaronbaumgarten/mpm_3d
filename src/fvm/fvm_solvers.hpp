//
// Created by aaron on 12/23/19.
// fvm_solvers.hpp
//

#ifndef MPM_V3_FVM_SOLVERS_HPP
#define MPM_V3_FVM_SOLVERS_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
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
class FiniteVolumeSolver : public MPMObject{
public:
    //functions that must be implemented by every finite volume solver (not many)
    virtual void init(Job*, FiniteVolumeDriver*) = 0;                                        //initialize from Job
    virtual void step(Job*, FiniteVolumeDriver*) = 0;                                        //perform single mpm step
};
 */

class FVMDefaultSolver : public FiniteVolumeSolver{
public:
    FVMDefaultSolver(){
        object_name = "FVMDefaultSolver";
    }

    Eigen::VectorXd density_fluxes;
    KinematicVectorArray momentum_fluxes;

    virtual void init(Job* job, FiniteVolumeDriver* driver);                                        //initialize from Job
    virtual void step(Job* job, FiniteVolumeDriver* driver);                                        //perform single mpm step
};

#endif //MPM_V3_FVM_SOLVERS_HPP

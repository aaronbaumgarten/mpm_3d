//
// Created by aaron on 12/23/19.
// fvm_materials.hpp
//

#ifndef MPM_V3_FVM_MATERIALS_HPP
#define MPM_V3_FVM_MATERIALS_HPP

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
//finite volume material class
//represents a single continuum body filling fluid domain
class FiniteVolumeMaterial : public MPMObject{
public:
    //initialize from job and driver
    virtual void init(Job*, FiniteVolumeDriver*) = 0;

    //functions to calculate fluid fields
    virtual MaterialTensor getStress(Job*, FiniteVolumeDriver*, KinematicVector x, KinematicTensor L, double rho, double theta) = 0;
    virtual MaterialTensor getShearStress(Job*, FiniteVolumeDriver*, KinematicVector x, KinematicTensor L, double rho, double theta) = 0;
    virtual double getPressure(Job*, FiniteVolumeDriver*, KinematicVector x, double rho, double theta) = 0;
    virtual double getSpeedOfSound(Job*, FiniteVolumeDriver*, KinematicVector x, double rho, double theta) = 0;
    virtual void calculateElementPressures(Job*, FiniteVolumeDriver*) = 0;
    virtual void calculateElementShearStresses(Job*, FiniteVolumeDriver*) = 0;
    virtual void solveMaterialEquations(Job*, FiniteVolumeDriver*) = 0;
};
 */

class FVMBarotropicViscousFluid : public FiniteVolumeMaterial{
public:
    FVMBarotropicViscousFluid(){
        object_name = "FVMBarotropicViscousFluid";
    }

    //material properties
    double kappa, eta, rho_0;

    //initialize from job and driver
    virtual void init(Job* job, FiniteVolumeDriver* driver);

    //functions to calculate fluid fields
    virtual MaterialTensor getStress(Job* job, FiniteVolumeDriver* driver, KinematicVector x, KinematicTensor L, double rho, double theta);
    virtual MaterialTensor getShearStress(Job* job, FiniteVolumeDriver* driver, KinematicVector x, KinematicTensor L, double rho, double theta);
    virtual double getPressure(Job* job, FiniteVolumeDriver* driver, KinematicVector x, double rho, double theta);
    virtual double getSpeedOfSound(Job* job, FiniteVolumeDriver* driver, KinematicVector x, double rho, double theta);
    virtual void calculateElementPressures(Job* job, FiniteVolumeDriver* driver);
    virtual void calculateElementShearStresses(Job* job, FiniteVolumeDriver* driver);
    virtual void solveMaterialEquations(Job* job, FiniteVolumeDriver* driver);
};

#endif //MPM_V3_FVM_MATERIALS_HPP

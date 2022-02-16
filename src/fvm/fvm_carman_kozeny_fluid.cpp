//
// Created by aaron on 12/23/19.
// fvm_barotropic_viscous_fluid.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <job.hpp>
#include <objects/bodies/bodies.hpp>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "mpm_objects.hpp"
#include "fvm_objects.hpp"
#include "fvm_materials.hpp"



KinematicVector FVMCarmanKozenyFluid::getInterphaseDrag(Job* job, FiniteVolumeDriver* driver,
                                                             double rho,
                                                             const KinematicVector& v_f,
                                                             const KinematicVector& v_s,
                                                             double n){
    //permeability
    double C = 0;
    //Carman-Kozeny
    if (n < 1e-10 || ((1-n) < 1e-10) || eta < 1e-10){
        C = 0;
    } else {
        C = 180.0 * (1 - n) * (1 - n) * eta / (n * grain_diam * grain_diam);
    }

    /*
    if (C > 1e-10) {
        std::cout << "C: " << C;
    }

    //approximate limit for C in explicit simulations
    if (n >= 1e-10 && ((1-n) >= 1e-10)){
        //std::cout << C << ", " << C / (1 + job->dt*C*(1.0/rho + n/(rho * (1.0 - n)))) << std::endl;
        C /= (1 + job->dt*C*(1.0/rho + n/(rho * (1.0 - n))));
    }

    if (C > 1e-10) {
        std::cout << " ?= " << C << std::endl;
    }
     */

    return C * (v_s - v_f);
}

double FVMCarmanKozenyFluid::getInterphaseDragCoefficient(Job* job, FiniteVolumeDriver* driver,
                                                         double rho,
                                                         const KinematicVector& v_f,
                                                         const KinematicVector& v_s,
                                                         double n,
                                                         int SPEC){
    //permeability
    double C = 0;
    //Carman-Kozeny
    if (n < 1e-10 || ((1-n) < 1e-10) || eta < 1e-10){
        return 0;
    } else {
        C = 180.0 * (1 - n) * (1 - n) * eta / (n * grain_diam * grain_diam);
    }

    if (SPEC == REGULAR_DRAG) {
        return C;
    } else {
        double rho_s = solid_rho * (1.0 - n);
        return C*rho*rho_s / (rho*rho_s + job->dt * C * (rho + rho_s));
    }
}
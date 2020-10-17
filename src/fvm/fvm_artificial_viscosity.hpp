//
// Created by aaron on 10/14/20.
//

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

#ifndef MPM_V3_FVM_ARTIFICIAL_VISCOSITY_HPP
#define MPM_V3_FVM_ARTIFICIAL_VISCOSITY_HPP

namespace ArtificialViscosityCalculator {
    double getAVCoeff(Job* job,
                      FiniteVolumeDriver* driver,
                      int f);

    struct fluxVector {
        double rho;
        KinematicVector p;
        double rhoE;
    };

    fluxVector getArtificialViscosityFlux(Job* job,
                                      FiniteVolumeDriver* driver,
                                      int f);

    void writeFrame(Job* job, FiniteVolumeDriver* driver);
};

#endif //MPM_V3_FVM_ARTIFICIAL_VISCOSITY_HPP

//
// Created by aaron on 10/28/16.
// loading.hpp
//

#ifndef MPM_3D_LOADING_HPP
#define MPM_3D_LOADING_HPP

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "particle.hpp"
#include "process.hpp"

void initial_loads(job_t *job);
void time_varying_loads(job_t *job);

#endif //MPM_3D_LOADING_HPP

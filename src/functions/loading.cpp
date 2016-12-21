//
// Created by aaron on 10/28/16.
// loading.cpp
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "particle.hpp"
#include "process.hpp"
#include "loading.hpp"

#define G_MAG 9.81
#define RAMP_TIME 0.0

void initial_loads(job_t *job){
    for (size_t b=0;b<job->num_bodies;b++){
        job->bodies[b].particles.bx.setZero();
        job->bodies[b].particles.by.setZero();
        job->bodies[b].particles.bz.setZero();
    }
    return;
}
void time_varying_loads(job_t *job){
    if (job->use_3d==1) {
        double gravity;

        if ((RAMP_TIME + job->step_start_time) > job->t) {
            gravity = -G_MAG * ((job->t - job->step_start_time) / RAMP_TIME);
        } else {
            gravity = -G_MAG;
        }

        for (size_t b = 0; b < job->num_bodies; b++) {
            for (size_t i = 0; i < job->bodies[b].p; i++) {
                job->bodies[b].particles.bx[i] = 0;
                job->bodies[b].particles.by[i] = 0;
                job->bodies[b].particles.bz[i] = gravity;
            }
        }

        return;
    } else {
        return time_varying_loads2D(job);
    }
}
void time_varying_loads2D(job_t *job){
    double gravity;

    if ((RAMP_TIME+job->step_start_time) > job->t){
        gravity = -G_MAG * ((job->t - job->step_start_time) / RAMP_TIME);
    } else {
        gravity = -G_MAG;
    }

    for (size_t b=0;b<job->num_bodies;b++){
        for (size_t i=0; i<job->bodies[b].p; i++){
            job->bodies[b].particles.bx[i] = 0;
            job->bodies[b].particles.by[i] = gravity;
            job->bodies[b].particles.bz[i] = 0;
        }
    }

    return;
}
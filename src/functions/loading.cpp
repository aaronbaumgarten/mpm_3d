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

#define G_MAG 0.0//9.81
#define RAMP_TIME 0.0

void initial_loads(job_t *job){
    for (size_t b=0;b<job->num_bodies;b++){
        for (size_t i=0;i<job->bodies[b].p;i++){
            job->bodies[b].particles[i].bx[0] = 0;
            job->bodies[b].particles[i].by[0] = 0;
            job->bodies[b].particles[i].bz[0] = 0;
        }
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
                job->bodies[b].particles[i].bx[0] = 0;
                job->bodies[b].particles[i].by[0] = 0;
                job->bodies[b].particles[i].bz[0] = gravity;
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
            job->bodies[b].particles[i].bx[0] = 0;
            job->bodies[b].particles[i].by[0] = gravity;
            job->bodies[b].particles[i].bz[0] = 0;
        }
    }

    return;
}
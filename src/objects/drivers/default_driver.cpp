//
// Created by aaron on 5/10/18.
// default_driver.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <time.h>

#include "job.hpp"
#include "mpm_objects.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"

#include "drivers.hpp"

/*----------------------------------------------------------------------------*/

void DefaultDriver::init(Job* job){
    if (fp64_props.size() < 1) {
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 1 property defined (stop_time).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        stop_time = fp64_props[0];
        gravity = KinematicVector(job->JOB_TYPE);

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ")." << std::endl;
    }

    std::cout << "Driver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void DefaultDriver::run(Job* job) {
    int stepCount = 0;
    int frameCount = 0;

    struct timespec timeStart, timeFrame, timeFinish;
    clock_gettime(CLOCK_MONOTONIC, &timeStart);
    timeFrame = timeStart;
    //clock_t clockSim = clock();
    //clock_t clockFrame = clock();
    double tSim = 0;
    double tFrame = 0;

    //initialize gravity
    generateGravity(job);
    applyGravity(job);

    //run simulation until stop_time
    while (job->t <= stop_time){
        //run solver
        job->solver->step(job);
        //std::cout << "Step Completed [" << ++stepCount << "]." << std::flush;
        if (job->serializer->writeFrame(job) == 1) {
            //successful frame written
            //tFrame = (double)(clock() - clockFrame)/CLOCKS_PER_SEC;
            //tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
            clock_gettime(CLOCK_MONOTONIC,&timeFinish);
            tFrame = (timeFinish.tv_sec - timeFrame.tv_sec) + (timeFinish.tv_nsec - timeFrame.tv_nsec)/1000000000.0;
            tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
            timeFrame = timeFinish;
            printf("\33[2K");
            std::cout << "Frame Written [" << ++frameCount << "]. Time/Frame [" << tFrame << " s]. Elapsed Time [" << tSim << " s]." << std::flush;
        }
        std::cout << "\r";
        job->t += job->dt;
    }
    //tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
    clock_gettime(CLOCK_MONOTONIC,&timeFinish);
    tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
    std::cout << std::endl << std::endl << "Simulation Complete. Elapsed Time [" << tSim << "s]." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void DefaultDriver::generateGravity(Job* job) {
    double g = -9.81;
    gravity.setZero();
    if (job->JOB_TYPE == Job::JOB_1D){
        gravity(0) = g;
    } else if (job->JOB_TYPE == Job::JOB_2D){
        gravity(1) = g;
    } else if (job->JOB_TYPE == Job::JOB_3D){
        gravity(2) = g;
    } else if (job->JOB_TYPE == Job::JOB_AXISYM){
        gravity(1) = g;
    } else if (job->JOB_TYPE == Job::JOB_2D_OOP){
        gravity(1) = g;
    }
    return;
}

/*----------------------------------------------------------------------------*/

void DefaultDriver::applyGravity(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (int i=0;i<job->bodies[b]->points->b.size();i++){
            job->bodies[b]->points->b(i) = gravity;
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

std::string DefaultDriver::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "error";
}

/*----------------------------------------------------------------------------*/

int DefaultDriver::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

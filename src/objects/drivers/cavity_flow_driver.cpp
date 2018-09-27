//
// Created by aaron on 7/9/18.
// cavity_flow_driver.cpp
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

void CavityFlowDriver::init(Job* job){
    if (job->DIM != 2){
        std::cerr << "ERROR: CavityFlowDriver requires job->DIM = 2, got " << job->DIM << ". Exiting." << std::endl;
        exit(0);
    }

    if (fp64_props.size() < 2) {
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined (stop_time, v_set).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        stop_time = fp64_props[0];
        v_set = fp64_props[1];

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ", v_set = " << v_set << ")." << std::endl;
    }

    std::cout << "Driver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void CavityFlowDriver::run(Job* job) {
    int stepCount = 0;
    int frameCount = 0;

    struct timespec timeStart, timeFrame, timeFinish;
    clock_gettime(CLOCK_MONOTONIC, &timeStart);
    timeFrame = timeStart;
    //clock_t clockSim = clock();
    //clock_t clockFrame = clock();
    double tSim = 0;
    double tFrame = 0;

    //initialize flow parameters (all points in top element have set velocity)
    Ly = 0;
    hy = -1;
    for (int i=0; i<job->bodies[0]->nodes->x.size(); i++){
        if (job->bodies[0]->nodes->x(i,1) > Ly){
            if (hy == -1 || hy > (job->bodies[0]->nodes->x(i,1)-Ly)){
                hy = job->bodies[0]->nodes->x(i,1)-Ly;
            }
            Ly = job->bodies[0]->nodes->x(i,1);
        }
    }

    //run simulation until stop_time
    while (job->t <= stop_time){
        //set point velocities:
        for (int b=0; b<job->bodies.size(); b++){
            if (job->activeBodies[b] == 1){
                for (int i=0; i<job->bodies[b]->points->x.size(); i++){
                    if (job->bodies[b]->points->x(i,1) > (Ly - hy)){
                        job->bodies[b]->points->x_t(i,0) = v_set;
                        job->bodies[b]->points->mx_t(i,0) = v_set * job->bodies[b]->points->m(i);
                    }
                }
            }
        }

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

void CavityFlowDriver::generateGravity(Job* job) {
    return;
}

/*----------------------------------------------------------------------------*/

void CavityFlowDriver::applyGravity(Job* job){
    return;
}

/*----------------------------------------------------------------------------*/

std::string CavityFlowDriver::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "error";
}

/*----------------------------------------------------------------------------*/

int CavityFlowDriver::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

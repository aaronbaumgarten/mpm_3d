//
// Created by aaron on 10/26/18.
// ballistic_driver.cpp
//

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

void BallisticDriver::init(Job* job){
    if (fp64_props.size() < 1 + job->grid->GRID_DIM || (int_props.size() < 1 && str_props.size() < 1)) {
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 3 properties defined (stop_time, <velocity>) and (ballistic body name).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        stop_time = fp64_props[0];
        gravity = KinematicVector(job->JOB_TYPE);
        velocity = KinematicVector(job->JOB_TYPE);
        ballistic_id = -1;

        //set velocity
        for (int i=0; i<job->grid->GRID_DIM; i++){
            velocity[i] = fp64_props[i+1];
        }

        for (int i=job->grid->GRID_DIM; i<job->DIM; i++){
            velocity[i] = 0;
        }

        //assign ballistic body
        //set body ids by name
        if (str_props.size() == 1){
            for (int b = 0; b < job->bodies.size(); b++) {
                if (str_props[0].compare(job->bodies[b]->name) == 0){
                    ballistic_id = b;
                    break;
                }
            }
        }

        // or set body ids by int
        if (ballistic_id < 0){
            ballistic_id = int_props[0];
        }

        //set ballistic body properties
        job->activeBodies[ballistic_id] = 0;
        for (int i=0; i<job->bodies[ballistic_id]->points->x_t.size(); i++){
            job->bodies[ballistic_id]->points->x_t[i] = velocity;
        }

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ", velocity = " << EIGEN_MAP_OF_KINEMATIC_VECTOR(velocity) << ", body = " << job->bodies[ballistic_id]->name << ")." << std::endl;
    }

    std::cout << "Driver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void BallisticDriver::run(Job* job) {
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
        //check ballistic status
        if (job->activeBodies[ballistic_id] == 0 && job->t >= 0){
            job->activeBodies[ballistic_id] = 1;
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

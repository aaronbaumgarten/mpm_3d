//
// Created by aaron on 5/23/18.
// column_collapse_driver.cpp
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

void ColumnCollapseDriver::setPressure(Job *job) {
    //check that contact properties are set
    if (fp64_props.size() < 4 || (str_props.size() < 2 && int_props.size() < 2)){
        //need bodies
        std::cout << "fp64: " << fp64_props.size() << " < 4 || int: " << int_props.size() << " < 2 && str: " << str_props.size() << " < 2\n";
        std::cout << "WARNING: Not enough information provided to column_collapse_driver.so. No hydrostatic stress will be applied!" << std::endl;
        return;
    } else {
        std::vector<int> bodyIDs = {-1,-1};

        //set body ids by name
        if (str_props.size() == 2) {
            for (int i = 0; i < bodyIDs.size(); i++) {
                for (int b = 0; b < job->bodies.size(); b++) {
                    if (str_props[i].compare(job->bodies[b]->name) == 0) {
                        bodyIDs[i] = b;
                        break;
                    }
                }
            }
        }

        // or set body ids by int
        for (size_t i = 0; i < bodyIDs.size(); i++) {
            if (bodyIDs[i] < 0) {
                if (int_props.size() == 2) {
                    bodyIDs = int_props;
                } else {
                    std::cout << "fp64: " << fp64_props.size() << " < 4 || int: " << int_props.size() << " < 2 && str: " << str_props.size() << " < 2\n";
                    std::cout << "WARNING: Not enough information provided to column_collapse_driver.so. No hydrostatic stress will be applied!" << std::endl;
                    return;
                }
                break;
            }
        }

        int solid_body_id = bodyIDs[0];
        int fluid_body_id = bodyIDs[1];

        double grain_density = fp64_props[1];
        double packing_fraction = fp64_props[2];
        double liquid_density = fp64_props[3];

        double height, pressure;
        double g = 9.81;

        //set liquid phase pressure
        height = 0; //job->bodies[liquid_body_id]->points->x.col(job->DIM - 1).maxCoeff();
        for (int i=0; i < job->bodies[fluid_body_id]->points->x.size(); i++){
            if (job->bodies[fluid_body_id]->points->x(i,job->DIM - 1) > height){
                height = job->bodies[fluid_body_id]->points->x(i,job->DIM - 1);
            }
        }

        for (int i = 0; i < job->bodies[fluid_body_id]->points->x.size(); i++) {
            pressure = (height - job->bodies[fluid_body_id]->points->x(i, job->DIM - 1)) * liquid_density * g;
            job->bodies[fluid_body_id]->material->assignPressure(job, job->bodies[fluid_body_id].get(), pressure, i, Material::UPDATE);
        }

        //set solid phase pressure
        height = 0; //job->bodies[solid_body_id]->points->x.col(job->DIM - 1).maxCoeff();
        for (int i=0; i < job->bodies[solid_body_id]->points->x.size(); i++){
            if (job->bodies[solid_body_id]->points->x(i,job->DIM - 1) > height){
                height = job->bodies[solid_body_id]->points->x(i,job->DIM - 1);
            }
        }

        for (int i = 0; i < job->bodies[solid_body_id]->points->x.size(); i++) {
            pressure = (height - job->bodies[solid_body_id]->points->x(i, job->DIM - 1)) * packing_fraction * g * (grain_density - liquid_density);
            job->bodies[solid_body_id]->material->assignPressure(job, job->bodies[solid_body_id].get(), pressure, i, Material::UPDATE);
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void ColumnCollapseDriver::init(Job* job){
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

void ColumnCollapseDriver::run(Job* job) {
    size_t stepCount = 0;
    size_t frameCount = 0;

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

    //initialize pressure field
    setPressure(job);

    //run simulation until stop_time
    while (job->t < stop_time){
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

void ColumnCollapseDriver::generateGravity(Job* job) {
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
    }
    return;
}

/*----------------------------------------------------------------------------*/

void ColumnCollapseDriver::applyGravity(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (size_t i=0;i<job->bodies[b]->points->b.size();i++){
            job->bodies[b]->points->b(i) = gravity;
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

std::string ColumnCollapseDriver::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "error";
}

/*----------------------------------------------------------------------------*/

int ColumnCollapseDriver::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

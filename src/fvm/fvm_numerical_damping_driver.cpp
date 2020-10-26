//
// Created by aaron on 12/12/19.
// finite_volume_driver.cpp
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

#include "fvm_objects.hpp"
#include "fvm_drivers.hpp"
#include "fvm_artificial_viscosity.hpp"
#include "registry.hpp"

/*----------------------------------------------------------------------------*/

void FVMNumericalDampingDriver::init(Job* job){

    //call parent initializer
    FiniteVolumeDriver::init(job);

    int GRID_DIM = 0;
    //assign grid dimension from job type
    if (job->JOB_TYPE == job->JOB_1D){
        GRID_DIM = 1;
    } else if (job->JOB_TYPE == job->JOB_2D){
        GRID_DIM = 2;
    } else if (job->JOB_TYPE == job->JOB_3D){
        GRID_DIM = 3;
    } else if (job->JOB_TYPE == job->JOB_2D_OOP){
        GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
    } else if (job->JOB_TYPE == job->JOB_AXISYM){
        GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
    } else {
        std::cerr << "Job doesn't have defined type for input " << job->JOB_TYPE << "." << std::endl;
        exit(0);
    }

    //check size of properties passed to driver object
    if (fp64_props.size() <= 2+GRID_DIM) {
        //use default properties
        std::cerr << "Need 2 properties defined for FVMNumericalDampingDriver. Given " << fp64_props.size()-GRID_DIM-1 << std::endl;
        exit(0);
    } else {
        C0 = fp64_props[1+GRID_DIM];
        t0 = fp64_props[2+GRID_DIM];
    }

    //loop over str-props and assign relevant flags
    std::vector<std::string> options = {"AVOID_FLUID_DAMPING",
                                        "LOAD_SOLID_DIRECTLY"};
    for (int i=0; i<str_props.size(); i++){
        switch (Parser::findStringID(options, str_props[i])){
            case 0:
                //AVOID_FLUID_DAMPING
                AVOID_FLUID_DAMPING = true;
                std::cout << "Avoiding fluid damping." << std::endl;
                break;
            case 1:
                //LOAD_SOLID_DIRECTLY
                LOAD_SOLID_DIRECTLY = true;
                std::cout << "Loading solid directly. AXIAL TEST ONLY!!!" << std::endl;
                if (fp64_props.size() <= 2+3*GRID_DIM){
                    //
                    std::cerr << "Need additional properties defined for FVMNumericalDampingDriver. Given " << fp64_props.size()-GRID_DIM-2 << std::endl;
                    exit(0);
                } else {
                    xmax = KinematicVector(job->JOB_TYPE);
                    load = KinematicVector(job->JOB_TYPE);
                    for (int pos=0; pos<GRID_DIM; pos++){
                        xmax[pos] = fp64_props[3+GRID_DIM+pos];
                        load[pos] = fp64_props[3+2*GRID_DIM+pos];
                    }
                }
                break;
            default:
                //do nothing
                break;
        }
    }

    std::cout << "Using numerical damping. t0: " << t0 << "s, C0: " << C0 << " 1/s." << std::endl;
    if (LOAD_SOLID_DIRECTLY){
        std::cout << "Loading points w/ x < " << EIGEN_MAP_OF_KINEMATIC_VECTOR(xmax).transpose() << " with load: ";
        std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(load).transpose() << std::endl;
    }

    return;
}

/*----------------------------------------------------------------------------*/

void FVMNumericalDampingDriver::run(Job* job) {
    //initialize FVM objects
    serializer->init(job, this);
    fluid_grid->init(job, this);
    solver->init(job, this);
    fluid_material->init(job, this);
    fluid_body->init(job, this);

    //set counters to zero
    int stepCount = 0;
    int frameCount = 0;

    struct timespec timeStart, timeFrame, timeFinish, timeStep;
    clock_gettime(CLOCK_MONOTONIC, &timeStart);
    timeFrame = timeStart;
    timeStep = timeStart;
    //clock_t clockSim = clock();
    //clock_t clockFrame = clock();
    double tSim = 0;
    double tFrame = 0;
    double tStep = 0;

    //initialize gravity
    generateGravity(job);
    applyGravity(job);

    //damping scale
    double C = C0;
    double scale = 1.0;

    //run simulation until stop_time
    while (job->t <= stop_time){

        //run solver
        solver->step(job,this);

        //apply increment of damping
        C = C0 * std::exp(-job->t/t0);
        scale = (2.0 - job->dt*C)/(2.0 + job->dt*C);
        if (!AVOID_FLUID_DAMPING) {
            for (int e = 0; e < fluid_grid->element_count; e++) {
                fluid_body->p(e) *= scale;
            }
        }
        for (int b=0; b<job->bodies.size(); b++){
            for (int i=0; i<job->bodies[b]->points->x.size(); i++){
                job->bodies[b]->points->x_t[i] *= scale;
                job->bodies[b]->points->mx_t[i] *= scale;
            }
        }

        //clock_gettime(CLOCK_MONOTONIC, &timeFinish);
        //tStep = (timeFinish.tv_sec - timeStep.tv_sec) + (timeFinish.tv_nsec - timeStep.tv_nsec)/1000000000.0;
        //timeStep = timeFinish;

        //std::cout << "Step Completed [" << ++stepCount << "]. Step Time [" << tStep << " s]." << std::flush;
        if (job->serializer->writeFrame(job) == 1) {
            //call fvm serializer to write frame as well:
            serializer->writeFrame(job,this);

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


KinematicVector FVMNumericalDampingDriver::getSolidLoading(Job *job, const KinematicVector &x){
    bool ADD_LOAD = true;
    if (!LOAD_SOLID_DIRECTLY){
        //do nothing
    } else {
        for (int pos=0; pos<job->grid->GRID_DIM; pos++){
            if (x[pos] > xmax[pos]){
                ADD_LOAD = false;
            }
        }
        if (ADD_LOAD){
            return load;
        }
    }
    return KinematicVector(job->JOB_TYPE);
}
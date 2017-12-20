//
// Created by aaron on 12/18/17.
// column_collapse_driver.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <time.h>
#include "signal_resolution.hpp"

#include "job.hpp"

#include "serializer.hpp"
#include "driver.hpp"
#include "solver.hpp"

#include "body.hpp"
#include "points.hpp"
#include "nodes.hpp"

double stop_time;
Eigen::VectorXd gravity;

extern "C" void driverInit(Job* job); //initialize driver
extern "C" void driverRun(Job* job); //run simulation
extern "C" void driverGenerateGravity(Job* job); //apply gravity to simulation
extern "C" void driverApplyGravity(Job* job); //apply gravity to simulation

extern "C" std::string driverSaveState(Job* job, Serializer* serializer, std::string filepath); //save driver state to returned filename in serializer folder
extern "C" int driverLoadState(Job* job, Serializer* serializer, std::string fullpath); //load state from given full path

/*----------------------------------------------------------------------------*/

void setPressure(Job* job){

    //check that contact properties are set
    if (job->driver.fp64_props.size() < 4 || (job->driver.str_props.size() < 2 && job->driver.int_props.size() < 2)){
        //need bodies
        std::cout << "fp64: " << job->driver.fp64_props.size() << " < 4 || int: " << job->driver.int_props.size() << " < 2 && str: " << job->driver.str_props.size() << " < 2\n";
        std::cout << "WARNING: Not enough information provided to column_collapse_driver.so. No hydrostatic stress will be applied!" << std::endl;
        return;
    } else {
        std::vector<int> bodyIDs = {-1,-1};

        //set body ids by name
        if (job->driver.str_props.size() == 2) {
            for (size_t i = 0; i < bodyIDs.size(); i++) {
                for (size_t b = 0; b < job->bodies.size(); b++) {
                    if (job->driver.str_props[i].compare(job->bodies[b].name) == 0) {
                        bodyIDs[i] = b;
                        break;
                    }
                }
            }
        }

        // or set body ids by int
        for (size_t i = 0; i < bodyIDs.size(); i++) {
            if (bodyIDs[i] < 0) {
                if (job->driver.int_props.size() == 2) {
                    bodyIDs = job->driver.int_props;
                } else {
                    std::cout << "fp64: " << job->driver.fp64_props.size() << " < 4 || int: " << job->driver.int_props.size() << " < 2 && str: " << job->driver.str_props.size() << " < 2\n";
                    std::cout << "WARNING: Not enough information provided to column_collapse_driver.so. No hydrostatic stress will be applied!" << std::endl;
                    return;
                }
                break;
            }
        }

        int solid_body_id = bodyIDs[0];
        int liquid_body_id = bodyIDs[1];

        double grain_density = job->driver.fp64_props[1];
        double packing_fraction = job->driver.fp64_props[2];
        double liquid_density = job->driver.fp64_props[3];

        double height, pressure;
        double g = 9.81;

        //set liquid phase pressure
        height = job->bodies[liquid_body_id].points.x.col(job->DIM - 1).maxCoeff();
        for (size_t i = 0; i < job->bodies[liquid_body_id].points.x.rows(); i++) {
            pressure = (height - job->bodies[liquid_body_id].points.x(i, job->DIM - 1)) * liquid_density * g;
            job->bodies[liquid_body_id].material.materialAssignPressure(job, &job->bodies[liquid_body_id], pressure, i, Material::UPDATE);
        }

        //set solid phase pressure
        height = job->bodies[solid_body_id].points.x.col(job->DIM - 1).maxCoeff();
        for (size_t i = 0; i < job->bodies[solid_body_id].points.x.rows(); i++) {
            pressure = (height - job->bodies[solid_body_id].points.x(i, job->DIM - 1)) * packing_fraction * g * (grain_density - liquid_density);
            job->bodies[solid_body_id].material.materialAssignPressure(job, &job->bodies[solid_body_id], pressure, i, Material::UPDATE);
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void driverInit(Job* job){
    if (job->driver.fp64_props.size() < 1) {
        std::cout << job->driver.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 1 property defined (stop_time).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        stop_time = job->driver.fp64_props[0];
        gravity = job->jobVector<double>(Job::ZERO);

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ")." << std::endl;
    }

    std::cout << "Driver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void driverRun(Job* job) {
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
    driverGenerateGravity(job);
    driverApplyGravity(job);

    //initialize pressure field
    setPressure(job);

    //run simulation until stop_time
    while (job->t < stop_time){
        if (check_interupt()){
            //sigint called, exit
            return;
        }

        //run solver
        job->solver.solverStep(job);
        //std::cout << "Step Completed [" << ++stepCount << "]." << std::flush;
        if (job->serializer.serializerWriteFrame(job) == 1) {
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

void driverGenerateGravity(Job* job) {
    double g = -9.81;
    gravity.setZero();
    if (job->DIM == 1){
        gravity(0) = g;
    } else if (job->DIM == 2){
        gravity(1) = g;
    } else if (job->DIM == 3){
        gravity(2) = g;
    }
    return;
}

/*----------------------------------------------------------------------------*/

void driverApplyGravity(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (size_t i=0;i<job->bodies[b].points.b.rows();i++){
            job->bodies[b].points.b.row(i) = gravity.transpose();
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

std::string driverSaveState(Job* job, Serializer* serializer, std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename
    std::stringstream s;
    s << "mpm_v2.driver." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    if (ffile.is_open()){
        ffile << "# mpm_v2 drivers/default_driver.so\n";
        ffile << stop_time; //save stop time only
        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filepath+filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Driver Saved." << std::endl;

    return filename;
}

/*----------------------------------------------------------------------------*/

int driverLoadState(Job* job, Serializer* serializer, std::string fullpath){
    //job object should be loaded first
    std::string line;
    std::stringstream ss;
    std::ifstream fin(fullpath);

    if(fin.is_open()){
        std::getline(fin,line); //first line

        std::getline(fin,line); //stoptime
        ss = std::stringstream(line);
        ss >> stop_time;

        gravity = job->jobVector<double>(Job::ZERO);

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << " )." << std::endl;

        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Driver Loaded." << std::endl;
    return 1;
}

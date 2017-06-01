//
// Created by aaron on 5/28/17.
// default_driver.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

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

    //initialize gravity
    driverGenerateGravity(job);
    driverApplyGravity(job);

    //run simulation until stop_time
    while (job->t < stop_time){
        //run solver
        job->solver.solverStep(job);
        std::cout << "Step Completed [" << ++stepCount << "]." << std::flush;
        if (job->serializer.serializerWriteFrame(job) == 1) {
            //successful frame written
            std::cout << " Frame Written [" << ++frameCount << "]." << std::flush;
        }
        std::cout << "\r";
        job->t += job->dt;
    }
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

    //create filename
    std::ostringstream s;
    s << "mpm_v2.driver." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    std::string filename = s.str();
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

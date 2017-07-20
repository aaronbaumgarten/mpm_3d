//
// Created by aaron on 7/16/17.
// user_def_grav_driver.cpp
//

//
// Created by aaron on 7/15/17.
// shear_driver.cpp
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

void driverInit(Job* job){
    if (job->driver.fp64_props.size() < (1+job->DIM)) {
        std::cout << job->driver.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i properties defined ({stop_time, <gravity vector>}).\n",
                __FILE__, __func__, job->DIM);
        exit(0);
    } else {
        //store stop_time
        stop_time = job->driver.fp64_props[0];
        gravity = job->jobVector<double>(Job::ZERO);
        for (size_t pos=0;pos<job->DIM;pos++){
            gravity(pos) = job->driver.fp64_props[1+pos];
        }

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ", gravity = " << gravity.transpose() << ")." << std::endl;
    }

    std::cout << "Driver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void driverRun(Job* job) {
    size_t stepCount = 0;
    size_t frameCount = 0;
    clock_t clockSim = clock();
    clock_t clockFrame = clock();
    double tSim = 0;
    double tFrame = 0;

    //initialize gravity
    driverGenerateGravity(job);
    driverApplyGravity(job);

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
            tFrame = (double)(clock() - clockFrame)/CLOCKS_PER_SEC;
            tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
            printf("\33[2K");
            std::cout << "Frame Written [" << ++frameCount << "]. Time/Frame [" << tFrame << " s]. Elapsed Time [" << tSim << " s]." << std::flush;
            clockFrame = clock();
        }
        std::cout << "\r";

        job->t += job->dt;
    }
    tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
    std::cout << std::endl << std::endl << "Simulation Complete. Elapsed Time [" << tSim << "s]." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void driverGenerateGravity(Job* job) {
    //do nothing
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
        ffile << stop_time << "\n"; //save stop time only
        for (size_t i=0;i<gravity.rows();i++){
            ffile << gravity(i) << "\n";
        }
        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ", gravity = " << gravity.transpose() << ")." << std::endl;
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

        for (size_t i=0;i<gravity.rows();i++){
            std::getline(fin,line);
            gravity(i) = std::stod(line);
        }

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

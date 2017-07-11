//
// Created by aaron on 6/15/17.
// ballistic_impactor.cpp
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

int impactorID = -1;
double stop_time;
Eigen::VectorXd gravity;

double settling_time;
double settling_dt;

bool simStart;
double sim_dt;

extern "C" void driverInit(Job* job); //initialize driver
extern "C" void driverRun(Job* job); //run simulation
extern "C" void driverGenerateGravity(Job* job); //apply gravity to simulation
extern "C" void driverApplyGravity(Job* job); //apply gravity to simulation

extern "C" std::string driverSaveState(Job* job, Serializer* serializer, std::string filepath); //save driver state to returned filename in serializer folder
extern "C" int driverLoadState(Job* job, Serializer* serializer, std::string fullpath); //load state from given full path

/*----------------------------------------------------------------------------*/

void driverInit(Job* job){
    //set id of impacting body
    impactorID = -1;

    if (job->driver.fp64_props.size() < 3 || (job->driver.str_props.size() < 1 && job->driver.int_props.size() < 1)) {
        std::cout << job->driver.fp64_props.size() << ", " << job->driver.int_props.size() << ", " << job->driver.str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 4 property defined ({stop_time, settling_time, settling_dt}, {impactor}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        stop_time = job->driver.fp64_props[0];
        settling_time = job->driver.fp64_props[1];
        settling_dt = job->driver.fp64_props[2];
        gravity = job->jobVector<double>(Job::ZERO);

        //set body ids by name
        if (job->driver.str_props.size() >= 1) {
            for (size_t b = 0; b < job->bodies.size(); b++) {
                if (job->driver.str_props[0].compare(job->bodies[b].name) == 0) {
                    impactorID = b;
                    break;
                }
            }
        }

        // or set body ids by int
        if (impactorID < 0){
            if (job->driver.int_props.size() >= 1) {
                impactorID = job->driver.int_props[0];
            } else {
                std::cout << job->driver.fp64_props.size() << ", " << job->driver.int_props.size() << ", " << job->driver.str_props.size() << "\n";
                fprintf(stderr,
                        "%s:%s: Need at least 4 property defined ({stop_time, settling_time, settling_dt}, {impactor}).\n",
                        __FILE__, __func__);
                exit(0);
            }
        }

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ", settling_time = " << settling_time << ", settling_dt = " << settling_dt << ", impactor = ? )." << std::endl;
    }

    simStart = false;
    sim_dt = job->dt;

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

    //initialize simulation details
    if (!simStart) {
        job->t = -settling_time;
        job->dt = settling_dt;
    }

    std::cout << "Initializing Implicit Solver." << std::endl;

    //set solver to implicit for settling materials
    job->solver = Solver();
    job->solver.filepath = "src/solvers/";
    job->solver.filename = "newton_bicgstab.so";
    job->solver.fp64_props = {1e-10, 1e-7};
    job->solver.int_props = {100};
    job->solver.str_props = {};
    job->solver.solverSetPlugin(job,
                                job->serializer.mainpath + job->solver.filepath,
                                job->solver.filename,
                                job->solver.fp64_props,
                                job->solver.int_props,
                                job->solver.str_props);
    job->solver.solverInit(job);

    std::cout << std::endl << "Settling Phase:" << std::endl;

    //set impactor to inactive
    job->activeBodies[impactorID] = 0;

    //run simulation to settle material (w/o impactor)
    while (job->t < 0) {
        if (check_interupt()){
            return;
        }
        //run solver
        job->solver.solverStep(job);

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

    std::cout << std::endl << std::endl << "Initializing Explicit Solver." << std::endl;

    //set solver to explicit for impact simulation
    job->solver = Solver();
    job->solver.filepath = "src/solvers/";
    job->solver.filename = "explicit_usl.so";
    job->solver.fp64_props = {};
    job->solver.int_props = {};
    job->solver.str_props = {};
    job->solver.solverSetPlugin(job,
                                job->serializer.mainpath + job->solver.filepath,
                                job->solver.filename,
                                job->solver.fp64_props,
                                job->solver.int_props,
                                job->solver.str_props);
    job->solver.solverInit(job);

    std::cout << std::endl << "Ballistic Phase:" << std::endl;

    //set impactor to active and fix dt etc.
    job->activeBodies[impactorID] = 1;
    job->dt = sim_dt;

    //run simulation until stop_time
    while (job->t < stop_time){
        if (check_interupt()){
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
        ffile << stop_time << "\n"; //save stop
        ffile << impactorID << "\n"; //save impactor id
        ffile << settling_time << "\n";
        ffile << settling_dt << "\n";
        ffile << (int)simStart << "\n"; //integer of boolean
        ffile << sim_dt << "\n";
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
    std::ifstream fin(fullpath);

    if(fin.is_open()){
        std::getline(fin,line); //first line

        std::getline(fin,line); //stoptime
        stop_time = std::stod(line);
        std::getline(fin,line); //impactorID
        impactorID = std::stoi(line);
        std::getline(fin,line); //settline time
        settling_time = std::stod(line);
        std::getline(fin,line); //settling dt
        settling_dt = std::stod(line);
        std::getline(fin,line); //simStart
        simStart = (bool)std::stoi(line);
        std::getline(fin,line); //sim_dt
        sim_dt = std::stod(line);

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

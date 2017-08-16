//
// Created by aaron on 8/14/17.
// launched_projectile_driver.cpp
//

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

std::string name;
int impactorID = -1;
double stop_time, launch_speed, launch_angle, sample_rate;
Eigen::VectorXd gravity;

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

    if (job->DIM < 2){
        std::cerr << "ERROR: launched_projectile_driver.so requires at least 2 dimensions!" << std::endl;
        exit(0);
    }

    if (job->driver.fp64_props.size() < 4 || (job->driver.str_props.size() < 2 && job->driver.int_props.size() < 1)) {
        std::cout << job->driver.fp64_props.size() << ", " << job->driver.int_props.size() << ", " << job->driver.str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 6 property defined ({stop_time, launch_speed, launch_angle, sample_rate}, {name, impactor}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        stop_time = job->driver.fp64_props[0];
        launch_speed = job->driver.fp64_props[1];
        launch_angle = job->driver.fp64_props[2];
        sample_rate = job->driver.fp64_props[3];
        gravity = job->jobVector<double>(Job::ZERO);
        name = job->driver.str_props[0];

        //set body ids by name
        if (job->driver.str_props.size() >= 2) {
            for (size_t b = 0; b < job->bodies.size(); b++) {
                if (job->driver.str_props[1].compare(job->bodies[b].name) == 0) {
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
                        "%s:%s: Need at least 6 property defined ({stop_time, launch_speed, launch_angle, sample_rate}, {name, impactor}).\n",
                        __FILE__, __func__);
                exit(0);
            }
        }

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ", launch_speed = " << launch_speed << ", launch_angle = " << launch_angle << ", sample_rate = " << sample_rate << ", name = " << name << ", impactor = " << job->bodies[impactorID].name << " )." << std::endl;
    }

    std::cout << "Driver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void driverRun(Job* job) {
    size_t stepCount = 0;
    size_t frameCount = 0;
    size_t sampleCount = 0;
    clock_t clockSim = clock();
    clock_t clockFrame = clock();
    double tSim = 0;
    double tFrame = 0;
    double pi = std::acos(-1.0);

    //initialize gravity
    driverGenerateGravity(job);
    driverApplyGravity(job);

    //initialize launch of projectile
    job->bodies[impactorID].points.x_t.setZero();
    job->bodies[impactorID].points.mx_t.setZero();
    for (size_t i=0;i<job->bodies[impactorID].points.x_t.rows();i++) {
        job->bodies[impactorID].points.x_t(i,0) = launch_speed*std::cos(launch_angle*pi/180.0);
        job->bodies[impactorID].points.x_t(i,job->DIM-1) = launch_speed*std::sin(launch_angle*pi/180.0);
        job->bodies[impactorID].points.mx_t(i,0) = job->bodies[impactorID].points.m(i) * job->bodies[impactorID].points.x_t(i,0);
        job->bodies[impactorID].points.mx_t(i,job->DIM-1) = job->bodies[impactorID].points.m(i) * job->bodies[impactorID].points.x_t(i,job->DIM-1);
    }

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

        if (sampleCount <= job->t*sample_rate) {
            sampleCount += 1;
            std::ofstream ffile(name + ".csv", std::ios::app);
            if (ffile.is_open()) {
                ffile << job->t;
                Eigen::VectorXd x_t = job->jobVector<double>();
                Eigen::VectorXd x = job->jobVector<double>();
                double m;
                for (size_t b = 0; b < job->bodies.size(); b++) {
                    if (job->activeBodies[b] == 0) {
                        continue;
                    }
                    m = 0;
                    x_t.setZero();
                    x.setZero();
                    for (size_t i = 0; i < job->bodies[b].points.m.rows(); i++) {
                        x_t = x_t + job->bodies[b].points.m(i)*job->bodies[b].points.x_t.row(i).transpose();
                        x = x + job->bodies[b].points.m(i)*job->bodies[b].points.x.row(i).transpose();
                        m += job->bodies[b].points.m(i);
                    }
                    x_t = x_t / m;
                    x = x / m;

                    for (size_t pos=0;pos<x.size();pos++){
                        ffile << ", " << x(pos);
                    }

                    for (size_t pos=0;pos<x_t.size();pos++){
                        ffile << ", " << x_t(pos);
                    }
                }
                ffile << "\n";
            }
            ffile.close();
        }

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
        ffile << launch_speed << "\n";
        ffile << launch_angle << "\n";
        ffile << sample_rate << "\n";
        ffile << name << "\n";
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
        launch_speed = std::stod(line);
        std::getline(fin,line); //settling dt
        launch_angle = std::stod(line);
        std::getline(fin,line); //sample_rate
        sample_rate = std::stod(line);
        std::getline(fin,line); //name
        name = line;

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

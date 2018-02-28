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
std::string output_name;
Eigen::VectorXd gravity;

extern "C" void driverInit(Job* job); //initialize driver
extern "C" void driverRun(Job* job); //run simulation
extern "C" void driverGenerateGravity(Job* job); //apply gravity to simulation
extern "C" void driverApplyGravity(Job* job); //apply gravity to simulation

extern "C" std::string driverSaveState(Job* job, Serializer* serializer, std::string filepath); //save driver state to returned filename in serializer folder
extern "C" int driverLoadState(Job* job, Serializer* serializer, std::string fullpath); //load state from given full path

/*----------------------------------------------------------------------------*/

void driverInit(Job* job){
    if (job->driver.fp64_props.size() < 1 || job->driver.str_props.size() < 1) {
        std::cout << job->driver.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 property defined ({stop_time},{name}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        stop_time = job->driver.fp64_props[0];
        gravity = job->jobVector<double>(Job::ZERO);
        output_name = job->driver.str_props[0];

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ", name = " + output_name + ")." << std::endl;
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

        //switch to 3D
        job->DIM = 3;

        std::ofstream ffile(output_name+".csv", std::ios::app);
        if (ffile.is_open()) {
            ffile << job->t;
            double p, tau, v;
            Eigen::MatrixXd T = job->jobTensor<double>();
            Eigen::MatrixXd T_avg = job->jobTensor<double>();
            Eigen::VectorXd tmpVec;
            for (size_t b = 0; b<job->bodies.size(); b++) {
                if (job->activeBodies[b] == 0){
                    continue;
                }
                T_avg.setZero();
                tau = 0;
                p = 0;
                v = 0;
                for (size_t i = 0; i < job->bodies[b].points.m.rows(); i++) {
                    tmpVec = job->bodies[b].points.T.row(i);
                    T = job->jobTensor<double>(tmpVec.data());
                    v += job->bodies[b].points.v(i);
                    //p -= T.trace() / T.rows() * job->bodies[b].points.v(i);
                    //tau += (T - T.trace()/T.rows()*job->jobTensor<double>(Job::IDENTITY)).norm() * job->bodies[b].points.v(i);
                    T_avg += T*job->bodies[b].points.v(i);
                }
                //p /= v;
                //tau /= v*std::sqrt(2.0);
                T_avg  = T_avg / v;
                tau = (T_avg - T_avg.trace()/T_avg.rows()*job->jobTensor<double>(Job::IDENTITY)).norm() / std::sqrt(2);
                p = -T_avg.trace() / T_avg.rows();

                ffile << ", " << p << ", " << tau << ", " << v;
            }
            ffile << "\n";
        }
        ffile.close();

        //switch to 2D
        job->DIM = 2;

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
    //if (job->DIM == 1){
    //    gravity(0) = g;
    //} else if (job->DIM == 2){
    //    gravity(1) = g;
    //} else if (job->DIM == 3){
    //    gravity(2) = g;
    //}
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
        ffile << "# mpm_v2 drivers/shear_driver.so\n";
        ffile << stop_time << "\n"; //save stop time only
        ffile << output_name;
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

        std::getline(fin,line);
        output_name = line;

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

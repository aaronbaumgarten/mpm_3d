//
// Created by aaron on 10/12/17.
// shear_traction_driver_subsample.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <time.h>
#include <ctime>
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
int traction_body_id;
double f1;
double f2;
double sample_rate;
int sampled_frames = 0;
std::vector<int> sample_ids(0);

extern "C" void driverInit(Job* job); //initialize driver
extern "C" void driverRun(Job* job); //run simulation
extern "C" void driverGenerateGravity(Job* job); //apply gravity to simulation
extern "C" void driverApplyGravity(Job* job); //apply gravity to simulation

extern "C" std::string driverSaveState(Job* job, Serializer* serializer, std::string filepath); //save driver state to returned filename in serializer folder
extern "C" int driverLoadState(Job* job, Serializer* serializer, std::string fullpath); //load state from given full path

/*----------------------------------------------------------------------------*/

void driverInit(Job* job){
    if (job->DIM != 2){
        std::cerr << "ERROR: shear_traction_driver.so requires DIM = 2." << std::endl;
        exit(0);
    }

    if (job->driver.fp64_props.size() < 4 || job->driver.int_props.size() < 1 || job->driver.str_props.size() < 2) {
        std::cout << job->driver.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 7 property defined ({stop_time, f1, f2, sample_rate},{<list of point ids to sample>},{name,traction_body}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        stop_time = job->driver.fp64_props[0];
        f1 = job->driver.fp64_props[1];
        f2 = job->driver.fp64_props[2];
        sample_rate = job->driver.fp64_props[3];

        sample_ids = job->driver.int_props;

        gravity = job->jobVector<double>(Job::ZERO);
        output_name = job->driver.str_props[0];
        for (size_t b = 0; b < job->bodies.size(); b++) {
            if (job->driver.str_props[1].compare(job->bodies[b].name) == 0){
                traction_body_id = b;
                break;
            }
        }

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ", body force = " << f1 << ", " << f2 << ", name = " + output_name + ")." << std::endl;
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

        if (job->t*sample_rate >= sampled_frames) {
            std::ofstream ffile(output_name + ".csv", std::ios::app);
            if (ffile.is_open()) {
                ffile << job->t;
                double p, tau, v, m, gdp;
                Eigen::MatrixXd T = job->jobTensor<double>();
                Eigen::MatrixXd T_avg = job->jobTensor<double>();
                Eigen::MatrixXd L = job->jobTensor<double>();
                Eigen::MatrixXd D = job->jobTensor<double>();
                Eigen::MatrixXd D_avg = job->jobTensor<double>();
                Eigen::VectorXd tmpVec;
                for (size_t b = 0; b < job->bodies.size(); b++) {
                    if (job->activeBodies[b] == 0) {
                        continue;
                    }
                    T_avg.setZero();
                    D_avg.setZero();
                    tau = 0;
                    p = 0;
                    v = 0;
                    m = 0;
                    gdp = 0;
                    for (size_t i = 0; i < sample_ids.size(); i++) {
                        if (sample_ids[i] >= job->bodies[b].points.x.rows()){
                            continue;
                        }
                        tmpVec = job->bodies[b].points.T.row(sample_ids[i]);
                        T = job->jobTensor<double>(tmpVec.data());
                        tmpVec = job->bodies[b].points.L.row(sample_ids[i]);
                        L = job->jobTensor<double>(tmpVec.data());
                        D = 0.5 * (L + L.transpose());
                        v += job->bodies[b].points.v(sample_ids[i]);
                        m += job->bodies[b].points.m(sample_ids[i]);
                        //p -= T.trace() / T.rows() * job->bodies[b].points.v(i);
                        //tau += (T - T.trace()/T.rows()*job->jobTensor<double>(Job::IDENTITY)).norm() * job->bodies[b].points.v(i);
                        T_avg += T * job->bodies[b].points.v(sample_ids[i]);
                        D_avg += D * job->bodies[b].points.v(sample_ids[i]);
                    }
                    //p /= v;
                    //tau /= v*std::sqrt(2.0);
                    T_avg = T_avg / v;
                    D_avg = D_avg / v;
                    tau = (T_avg - T_avg.trace() / T_avg.rows() * job->jobTensor<double>(Job::IDENTITY)).norm() /
                          std::sqrt(2.0);
                    p = -T_avg.trace() / T_avg.rows();
                    gdp = (D_avg - D_avg.trace() / D_avg.rows() * job->jobTensor<double>(Job::IDENTITY)).norm() *
                          std::sqrt(2.0);

                    ffile << ", " << p << ", " << tau << ", " << v << ", " << m << ", " << gdp;
                }
                ffile << "\n";
            }
            ffile.close();
            sampled_frames += 1;
        }

        job->t += job->dt;
    }
    tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
    std::cout << std::endl << std::endl << "Simulation Complete. Elapsed Time [" << tSim << "s]." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void driverGenerateGravity(Job* job) {
    double m_tot = job->bodies[traction_body_id].points.m.sum();
    gravity(0) = f1/m_tot;
    gravity(1) = f2/m_tot;
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
        job->bodies[b].points.b.setZero();
    }

    for (size_t i=0;i<job->bodies[traction_body_id].points.b.rows();i++){
        job->bodies[traction_body_id].points.b.row(i) = gravity.transpose();
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
        ffile << f1 << "\n" << f2 << "\n" << sample_rate << "\n" << sampled_frames << "\n" << traction_body_id << "\n";
        ffile << output_name << "\n";
        ffile << sample_ids.size();
        for (size_t i=0;i<sample_ids.size();i++){
            ffile << "\n" << sample_ids[i];
        }
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
    int len;

    if(fin.is_open()){
        std::getline(fin,line); //first line

        std::getline(fin,line); //stoptime
        ss = std::stringstream(line);
        ss >> stop_time;

        std::getline(fin,line);
        f1 = std::stod(line);

        std::getline(fin,line);
        f2 = std::stod(line);

        std::getline(fin,line);
        sample_rate = std::stod(line);

        std::getline(fin,line);
        sampled_frames = std::stoi(line);

        std::getline(fin,line);
        traction_body_id = std::stoi(line);

        std::getline(fin,line);
        output_name = line;

        std::getline(fin,line);
        len = std::stoi(line);
        sample_ids.resize(len);

        for (size_t i=0;i<len;i++){
            std::getline(fin,line);
            sample_ids[i] = std::stoi(line);
        }

        gravity = job->jobVector<double>(Job::ZERO);

        //print grid properties
        std::cout << "Driver properties (stop_time = " << stop_time << ", body force = " << f1 << ", " << f2 << ", name = " + output_name + ")." << std::endl;

        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Driver Loaded." << std::endl;
    return 1;
}

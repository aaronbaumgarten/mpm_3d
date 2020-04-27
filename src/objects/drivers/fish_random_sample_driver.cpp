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
#include <ctime>
#include <ostream>

#include "job.hpp"
#include "mpm_objects.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"

#include "drivers.hpp"

/*----------------------------------------------------------------------------*/

void FishRandomSampleDriver::init(Job* job){
    if (fp64_props.size() < 7 || int_props.size() < 1 || (str_props.size() < 2 && int_props.size() < 3)) {
        std::cout << fp64_props.size() << ", " << int_props.size() << ", " << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 11 properties defined ({stop_time, alpha_min, alpha_max, beta_min, beta_max, gamma_min, gamma_max}, {num_samples}, {filename, fish_body}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        stop_time = fp64_props[0];
        gravity = KinematicVector(job->JOB_TYPE);

        alpha_min = fp64_props[1];
        alpha_max = fp64_props[2];
        beta_min = fp64_props[3];
        beta_max = fp64_props[4];
        gamma_min = fp64_props[5];
        gamma_max = fp64_props[6];

        num_samples = int_props[0];

        output_file = str_props[0];

        //set body ids by name
        fish_body = -1;
        if (str_props.size() >= 2){
            for (int b = 0; b < job->bodies.size(); b++) {
                if (str_props[1].compare(job->bodies[b]->name) == 0){
                    fish_body = b;
                    break;
                }
            }
        }

        // or set body ids by int
        if (fish_body < 0){
            if (int_props.size() >= 2) {
                fish_body = int_props[1];
            } else {
                std::cout << fp64_props.size() << ", " << int_props.size() << ", " << str_props.size() << "\n";
                fprintf(stderr,
                        "%s:%s: Need at least 10 properties defined ({stop_time, alpha_min, alpha_max, beta_min, beta_max, gamma_min, gamma_max}, {num_samples}, {fish_body}).\n",
                        __FILE__, __func__);
                exit(0);
            }
        }

        //print grid properties
        std::cout << "Driver properties:" << std::endl;
        std::cout << "    filename   : " << output_file << std::endl;
        std::cout << "    stop_time  : " << stop_time << std::endl;
        std::cout << "    num_samples: " << num_samples << std::endl;
        std::cout << "    fish_body  : " << fish_body << std::endl;
        std::cout << "    alpha_min  : " << alpha_min << std::endl;
        std::cout << "    alpha_max  : " << alpha_max << std::endl;
        std::cout << "    beta_min   : " << beta_min << std::endl;
        std::cout << "    beta_max   : " << beta_max << std::endl;
        std::cout << "    gamma_min  : " << gamma_min << std::endl;
        std::cout << "    gamma_max  : " << gamma_max << std::endl;

    }

    std::cout << "Driver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void FishRandomSampleDriver::run(Job* job) {
    //open output file
    std::ofstream file;
    file.open(output_file, std::ios::trunc);
    if (!file.is_open()){
        std::cerr << "ERROR! Cannot open output file " << output_file << "! Exiting." << std::endl;
        exit(0);
    } else {
        file.close();
    }

    int frameCount = 0;
    int sampleCount = 0;

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

    //seed randomizer
    std::srand(std::time(0));

    //generate random samples
    while (sampleCount < num_samples) {
        //iterate counter
        sampleCount++;

        //reset job
        job->t = 0;

        //random sample of fish properties
        std::vector<double> fish_fp64_props(11);
        fish_fp64_props[0] = job->bodies[fish_body]->material->fp64_props[0];
        fish_fp64_props[1] = job->bodies[fish_body]->material->fp64_props[1];
        fish_fp64_props[2] = alpha_min + (double)std::rand()/RAND_MAX * (alpha_max - alpha_min);
        fish_fp64_props[3] = alpha_min + (double)std::rand()/RAND_MAX * (alpha_max - alpha_min);
        fish_fp64_props[4] = alpha_min + (double)std::rand()/RAND_MAX * (alpha_max - alpha_min);
        fish_fp64_props[5] = beta_min + (double)std::rand()/RAND_MAX * (beta_max - beta_min);
        fish_fp64_props[6] = beta_min + (double)std::rand()/RAND_MAX * (beta_max - beta_min);
        fish_fp64_props[7] = beta_min + (double)std::rand()/RAND_MAX * (beta_max - beta_min);
        fish_fp64_props[8] = gamma_min + (double)std::rand()/RAND_MAX * (gamma_max - gamma_min);
        fish_fp64_props[9] = gamma_min + (double)std::rand()/RAND_MAX * (gamma_max - gamma_min);
        fish_fp64_props[10] = gamma_min + (double)std::rand()/RAND_MAX * (gamma_max - gamma_min);

        //assign properties to fish
        job->bodies[fish_body]->material->fp64_props = fish_fp64_props;

        //re-initialize bodies
        for (int id=0; id<job->bodies.size(); id++) {
            job->bodies[id]->points->readFromFile(job, job->bodies[id].get(), job->bodies[id]->points->file);
            job->bodies[id]->init(job);
        }

        //run simulation until stop_time
        int stepCount = 0;
        frameCount = 0;
        std::cout << "Simulation Number: " << sampleCount << "/" << num_samples << std::endl;
        while (job->t <= stop_time) {
            //run solver
            job->solver->step(job);
            stepCount++;
            if (stepCount % 100 == 0) {
                clock_gettime(CLOCK_MONOTONIC, &timeFinish);
                tFrame = (timeFinish.tv_sec - timeFrame.tv_sec) +
                         (timeFinish.tv_nsec - timeFrame.tv_nsec) / 1000000000.0;
                tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec) / 1000000000.0;
                timeFrame = timeFinish;
                printf("\33[2K");
                std::cout << "Steps Complete [" << stepCount << "/" << (int)(stop_time/job->dt) << "]. Elapsed Time ["
                          << tSim << " s]." << std::flush;

                std::cout << "\r";
            }
            job->t += job->dt;
        }

        //wrap up
        std::cout << "\n" << std::flush;

        //determine success metric
        double m = 0;
        KinematicVector cm = KinematicVector(job->JOB_TYPE);
        for (int ii=0; ii<job->bodies[fish_body]->points->x.size(); ii++){
            m += job->bodies[fish_body]->points->m(ii);
            cm += job->bodies[fish_body]->points->m(ii) * job->bodies[fish_body]->points->u[ii];
        }

        //open file
        file.open(output_file,std::ios::app);

        if (file.is_open()) {
            //write metrics to file
            file << sampleCount << ", ";
            file << fish_fp64_props[2] << ", ";
            file << fish_fp64_props[3] << ", ";
            file << fish_fp64_props[4] << ", ";
            file << fish_fp64_props[5] << ", ";
            file << fish_fp64_props[6] << ", ";
            file << fish_fp64_props[7] << ", ";
            file << fish_fp64_props[8] << ", ";
            file << fish_fp64_props[9] << ", ";
            file << fish_fp64_props[10] << ", ";
            file << cm(0) / m << ", ";
            file << cm(1) / m << "\n";

            //close file
            file.close();
        } else {
            std::cout << "ERROR! Failed to write to file!" << std::endl;
        }

    }

    //close file
    file.close();

    //write out final simulation time
    clock_gettime(CLOCK_MONOTONIC,&timeFinish);
    tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
    std::cout << std::endl << std::endl << "Simulation Complete. Elapsed Time [" << tSim << "s]." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void FishRandomSampleDriver::generateGravity(Job* job) {
    gravity.setZero();
    return;
}

/*----------------------------------------------------------------------------*/

void FishRandomSampleDriver::applyGravity(Job* job){
    //do nothing
    return;
}

/*----------------------------------------------------------------------------*/

std::string FishRandomSampleDriver::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "error";
}

/*----------------------------------------------------------------------------*/

int FishRandomSampleDriver::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

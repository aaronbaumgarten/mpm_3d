//
// Created by aaron on 5/8/20.
// fish_mixture_model_random_sample_driver.cpp.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <time.h>
#include <ctime>
#include <ostream>
#include <istream>

#include "job.hpp"
#include "mpm_objects.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "parser.hpp"

#include "drivers.hpp"

/*----------------------------------------------------------------------------*/

void FishMixtureModelRandomSampleDriver::init(Job* job){
    if (fp64_props.size() < 7 || int_props.size() < 1 || (str_props.size() < 2 && int_props.size() < 3)) {
        std::cout << fp64_props.size() << ", " << int_props.size() << ", " << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 12 properties defined ({stop_time, alpha_min, alpha_max, beta_min, beta_max, gamma_min, gamma_max}, {num_samples}, {output_file, input_file, fish_body}).\n",
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
        input_file = str_props[1];

        //set body ids by name
        fish_body = -1;
        if (str_props.size() >= 3){
            for (int b = 0; b < job->bodies.size(); b++) {
                if (str_props[2].compare(job->bodies[b]->name) == 0){
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
                        "%s:%s: Need at least 12 properties defined ({stop_time, alpha_min, alpha_max, beta_min, beta_max, gamma_min, gamma_max}, {num_samples}, {output_file, input_file, fish_body}).\n",
                        __FILE__, __func__);
                exit(0);
            }
        }

        //print grid properties
        std::cout << "Driver properties:" << std::endl;
        std::cout << "    output_file: " << output_file << std::endl;
        std::cout << "    input_file : " << input_file << std::endl;
        std::cout << "    stop_time  : " << stop_time << std::endl;
        std::cout << "    num_samples: " << num_samples << std::endl;
        std::cout << "    fish_body  : " << fish_body << std::endl;
        std::cout << "    alpha_min  : " << alpha_min << std::endl;
        std::cout << "    alpha_max  : " << alpha_max << std::endl;
        std::cout << "    beta_min   : " << beta_min << std::endl;
        std::cout << "    beta_max   : " << beta_max << std::endl;
        std::cout << "    gamma_min  : " << gamma_min << std::endl;
        std::cout << "    gamma_max  : " << gamma_max << std::endl;

        //initialize normal distributions
        normal_distributions = std::vector<MultivariateNormalRandomVariable>();
        probability_distribution = std::vector<double>();

        //code to load random number generators from input file
        int k = 0;
        int d = 0;
        bool fail = false;
        std::string line;
        std::vector<std::string> lvec;
        std::ifstream fin(input_file);
        if (fin.is_open()){
            std::getline(fin,line);
            line = Parser::removeSpaces(line);
            lvec = Parser::splitString(line,',');
            //read number of clusters and number of dimensions (should be 9)
            if (lvec.size() > 1) {
                k = std::stoi(lvec[0]);
                d = std::stoi(lvec[1]);
            } else {
                fail = true;
            }

            //read header line that follows
            std::getline(fin,line);

            Eigen::VectorXd mu = Eigen::VectorXd(d);
            Eigen::MatrixXd Sigma = Eigen::MatrixXd(d,d);

            for (int j = 0; j < k; j++){
                //add distribution to list
                probability_distribution.push_back(0.0);
                normal_distributions.push_back(MultivariateNormalRandomVariable());

                //read line from file
                std::getline(fin,line);
                //remove spaces
                line = Parser::removeSpaces(line);
                //break line into parts
                lvec = Parser::splitString(line,',');

                //check lvec length
                if (lvec.size() < 1 + d + d*d){
                    fail = true;
                }

                //assign p
                probability_distribution[j] = std::stod(lvec[0]);

                //assign mu and Sigma from line info
                for (int ii = 0; ii < d; ii++){
                    mu(ii) = std::stod(lvec[1 + ii]);
                    for (int jj = 0; jj < d; jj++){
                        Sigma(ii,jj) = std::stod(lvec[1 + d + ii*d + jj]);
                    }
                }

                //update normal distribution
                normal_distributions[j].setMean(mu);
                normal_distributions[j].setCovar(Sigma);
            }

            fin.close();
        } else {
            fail = true;
        }

        if (fail){
            std::cout << "ERROR: Unable to open file: " << input_file << std::endl;
            exit(0);
        }

    }

    std::cout << "Driver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void FishMixtureModelRandomSampleDriver::run(Job* job) {
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
        /*
        fish_fp64_props[2] = alpha_min + (double)std::rand()/RAND_MAX * (alpha_max - alpha_min);
        fish_fp64_props[3] = alpha_min + (double)std::rand()/RAND_MAX * (alpha_max - alpha_min);
        fish_fp64_props[4] = alpha_min + (double)std::rand()/RAND_MAX * (alpha_max - alpha_min);
        fish_fp64_props[5] = beta_min + (double)std::rand()/RAND_MAX * (beta_max - beta_min);
        fish_fp64_props[6] = beta_min + (double)std::rand()/RAND_MAX * (beta_max - beta_min);
        fish_fp64_props[7] = beta_min + (double)std::rand()/RAND_MAX * (beta_max - beta_min);
        fish_fp64_props[8] = gamma_min + (double)std::rand()/RAND_MAX * (gamma_max - gamma_min);
        fish_fp64_props[9] = gamma_min + (double)std::rand()/RAND_MAX * (gamma_max - gamma_min);
        fish_fp64_props[10] = gamma_min + (double)std::rand()/RAND_MAX * (gamma_max - gamma_min);
         */
        double r_value = (double)std::rand()/RAND_MAX; //random number between 0 and 1
        Eigen::VectorXd r_vector = Eigen::VectorXd::Zero(9);
        //determine which distribution to select from
        for (int j=0; j<probability_distribution.size(); j++){
            if (r_value < probability_distribution[j] || j+1 == probability_distribution.size()){
                //use jth normal distribution
                std::cout << "[" << sampleCount << "]: Picked " << j << "th distribution." << std::endl;
                //get vector of values from distribution
                r_vector = normal_distributions[j]();
            }
        }
        fish_fp64_props[2] = r_vector(0);
        fish_fp64_props[3] = r_vector(1);
        fish_fp64_props[4] = r_vector(2);
        fish_fp64_props[5] = r_vector(3);
        fish_fp64_props[6] = r_vector(4);
        fish_fp64_props[7] = r_vector(5);
        fish_fp64_props[8] = r_vector(6);
        fish_fp64_props[9] = r_vector(7);
        fish_fp64_props[10] = r_vector(8);


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

void FishMixtureModelRandomSampleDriver::generateGravity(Job* job) {
    gravity.setZero();
    return;
}

/*----------------------------------------------------------------------------*/

void FishMixtureModelRandomSampleDriver::applyGravity(Job* job){
    //do nothing
    return;
}

/*----------------------------------------------------------------------------*/

std::string FishMixtureModelRandomSampleDriver::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "error";
}

/*----------------------------------------------------------------------------*/

int FishMixtureModelRandomSampleDriver::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

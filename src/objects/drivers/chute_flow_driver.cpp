//
// Created by aaron on 10/26/18.
// ballistic_driver.cpp
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

#include "drivers.hpp"

/*----------------------------------------------------------------------------*/

void ChuteFlowDriver::init(Job* job){
    if (fp64_props.size() < 10 || str_props.size() < 1) {
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 11 properties defined (stop_time, g, mu_1, mu_2, b, d, rho_s, h, theta, phi) and filename.\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        stop_time = fp64_props[0];
        gravity = KinematicVector(job->JOB_TYPE);
        gravity.setZero();

        //store output filename
        output_filename = str_props[0];

        //store chute flow variables
        g = fp64_props[1];
        mu_1 = fp64_props[2];
        mu_2 = fp64_props[3];
        b = fp64_props[4];
        d = fp64_props[5];
        rho_s = fp64_props[6];
        h = fp64_props[7];
        theta = fp64_props[8];
        phi = fp64_props[9];

        if (job->JOB_TYPE != 2){
            std::cerr << "ChuteFlowDriver implemented for 2D sand flows only! Exiting." << std::endl;
            exit(0);
        }

        //assign gravity
        gravity[0] = g*std::sin(M_PI*theta/180.0);
        gravity[1] = -g*std::cos(M_PI*theta/180.0);

        //open and start output file
        std::ofstream file (output_filename,std::ios::trunc);
        if (file.is_open()){
            //success!
            //write file header
            file << "t, ||e||_\\infty, ||v - v_ss||_L2, max(y), avg(y), max(v_x), avg(v_x)\n";

            file.close();
        } else {
            std::cerr << "ERROR! Cannot open " << output_filename << "! Exiting." << std::endl;
            exit(0);
        }

        //print grid properties
        std::cout << "ChuteFlowDriver properties:" << std::endl;
        std::cout << "    stop_time = " << stop_time << "\n";
        std::cout << "    gravity   = " << g << "\n";
        std::cout << "    theta     = " << theta << " deg\n";
        std::cout << "    mu_1      = " << mu_1 << "\n";
        std::cout << "    mu_2      = " << mu_2 << "\n";
        std::cout << "    b         = " << b << "\n";
        std::cout << "    d         = " << d << "\n";
        std::cout << "    rho_s     = " << rho_s << "\n";
        std::cout << "    h         = " << h << "\n";
        std::cout << "    phi       = " << phi << "\n";

    }

    std::cout << "Driver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void ChuteFlowDriver::generateGravity(Job* job) {
    //already assigned
    return;
}

/*----------------------------------------------------------------------------*/

void ChuteFlowDriver::applyGravity(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (int i=0;i<job->bodies[b]->points->b.size();i++){
            job->bodies[b]->points->b(i) = gravity;
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void ChuteFlowDriver::run(Job* job) {
    int stepCount = 0;
    int frameCount = 0;

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

    //assign lithostatic pressure
    if (job->bodies.size() > 1){
        std::cout << "Hmmmm... Too many bodies in simulation for ChuteFlowDriver..." << std::endl;
    }
    for (int i=0; i<job->bodies[0]->points->x.size(); i++){
        job->bodies[0]->points->T[i] = rho_s*phi*gravity[1]*(h - job->bodies[0]->points->x(i,1))*MaterialTensor::Identity();
    }

    //run simulation until stop_time
    while (job->t <= stop_time){
        //run solver
        job->solver->step(job);
        //std::cout << "Step Completed [" << ++stepCount << "]." << std::flush;
        if (job->serializer->writeFrame(job) == 1) {
            //successful frame written
            //write error measures to file
            writeErrorInfo(job);

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


KinematicVector ChuteFlowDriver::getVelocity(Job* job, KinematicVector const &x){
    KinematicVector result = KinematicVector(job->JOB_TYPE);
    result[1] = 0;
    result[0] = 2.0/3.0 * ((std::tan(M_PI*theta/180.0) - mu_1) / (mu_2 - std::tan(M_PI*theta/180.0)))
                * (b / (d * sqrt(rho_s))) * std::sqrt(-rho_s*phi*gravity[1])
                * (std::pow(h, 3.0/2.0) - (std::pow(h - x[1], 3.0/2.0)));
    return result;
}

void ChuteFlowDriver::writeErrorInfo(Job* job){
    //calculate and write error measures to output file for first body in simulation
    //BE CAREFUL!!!

    //initialize arrays
    Eigen::VectorXd V_i = Eigen::VectorXd(job->grid->node_count);
    Eigen::VectorXd v_i = Eigen::VectorXd(job->grid->node_count);
    Eigen::VectorXd e = Eigen::VectorXd(job->grid->node_count);

    //initialize figures of merit
    double H_norm = 0;
    double e_norm = 0;
    double v_L2 = 0;
    double y_max = 0;
    double y_avg = 0;
    double v_max = 0;
    double v_avg = 0;

    //get exact node volumes form grid
    for (int i=0; i<job->grid->node_count;i++){
        V_i(i) = job->grid->nodeVolume(job,i);
    }

    //get integrated node volume from points
    v_i = job->bodies[0]->S * job->bodies[0]->points->v;

    //calculate arrays
    double tmpNum;
    for (int i = 0; i < V_i.rows(); i++) {
        tmpNum = (v_i(i) - V_i(i));
        e(i) = std::max(0.0, tmpNum) / V_i(i);
    }

    //calculate error measures
    double m_sum = 0;
    KinematicVector tmpVel = KinematicVector(job->JOB_TYPE);
    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    for (int i=0; i<V_i.rows(); i++){
        if (job->bodies[0]->nodes->m(i) > 0) {
            //||e||_\infty
            if (e(i) > e_norm) {
                e_norm = e(i);
            }

            //||v - v_ss||_L2
            tmpVel = getVelocity(job, job->bodies[0]->nodes->x[i]);
            tmpVec = (tmpVel - job->bodies[0]->nodes->x_t[i]);
            //std::cout << "[" << i << "] " << tmpVel[0] << " ?= " << job->bodies[0]->nodes->x_t(i,0) << std::endl;
            v_L2 += tmpVec.dot(tmpVec) * V_i(i);

            m_sum += job->bodies[0]->nodes->m(i);
            v_avg += job->bodies[0]->nodes->m(i) * job->bodies[0]->nodes->x_t(i,0);
            if (job->bodies[0]->nodes->x_t(i,0) > v_max){
                v_max = job->bodies[0]->nodes->x_t(i,0);
            }
        }
    }
    //sqrt of ||a_err||_L2^2
    v_L2 = std::sqrt(v_L2);
    //v_avg = mx_t_sum / m
    v_avg /= m_sum;

    m_sum = 0;
    for (int p=0; p<job->bodies[0]->points->x.size(); p++){
        if (job->bodies[0]->points->active(p) != 0){
            m_sum += job->bodies[0]->points->m(p);
            y_avg += job->bodies[0]->points->m(p) * job->bodies[0]->points->x(p,1);

            if (job->bodies[0]->points->x(p,1) > y_max){
                y_max = job->bodies[0]->points->x(p,1);
            }
        }
    }
    y_avg /= m_sum;

    //open and write to file
    std::ofstream file (output_filename,std::ios::app);
    if (file.is_open()){
        //success!
        //write to file
        file << job->t << ", ";
        file << e_norm << ", ";
        file << v_L2 << ", ";
        file << y_max << ", ";
        file << y_avg << ", ";
        file << v_max << ", ";
        file << v_avg << "\n";

        file.close();
    } else {
        std::cerr << "ERROR! Cannot open " << output_filename << "!" << std::endl;;
    }
    return;
}
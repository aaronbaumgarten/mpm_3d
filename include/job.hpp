//
// Created by aaron on 5/10/18.
// job.hpp
//

#ifndef MPM_V3_JOB_HPP
#define MPM_V3_JOB_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <iostream>
#include <fstream>

#include "mpm_objects.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"

class Job{
public:
    //i'm not sure if these are ever used, but who cares
    static const int ZERO = 0;
    static const int ONES = 1;
    static const int IDENTITY = 2;

    //definitions of job types
    static const int JOB_1D = KinematicTensor::TENSOR_1D;
    static const int JOB_2D = KinematicTensor::TENSOR_2D;
    static const int JOB_3D = KinematicTensor::TENSOR_3D;
    static const int JOB_AXISYM = KinematicTensor::TENSOR_AXISYM;
    static const int JOB_2D_OOP = KinematicTensor::TENSOR_2D_OOP;

    //time and job type
    double t0, t, dt;
    int JOB_TYPE = JOB_3D;
    int DIM = 3;

    //list of bodies and conacts
    std::vector<int> activeBodies;
    std::vector<std::unique_ptr<Body>> bodies;
    std::vector<int> activeContacts;
    std::vector<std::unique_ptr<Contact>> contacts;

    //job objects
    std::unique_ptr<Serializer> serializer;
    std::unique_ptr<Driver> driver;
    std::unique_ptr<Solver> solver;
    std::unique_ptr<Grid> grid;

    void assignJobType(int input){
        JOB_TYPE = input;
        if (JOB_TYPE == JOB_1D){
            DIM = 1;
        } else if (JOB_TYPE == JOB_2D){
            DIM = 2;
        } else if (JOB_TYPE == JOB_3D){
            DIM = 3;
        } else if (JOB_TYPE == JOB_2D_OOP){
            DIM = 3;
        } else if (JOB_TYPE == JOB_AXISYM){
            DIM = 3;
        } else {
            std::cerr << "Job doesn't have defined type for input " << JOB_TYPE << "." << std::endl;
        }
    }

    //initialize job
    void init(){
        //called by config object after properties have been assigned
        std::cout << "Job properties (JOB_TYPE = " << JOB_TYPE << ", t = " << t << ", dt = " << dt << ")." << std::endl;
        std::cout << "Job Initialized." << std::endl;

        //call initialization on sub-objects in order of precedence
        serializer->init(this);
        driver->init(this);
        solver->init(this);

        grid->init(this);

        for (int i=0;i<bodies.size();i++){
            bodies[i]->init(this);
        }

        for (int i=0;i<contacts.size();i++){
            contacts[i]->init(this);
        }

        return;
    }
};

#endif //MPM_V3_JOB_HPP

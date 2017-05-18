//
// Created by aaron on 5/5/17.
// config.hpp
//

#ifndef MPM_3D_CONFIG_HPP
#define MPM_3D_CONFIG_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

class Job;

class Config{
public:
    //config properties here
    std::string filename; //name of file
    std::string filepath; //directory of file for access

    //objects here

    //config object specific functions
    Config();
    int configInit(std::string,std::string); //initialize config object from filepath and filename

    int configConfigureJob(Job*); //configure job object
};

#endif //MPM_3D_CONFIG_HPP

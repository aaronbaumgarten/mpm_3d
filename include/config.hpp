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
    std::string file; //name of file
    std::string mainpath; //director of main function

    //objects here

    //config object specific functions
    Config();
    int configInit(std::string); //initialize config object from filepath and filename

    void configSetMainPath(std::string); //set path to main
    void configCheckConfigFile(std::string); //check headers in config file

    int configConfigureJob(Job*); //configure job object
};

#endif //MPM_3D_CONFIG_HPP

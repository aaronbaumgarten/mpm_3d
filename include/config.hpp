//
// Created by aaron on 5/10/18.
// config.hpp
//

#ifndef MPM_V3_CONFIG_HPP
#define MPM_V3_CONFIG_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "job.hpp"

class Configurater{
public:
    std::string file; //name of config file
    std::string mainpath; //directory of main function

    void init(); //initialization
    void setMainPath(std::string); //set path to main
    void checkConfigFile(std::string); //check headers and report problems

    int configureJob(Job*); //configure job object;
};

#endif //MPM_V3_CONFIG_HPP

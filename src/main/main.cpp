//
// Created by aaron on 5/5/17.
// main.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <memory>
#include <Eigen/Core>

#include "job.hpp"
#include "config.hpp"

void usage(char* program_name){
    std::cout << program_name << " [OPTION]" << std::endl;
    std::cout << "    OPTION is any of the following (only one expected)." << std::endl;
    std::cout << "        -c CFGFILE, run simulation using specified configuration file." << std::endl;
    std::cout << "        -d, run default simulation using default files." << std::endl;
    std::cout << "        -r SIMFILE, **TO BE ADDED**." << std::endl;
    return;
}

int main(int argc, char *argv[]) {

    //set threads for eigen
    Eigen::setNbThreads(0);

    std::cout << "Hello, World!" << std::endl;
    std::cout << argv[0] << std::endl;

    Job *job(new Job);
    Config config();

    delete(job);
    return 0;
}
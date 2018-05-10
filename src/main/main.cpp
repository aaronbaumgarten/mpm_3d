//
// Created by aaron on 4/22/18.
// main.cpp
//

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "test.hpp"
#include "mpm_objects.hpp"

void usage(char* program_name){
    std::cout << program_name << " [OPTION]" << std::endl;
    std::cout << "    OPTION is any of the following (only one expected)." << std::endl;
    std::cout << "        -c CFGFILE, run simulation using specified configuration file." << std::endl;
    std::cout << "        -d, run default simulation using default files." << std::endl;
    std::cout << "        -r SIMFILE, **TO BE ADDED**." << std::endl;
    return;
}

int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;
    std::cout << argv[0] << std::endl;

    //read command
    std::cout << argc << " arguments given." << std::endl;

    algebra_test();

    std::cout << "Exiting." << std::endl;
    return 0;
}
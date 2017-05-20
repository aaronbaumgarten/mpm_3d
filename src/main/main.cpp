//
// Created by aaron on 5/5/17.
// main.cpp
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <signal.h>
#include <string>
#include <vector>
#include <memory>
#include <Eigen/Core>

#include "stringparser.hpp"

#include "job.hpp"
#include "config.hpp"

#include "serializer.hpp"
#include "driver.hpp"
#include "solver.hpp"

#include "grid.hpp"
#include "contact.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "points.hpp"

#include "material.hpp"

void usage(char* program_name){
    std::cout << program_name << " [OPTION]" << std::endl;
    std::cout << "    OPTION is any of the following (only one expected)." << std::endl;
    std::cout << "        -c CFGFILE, run simulation using specified configuration file." << std::endl;
    std::cout << "        -d, run default simulation using default files." << std::endl;
    std::cout << "        -r SIMFILE, **TO BE ADDED**." << std::endl;
    return;
}

std::vector<Serializer*> master_list_serializers;
std::vector<Job*> master_list_jobs;

void sigint_handler(int s){
    std::cout << "SIGINT received." << std::endl;
    for (size_t i=0; i<master_list_serializers.size(); i++){
        std::cout << "Saving job " << i << "." << std::endl;
        master_list_serializers[i]->serializerSaveState(master_list_jobs[i]);
        delete(master_list_jobs[i]);
    }
    std::cout << "Exiting." << std::endl;
    exit(0);
}

int main(int argc, char *argv[]) {

    //set threads for eigen
    Eigen::setNbThreads(0);

    std::cout << "Hello, World!" << std::endl;
    std::cout << argv[0] << std::endl;

    //create job and config objects
    Job *job(new Job);
    Config config;

    //set main program path in config
    config.configSetMainPath(std::string(argv[0]));

    //read command line args and initialize sim
    if (argc < 2){
        std::cout << "No arguments passed, expected 1. Exiting." << std::endl;
        usage(argv[0]);
        delete(job);
        exit(0);
    } else {
        std::vector<std::string> inputOptions = {"-c","-d","-r"};
        std::string inputArg(argv[1]);
        switch (StringParser::stringFindStringID(inputOptions,inputArg)){
            case -1:
                //unexpected arg
                std::cout << "Unexpected argument \"" << inputArg << "\". Exiting." << std::endl;
                usage(argv[0]);
                delete(job);
                exit(0);
                break;
            case 0:
                // -c CFGFILE
                config.configInit(argv[2]);
                config.configCheckConfigFile(argv[2]);
                if(!(config.configConfigureJob(job))){delete(job); exit(0);}
                break;
            case 1:
                // -d (default)
                job->dt = 1e-3;
                job->t = 0;

                std::cout << "\"" << inputArg << "\" not yet implemented. Exiting." << std::endl;
                usage(argv[0]);
                delete(job);
                exit(0);

                break;
            case 2:
                // -r SIMFILE (to be added)
                std::cout << "\"" << inputArg << "\" not yet implemented. Exiting." << std::endl;
                usage(argv[0]);
                delete(job);
                exit(0);
                break;

        }
    }

    //initialize job (it will call all internal initializers
    job->jobInit();

    size_t num_points;
    std::cout << "Job Initialized." << std::endl;
    std::cout << "  Bodies: " << job->bodies.size() << std::endl;
    for (size_t b=0;b<job->bodies.size();b++){
        num_points += job->bodies[b].points.x.size();
    }
    std::cout << "  Points: " << num_points << std::endl;
    std::cout << "  Nodes: " << job->grid.node_count << std::endl;
    std::cout << "  Contacts: " << job->contacts.size() << std::endl;

    //setup sigint handling
    master_list_serializers.push_back(&(job->serializer));
    master_list_jobs.push_back(job);
    signal(SIGINT, sigint_handler);

    //start simulation
    job->driver.driverRun(job);

    delete(job);
    return 0;
}
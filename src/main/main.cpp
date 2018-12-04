//
// Created by aaron on 4/22/18.
// main.cpp
//

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "test.hpp"
#include "mpm_objects.hpp"
#include "registry.hpp"
#include "job.hpp"
#include "config.hpp"
#include "signal_handler.hpp"

void usage(char* program_name){
    std::cout << program_name << " [OPTION]" << std::endl;
    std::cout << "    OPTION is any of the following (only one expected)." << std::endl;
    std::cout << "        -c CFGFILE, run simulation using specified configuration file." << std::endl;
    //std::cout << "        -d, run default simulation using default files." << std::endl;
    std::cout << "        -r SIMFILE, **TO BE ADDED**." << std::endl;
    return;
}

int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;
    std::cout << argv[0] << std::endl;

    //algebra_test();
    //map_test();

    //create job and configurator objects
    Job *job = new(Job);
    Configurator config;

    //set main program path in config
    config.setMainPath(std::string(argv[0]));

    //read command line args and initialize sim
    if (argc < 2){
        std::cout << "No arguments given. Exiting." << std::endl;
        usage(argv[0]);
        delete(job);
        exit(0);
    } else {
        std::vector<std::string> inputOptions = {"-c","-r"};
        std::string inputArg(argv[1]);
        switch (Parser::findStringID(inputOptions,inputArg)){
            case -1:
                //unexpected arg
                std::cout << "Unexpected argument \"" << inputArg << "\". Exiting." << std::endl;
                usage(argv[0]);
                delete(job);
                exit(0);
                break;
            case 0:
                // -c CFGFILE
                config.init(argv[2]);
                config.checkConfigFile(argv[2]);
                if(!(config.configureJob(job))){delete(job); exit(0);}

                //initialize job (it will call all internal initializers
                job->init();
                std::cout << "Job Initialized." << std::endl;

                break;
            case 1:
                // -r SIMFILE
                std::cout << "\'-r\' not implemented. Exiting." << std::endl;
                usage(argv[0]);
                delete(job);
                exit(0);
                break;
        }
    }

    //print job info
    int num_points = 0;
    std::cout << "\n";
    std::cout << "  Bodies: " << job->bodies.size() << std::endl;
    for (int b=0;b<job->bodies.size();b++){
        num_points += job->bodies[b]->points->x.size();
    }
    std::cout << "  Points: " << num_points << std::endl;
    std::cout << "  Nodes: " << job->grid->node_count << std::endl;
    std::cout << "  Elements: " << job->grid->element_count << std::endl;
    std::cout << "  Contacts: " << job->contacts.size() << std::endl;
    std::cout << "\n";

    //setup sigint handling
    signal(SIGINT, sigint_handler);

    //start simulation
    job->driver->run(job);

    //save simulation
    std::cout << std::endl << "Simulation Complete. Saving." << std::endl;
    job->serializer->saveState(job);

    std::cout << "Exiting." << std::endl;

    delete(job);
    return 0;
}
//
// Created by aaron on 7/11/17.
// signal_resolution.hpp
//

#ifndef MPM_V2_SIGNAL_RESOLUTION_HPP
#define MPM_V2_SIGNAL_RESOLUTION_HPP

#include <stdlib.h>
#include <signal.h>

static volatile bool INTERUPT_JOB = false;

void sigint_handler(int s){
    //ask to save and return to standard sigint handler
    signal(SIGINT, SIG_DFL);
    std::cout << std::endl << std::endl << "SIGINT received." << std::endl;
    std::cout << "Save? [y/n] ";
    char c = getchar();
    if (c == 'y' || c == 'Y'){
        std::cout << "Attempting to Exit Safely..." << std::endl;
        INTERUPT_JOB = true;
        //hope it shuts down on its own;
    } else {//if (c == 'n' || c == 'N') {
        std::cout << "Exiting." << std::endl;
        exit(0);
    }
}

bool check_interupt(){
    return INTERUPT_JOB;
}

#endif //MPM_V2_SIGNAL_RESOLUTION_HPP

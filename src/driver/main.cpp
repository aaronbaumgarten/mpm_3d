#include <iostream>
#include <stdlib.h>

#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "element.hpp"
#include "process.hpp"

//using namespace std;

int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;

    //initialize job
    job_t *job;

    //parse configuration files
    char *fileParticle = "s.particles";
    job->readParticles(fileParticle);

    //initialize and allocate memory

    //colorize for threading

    //pre-run setup

    //process_usl

    //serialize??

    //kill threads and cleanup

    return 0;
}
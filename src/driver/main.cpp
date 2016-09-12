#include <iostream>
#include <stdlib.h>
#include <memory>

#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "element.hpp"
#include "process.hpp"

//using namespace std;

int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;

    //initialize job
    job_t *job(new job_t);

    //parse configuration files
    char *fileParticle = "s.particles";
    char *fileNodes = "s.grid";
    job->importNodesandParticles(fileNodes,fileParticle);

    //initialize and allocate memory
    /* hard-coding material properties for initial run */
    double fp64_props[2] = {1e7,0.3};
    int *int_props = NULL;
    job->assignMaterials();
    for (size_t i=0;i<job->num_bodies;i++){
        job->bodies[i].defineMaterial(fp64_props,int_props);
    }

    //colorize for threading

    //pre-run setup

    //process_usl

    //serialize??

    //testing
    double sumX = 0;
    for (size_t i=0; i<job->num_nodes; i++){
        sumX += job->bodies[0].nodes[i].x[0];
    }
    std::cout << "test: " << sumX << " =? 515150\n";
    sumX = 0;
    for (size_t i=0; i<job->num_particles; i++){
        sumX += job->bodies[0].particles[i].id;
    }
    std::cout << "test: " << sumX << " =? " << 7999*8000/2 <<"\n";

    //kill threads and cleanup
    delete(job);

    return 0;
}
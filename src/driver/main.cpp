#include <iostream>
#include <stdlib.h>
#include <memory>

#include "test.hpp"

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

    job->createMappings();
    std::cout << "Mapping created (" << job->bodies[0].Sip.nonZeros() << ").\n";

    //colorize for threading

    //pre-run setup

    //process_usl
    std::cout << "node[515151].m: " << job->bodies[0].nodes[515151].m[0] << "\n";
    job->mapParticles2Grid();
    std::cout << "node[515151].m: " << job->bodies[0].nodes[515151].m[0] << "\n";

    //serialize??

    //testing

    //kill threads and cleanup
    delete(job);

    return 0;
}
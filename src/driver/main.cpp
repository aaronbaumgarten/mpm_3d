#include <iostream>
#include <stdlib.h>
#include <string>
#include <memory>
#include <Eigen/Core>

#include "test.hpp"

#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "element.hpp"
#include "process.hpp"
#include "mpmio.hpp"

#define T_STOP 1.0

int main(int argc, char *argv[]) {

    //set threads for eigen
    Eigen::setNbThreads(4);

    std::cout << "Hello, World!" << std::endl;

    //initialize job
    job_t *job(new job_t);

    //parse configuration files
    //char *fileParticle = "s.particles";
    //char *fileNodes = "s.grid";
    std::string fileParticle = "s.particles";
    std::string fileNodes = "s.grid";
    job->importNodesandParticles(fileNodes,fileParticle);

    //initialize and allocate memory
    /* hard-coding material properties for initial run */
    double fp64_props[2] = {1e7,0.3};
    int *int_props = NULL;
    job->assignMaterials();
    for (size_t i=0;i<job->num_bodies;i++){
        job->bodies[i].defineMaterial(fp64_props,int_props);
    }
    job->assignBoundaryConditions();

    job->createMappings();
    std::cout << "Mapping created (" << job->bodies[0].Sip.nonZeros() << ").\n";

    //colorize for threading

    //pre-run setup
    MPMio mpmOut;
    mpmOut.setDefaultFiles();
    mpmOut.setJob(job);
    mpmOut.writeFrameOutputHeader();

    //process_usl
    /*job->mpmStepUSLExplicit();
    std::cout << "Step completed (1).\n";*/
    while (job->t < T_STOP) {
        job->mpmStepUSLExplicit();
        mpmOut.writeFrame();
        std::cout << "Step completed (" << job->stepcount << ").\r" << std::flush;
    }
    std::cout << "\n";

    //serialize??

    //testing

    //kill threads and cleanup
    delete(job);

    return 0;
}
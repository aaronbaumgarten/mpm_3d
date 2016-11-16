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
    Eigen::setNbThreads(0);

    std::cout << "Hello, World!" << std::endl;

    //initialize job
    job_t *job(new job_t);
    job->dt = 1e-4;
    job->use_3d = 1;

    //parse configuration files
    //char *fileParticle = "s.particles";
    //char *fileNodes = "s.grid";
    std::string fileParticle = "s.particles";
    std::string fileNodes = "s.grid";
    job->importNodesandParticles(fileNodes,fileParticle);

    //initialize and allocate memory
    /* hard-coding material properties for initial run */
    double fp64_props[2] = {1e6,0.3};
    int *int_props = NULL;
    job->assignMaterials();
    for (size_t i=0;i<job->num_bodies;i++){
        job->bodies[i].defineMaterial(fp64_props,2,int_props,0);
    }
    job->assignBoundaryConditions();

    job->createMappings();
    std::cout << "Mapping created (" << job->bodies[0].Phi.nonZeros() << ").\n";

    //testMappingGradient(job);

    //colorize for threading

    //pre-run setup
    MPMio mpmOut;
    mpmOut.setDefaultFiles();
    mpmOut.setJob(job);
    mpmOut.setSampleRate(120.0);

    //process_usl
    while (job->t < T_STOP) {
        //job->mpmStepUSLExplicit();
        job->mpmStepUSLExplicitDebug();
        std::cout << "\rStep completed (" << job->stepcount << ")." << std::flush;
        if (job->t * mpmOut.sampleRate > mpmOut.sampledFrames) {
            mpmOut.writeFrame();
            mpmOut.sampledFrames += 1;
            std::cout << " Frame captured (" << mpmOut.sampledFrames << ")." << std::flush;
        }
    }
    std::cout << "\n";

    //serialize??

    //testing

    //kill threads and cleanup
    delete(job);

    return 0;
}
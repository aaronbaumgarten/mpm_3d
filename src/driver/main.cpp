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
    job->dt = 1e-3;
    //job->dt_base = job->dt;
    //job->dt_minimum = 1e-6;
    job->use_3d = 1;
    job->use_implicit = 1;
    job->use_cpdi = 1;
    job->newtonTOL = 1e-10;
    job->linearStepSize = 1e-5;

    //parse configuration files
    //char *fileParticle = "s.particles";
    //char *fileNodes = "s.grid";
    std::string fileParticle = "test.particles";
    std::string fileNodes = "test.grid";
    if(!(job->importNodesandParticles(fileNodes,fileParticle))){
        std::cout << "failed to import nodes and particles" << std::endl;
        exit(0);
    }

    //initialize and allocate memory
    /* hard-coding material properties for initial run */
    if(!(job->assignDefaultMaterials())){
        std::cout << "failed to assign materials" << std::endl;
        exit(0);
    }

    if(!(job->assignBoundaryConditions())){
        std::cout << "failed to assign boundary conditions" << std::endl;
        exit(0);
    }

    if(!(job->assignDefaultContacts())){
        std::cout << "failed to assign contact rules" << std::endl;
        exit(0);
    }

    if(!(job->createMappings())){
        std::cout << "failed to create mappings" << std::endl;
        exit(0);
    }
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
        if (job->use_3d==1) {
            if (job->use_implicit==1){
                job->mpmStepUSLImplicit();
            } else {
                job->mpmStepUSLExplicit();
            }
        } else {
            job->mpmStepUSLExplicit2D();
        }
        std::cout << "\rStep completed (" << job->stepcount << ")." << std::flush;
        if (job->t * mpmOut.sampleRate > mpmOut.sampledFrames) {
            mpmOut.writeFrame();
            mpmOut.sampledFrames += 1;
            std::cout << " Frame captured (" << mpmOut.sampledFrames-1 << ")." << std::flush;
        }

        /*std::ostringstream s;
        s << "pressure.csv";
        std::ofstream ffile(s.str(), std::ios::app);
        if (ffile.is_open()) {
            double P = 0;
            for (size_t i=0; i<job->bodies[0].p; i++){
                double pP;
                pP = (job->bodies[0].particles[i].T[XX]+job->bodies[0].particles[i].T[YY]+job->bodies[0].particles[i].T[ZZ])/3;
                P += pP*pP;
                //if p is nan set to zero
            }
            std::ostringstream line;
            line << std::sqrt(P) << "\n";
            ffile << line.str();
        }
        ffile.close();*/
    }
    std::cout << "\n";

    //serialize??

    //testing

    //kill threads and cleanup
    delete(job);

    return 0;
}
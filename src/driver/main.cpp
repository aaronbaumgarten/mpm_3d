#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <memory>
#include <Eigen/Core>

#include "test.hpp"

#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "element.hpp"
#include "process.hpp"
#include "mpmio.hpp"
#include "mpmconfig.hpp"


void usage(char* program_name){
    std::cout << program_name << " [OPTION]" << std::endl;
    std::cout << "    OPTION is any of the following (only one expected)." << std::endl;
    std::cout << "        -c CFGFILE, run simulation using specified configuration file." << std::endl;
    std::cout << "        -d, run default simulation using default files." << std::endl;
    std::cout << "        -r SIMFILE, **TO BE ADDED**." << std::endl;
    return;
}

int main(int argc, char *argv[]) {

    //set threads for eigen
    Eigen::setNbThreads(0);

    std::cout << "Hello, World!" << std::endl;
    std::cout << argv[0] << std::endl;

    //initialize job and objects
    job_t *job(new job_t);
    MPMio mpmOut;
    MPMconfig config;
    config.setMainPath(std::string(argv[0]));

    //read command line args and initialize simulation
    if (argc < 2) {
        std::cout << "No arguments passed, expected 1. Exiting." << std::endl;
        usage(argv[0]);
        delete(job);
        exit(0);
    } else {
        std::vector<std::string> inputOptions = {"-c","-d","-r"};
        std::string inputArg(argv[1]);
        switch (config.findStringID(inputOptions,inputArg)){
            case -1:
                //unexpected arg
                std::cout << "Unexpected argument \"" << inputArg << "\". Exiting." << std::endl;
                usage(argv[0]);
                delete(job);
                exit(0);
                break;
            case 0:
                // -c CFGFILE
                config.setConfigFile(argv[2]);
                config.checkConfigFile(argv[2]);
                if(!(config.configJob(job))){delete(job); exit(0);}
                if(!(config.configInput(job))){delete(job); exit(0);}
                if(!(config.configBoundary(job))){delete(job); exit(0);}
                if(!(config.configMaterial(job))){delete(job); exit(0);}
                if(!(config.configContact(job))){delete(job); exit(0);}
                if(!(config.configOutput(job,&(mpmOut)))){delete(job); exit(0);}
                break;
            case 1:
                // -d (default)
                job->dt = 1e-3;
                job->use_3d = 1;
                job->use_implicit = 1;
                job->use_cpdi = 1;
                job->use_smoothing = 0;
                job->newtonTOL = 1e-10;
                job->linearStepSize = 1e-5;

                //parse configuration files
                {
                    std::string fileParticle = "test.particles";
                    std::string fileNodes = "test.grid";
                    if (!(job->importNodesandParticles(fileNodes, fileParticle))) {
                        std::cout << "failed to import nodes and particles" << std::endl;
                        delete(job);
                        exit(0);
                    }
                }

                //initialize and allocate memory
                if (!(job->assignDefaultMaterials())) {
                    std::cout << "failed to assign materials" << std::endl;
                    delete(job);
                    exit(0);
                }
                if (!(job->assignDefaultBoundaryConditions())) {
                    std::cout << "failed to assign boundary conditions" << std::endl;
                    delete(job);
                    exit(0);
                }
                if (!(job->assignDefaultContacts())) {
                    std::cout << "failed to assign contact rules" << std::endl;
                    delete(job);
                    exit(0);
                }
                if (!(job->createMappings())) {
                    std::cout << "failed to create mappings" << std::endl;
                    delete(job);
                    exit(0);
                }
                std::cout << "Mapping created (" << job->bodies[0].Phi.nonZeros() << ").\n";

                //pre-run setup
                mpmOut.setDefaultFiles();
                mpmOut.setJob(job);
                mpmOut.setSampleRate(120.0);
                mpmOut.setSampleTime(1.0);
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


    //process_usl
    while (job->t < mpmOut.sampleTime) {
        //job->mpmStepUSLExplicit();
        if (job->use_3d == 1) {
            if (job->use_implicit == 1) {
                job->mpmStepUSLImplicit();
            } else {
                job->mpmStepUSLExplicit();
            }
        } else {
            if (job->use_implicit == 1){
                job->mpmStepUSLImplicit2D();
            } else {
                job->mpmStepUSLExplicit2D();
            }
        }

        std::cout << "\rStep completed (" << job->stepcount << ")." << std::flush;
        if (job->t * mpmOut.sampleRate > mpmOut.sampledFrames){// || job->stepcount > 1000) {
            mpmOut.writeFrame();
            mpmOut.sampledFrames += 1;
            std::cout << " Frame captured (" << mpmOut.sampledFrames - 1 << ")." << std::flush;
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

    //kill threads and cleanup
    delete (job);

    return 0;
}
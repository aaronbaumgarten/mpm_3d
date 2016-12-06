//
// Created by aaron on 10/29/16.
// mpmio.hpp
//

#ifndef MPM_3D_MPMIO_HPP
#define MPM_3D_MPMIO_HPP

#include <iostream>
#include <fstream>
#include <sstream>
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

class MPMio{
public:
    //files for reading and writing
    std::string inputFile; //read state
    std::string outputFile; //write state
    std::string frameDirectory; //output directory
    std::string frameFile; //write output frame for visualizers

    //sampling info
    double sampleRate;
    size_t sampledFrames;

    //job
    job_t* job;

    //x,y,z limits
    double xLimit;
    double yLimit;
    double zLimit;

    //functions
    MPMio(){}
    MPMio(std::string ifile,std::string ofile,std::string ffile, job_t* jobIn){
        inputFile = ifile;
        outputFile = ofile;
        frameFile = ffile;

        job = jobIn;

        xLimit = 2*jobIn->Lx;
        yLimit = 2*jobIn->Ly;
        zLimit = 2*jobIn->Lz;
    }

    void setDefaultFiles();
    void defineOutputFile(std::string);
    void defineInputFile(std::string);
    void defineFrameFile(std::string);
    void defineFrameDirectory(std::string);

    void setSampleRate(double);

    void setJob(job_t*);

    void readInput(); //read state saved in file to job object
    void readInput(job_t*);
    void writeOutput(); //write state to file
    void writeOutput(job_t*);

    void writeFrame(); //write frame to file
    void writeFrame(job_t*);
    void writeParticles(job_t*);
    void writeNodes(job_t *);
    void writeCorner(job_t*,size_t);

};


#endif //MPM_3D_MPMIO_HPP

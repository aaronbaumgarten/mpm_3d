//
// Created by aaron on 10/29/16.
// mpmio.cpp
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <memory>
#include <Eigen/Core>

#include "mpmio.hpp"

#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "element.hpp"
#include "process.hpp"

void MPMio::setDefaultFiles(){
    this->inputFile = "input.mpm";
    this->outputFile = "output.mpm";
    this->frameFile = "frames.mpm";
    return;
}

void MPMio::defineOutputFile(std::string ofile){
    this->outputFile = ofile;
    return;
}

void MPMio::defineInputFile(std::string ifile){
    this->inputFile = ifile;
    return;
}

void MPMio::defineFrameFile(std::string ffile){
    this->frameFile = ffile;
    return;
}

void MPMio::setJob(job_t *jobIn) {
    this->job = jobIn;
    return;
}

/*********************************************/
void MPMio::readInput(){
    this->readInput(this->job);
    return;
}

void MPMio::readInput(job_t* jobIn){
    //read input file to job object (set state of job to saved state)
    return;
}
/**********************************************/

/**********************************************/
void MPMio::writeOutput(){
    this->readInput(this->job);
    return;
}

void MPMio::writeOutput(job_t* jobIn){
    //write state of job to file
    return;
}
/***********************************************/

/***********************************************/
void MPMio::writeFrameOutputHeader(){
    this->writeFrameOutputHeader(this->job);
    return;
}

void MPMio::writeFrameOutputHeader(job_t* jobIn){
    //open frame file and overwrite
    std::ofstream ffile(this->frameFile, std::ios::trunc);

    //write header
    if (ffile.is_open()){
        ffile << "Job:\n";
        ffile << "Bodies: " << jobIn->num_bodies << "\n";
        ffile << "Particles: " << jobIn->num_particles << "\n";
        ffile << "Nodes: " << jobIn->num_nodes << "\n";
        ffile << "Elements: " << jobIn->num_elements << "\n";
        ffile << "t, particle[i].x, particle[i].y, particle[i].z\n";
        ffile.close();
    } else {
        std::cout << "Unable to open frame output file in MPMio object.\n";
    }

    return;
}
/************************************************/

/************************************************/
void MPMio::writeFrame(){
    this->writeFrame(this->job);
    return;
}

void MPMio::writeFrame(job_t* jobIn){
    //open frame file to write at end of file
    std::ofstream ffile(this->frameFile, std::ios::app);

    //write frame data
    if (ffile.is_open()){
        ffile << jobIn->t;
        for (size_t b=0; b<jobIn->num_bodies;b++){
            for (size_t i=0;i<jobIn->bodies[b].p;i++){
                ffile << ", " << jobIn->bodies[b].particles[i].x[0] << ", " << jobIn->bodies[b].particles[i].y[0] << ", " << jobIn->bodies[b].particles[i].z[0];
            }
        }
        ffile << "\n";
        ffile.close();
    } else {
        std::cout << "Unable to open frame output file in MPMio object.\n";
    }
}
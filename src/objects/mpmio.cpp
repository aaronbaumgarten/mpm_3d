//
// Created by aaron on 10/29/16.
// mpmio.cpp
//

#include <iostream>
#include <iomanip>
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
#include "tensor.hpp"

void MPMio::setDefaultFiles(){
    this->inputFile = "input.mpm";
    this->outputFile = "output.mpm";
    this->frameFile = "frames.mpm";
    this->frameDirectory = "output";
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

void MPMio::defineFrameDirectory(std::string fdir) {
    this->frameDirectory = fdir;
    return;
}

void MPMio::setSampleRate(double rateIn) {
    this->sampleRate = rateIn;
    this->sampledFrames = 0;
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

void MPMio::writeFrame(){
    this->writeFrame(this->job);
    return;
}

void MPMio::writeFrame(job_t* jobIn){
    //open frame file to write
    std::ostringstream s;
    s << this->frameDirectory << "/" << this->frameFile << "." << std::setw(10) << std::setfill('0') << this->sampledFrames << ".vtk";
    std::ofstream ffile(s.str(), std::ios::trunc);

    //write frame data
    if (ffile.is_open()){
        std::ostringstream fheader;
        fheader << "Frame: " << this->sampledFrames << ", Time: " << jobIn->t << "\n";
        std::ostringstream numPoints;
        numPoints << jobIn->num_particles;
        std::ostringstream strSize;
        strSize << (jobIn->num_particles*2);
        ffile << "# vtk DataFile Version 3.0\n";
        ffile << fheader.str();
        ffile << "ASCII\n";
        ffile << "DATASET UNSTRUCTURED_GRID\n";

        ffile << "POINTS " << numPoints.str() << " double\n";
        for (size_t b=0; b<jobIn->num_bodies;b++){
            for (size_t i=0;i<jobIn->bodies[b].p;i++){
                //position
                ffile << jobIn->bodies[b].particles[i].x[0] << " ";
                ffile << jobIn->bodies[b].particles[i].y[0] << " ";
                ffile << jobIn->bodies[b].particles[i].z[0] << "\n";
            }
        }

        ffile << "CELLS " << numPoints.str() << " " << strSize.str() << "\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << "1 " << i << "\n";
                ffile << line.str();
            }
        }

        ffile << "CELL_TYPES " << numPoints.str() << "\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                ffile << "1\n";
            }
        }

        ffile << "POINT_DATA " << numPoints.str() << "\n";
        ffile << "SCALARS mass double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << jobIn->bodies[b].particles[i].m[0] << "\n";
                ffile << line.str();
            }
        }

        ffile << "SCALARS volume double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << jobIn->bodies[b].particles[i].v[0] << "\n";
                ffile << line.str();
            }
        }

        ffile << "SCALARS density double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << jobIn->bodies[b].particles[i].m[0]/jobIn->bodies[b].particles[i].v[0] << "\n";
                ffile << line.str();
            }
        }

        ffile << "VECTORS velocity double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << jobIn->bodies[b].particles[i].x_t[0] << " " << jobIn->bodies[b].particles[i].y_t[0] << " " << jobIn->bodies[b].particles[i].z_t[0] << "\n";
                ffile << line.str();
            }
        }

        ffile << "VECTORS force double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << jobIn->bodies[b].particles[i].bx[0] << " " << jobIn->bodies[b].particles[i].by[0] << " " << jobIn->bodies[b].particles[i].bz[0] << "\n";
                ffile << line.str();
            }
        }

        ffile << "SCALARS pressure double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                double P = 0;
                tensor_trace3(&P,jobIn->bodies[b].particles[i].T);
                line << P << "\n";
                ffile << line.str();
            }
        }

        ffile << "TENSORS stress double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << jobIn->bodies[b].particles[i].T[XX] << " " << jobIn->bodies[b].particles[i].T[XY] << " " << jobIn->bodies[b].particles[i].T[XZ] << "\n";
                line << jobIn->bodies[b].particles[i].T[YX] << " " << jobIn->bodies[b].particles[i].T[YY] << " " << jobIn->bodies[b].particles[i].T[YZ] << "\n";
                line << jobIn->bodies[b].particles[i].T[ZX] << " " << jobIn->bodies[b].particles[i].T[ZY] << " " << jobIn->bodies[b].particles[i].T[ZZ] << "\n";
                ffile << line.str() << std::endl;
            }
        }

        ffile << "TENSORS F double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << jobIn->bodies[b].particles[i].F[XX] << " " << jobIn->bodies[b].particles[i].F[XY] << " " << jobIn->bodies[b].particles[i].F[XZ] << "\n";
                line << jobIn->bodies[b].particles[i].F[YX] << " " << jobIn->bodies[b].particles[i].F[YY] << " " << jobIn->bodies[b].particles[i].F[YZ] << "\n";
                line << jobIn->bodies[b].particles[i].F[ZX] << " " << jobIn->bodies[b].particles[i].F[ZY] << " " << jobIn->bodies[b].particles[i].F[ZZ] << "\n";
                ffile << line.str() << std::endl;
            }
        }

        ffile << "TENSORS Fp double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << jobIn->bodies[b].particles[i].Fp[XX] << " " << jobIn->bodies[b].particles[i].Fp[XY] << " " << jobIn->bodies[b].particles[i].Fp[XZ] << "\n";
                line << jobIn->bodies[b].particles[i].Fp[YX] << " " << jobIn->bodies[b].particles[i].Fp[YY] << " " << jobIn->bodies[b].particles[i].Fp[YZ] << "\n";
                line << jobIn->bodies[b].particles[i].Fp[ZX] << " " << jobIn->bodies[b].particles[i].Fp[ZY] << " " << jobIn->bodies[b].particles[i].Fp[ZZ] << "\n";
                ffile << line.str() << std::endl;
            }
        }

        ffile.close();
    } else {
        std::cout << "Unable to open frame output file in MPMio object.\n";
    }
}
/************************************************/
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

void MPMio::setSampleTime(double timeIn) {
    this->sampleTime = timeIn;
    return;
}

void MPMio::setJob(job_t *jobIn) {
    this->job = jobIn;

    this->xLimit = jobIn->Lx;
    this->yLimit = jobIn->Ly;
    this->zLimit = jobIn->Lz;

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
    /*for (size_t i=0;i<8;i++){
        this->writeCorner(this->job,i);
    }*/
    return;
}

void MPMio::writeFrame(job_t* jobIn) {
    writeParticles(jobIn);
    writeNodes(jobIn);
    return;
}

void MPMio::writeParticles(job_t* jobIn) {
    // set limits on position
    this->xLimit = jobIn->Lx;
    this->yLimit = jobIn->Ly;
    this->zLimit = jobIn->Lz;

    /*double xTest = -2.20021e46;
    std::cout << std::fabs(xTest) << "\n" << xTest << "\n";
    if (std::isnan(xTest) || std::isinf(xTest) || std::fabs(xTest)>this->xLimit){
        xTest = this->xLimit;
    }
    std::cout << xTest << "\n";*/

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
                double x = jobIn->bodies[b].particles.x[i];
                double y = jobIn->bodies[b].particles.y[i];
                double z = jobIn->bodies[b].particles.z[i];
                //ffile << jobIn->bodies[b].particles[i].x[0] << " ";
                //ffile << jobIn->bodies[b].particles[i].y[0] << " ";
                //ffile << jobIn->bodies[b].particles[i].z[0] << "\n";

                if (std::isnan(x) || std::isinf(x) || x<0){
                    x = 0;
                } else if (x>this->xLimit){
                    x = this->xLimit;
                }
                if (std::isnan(y) || std::isinf(y) || y<0) {
                    y = 0;
                } else if (y>this->yLimit) {
                    y = this->yLimit;
                }
                if (std::isnan(z) || std::isinf(z) || job->use_3d!=1 || z<0) {
                    z = 0; // if 2d set z to 0
                } else if (z>this->zLimit) {
                    z = this->zLimit;
                }
                ffile << x << " " << y << " " << z << "\n";
            }
        }

        double countCell = 0;
        ffile << "CELLS " << numPoints.str() << " " << strSize.str() << "\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << "1 " << countCell << "\n";
                ffile << line.str();
                countCell += 1;
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
                line << jobIn->bodies[b].particles.m[i] << "\n";
                ffile << line.str();
            }
        }

        ffile << "SCALARS body double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << b << "\n";
                ffile << line.str();
            }
        }

        ffile << "SCALARS volume double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                //set volume to 0 if nan
                if (std::isnan(jobIn->bodies[b].particles.v[i]) || std::isinf(jobIn->bodies[b].particles.v[i])){
                    line << 0 << "\n";
                } else {
                    line << jobIn->bodies[b].particles.v[i] << "\n";
                }
                ffile << line.str();
            }
        }

        ffile << "SCALARS density double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                //set density to zero if nan
                if (std::isnan(jobIn->bodies[b].particles.m[i]/jobIn->bodies[b].particles.v[i]) || std::isinf(jobIn->bodies[b].particles.m[i]/jobIn->bodies[b].particles.v[i])){
                    line << 0 << "\n";
                } else {
                    line << jobIn->bodies[b].particles.m[i] / jobIn->bodies[b].particles.v[i] << "\n";
                }
                ffile << line.str();
            }
        }

        ffile << "VECTORS velocity double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                //set velocity to zeros if nan
                double x_t = jobIn->bodies[b].particles.x_t[i];
                double y_t = jobIn->bodies[b].particles.y_t[i];
                double z_t = jobIn->bodies[b].particles.z_t[i];
                if (std::isnan(x_t) || std::isinf(x_t)){
                    x_t = 0;
                }
                if (std::isnan(y_t) || std::isinf(y_t)){
                    y_t = 0;
                }
                if (std::isnan(z_t) || std::isinf(z_t)){
                    z_t = 0;
                }
                line << x_t << " " << y_t << " " << z_t << "\n";
                ffile << line.str();
            }
        }

        ffile << "VECTORS bodyForce double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                //set body force to zero is nan
                double bx = jobIn->bodies[b].particles.bx[i];
                double by = jobIn->bodies[b].particles.by[i];
                double bz = jobIn->bodies[b].particles.bz[i];
                if (std::isnan(bx) || std::isinf(bx)){
                    bx = 0;
                }
                if (std::isnan(by) || std::isinf(by)){
                    by = 0;
                }
                if (std::isnan(bz) || std::isinf(bz)){
                    bz = 0;
                }
                line << bx << " " << by << " " << bz << "\n";
                ffile << line.str();
            }
        }

        ffile << "SCALARS pressure double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                double P = (jobIn->bodies[b].particles.T(i,XX)+jobIn->bodies[b].particles.T(i,YY)+jobIn->bodies[b].particles.T(i,ZZ));
                if (jobIn->use_3d==1){
                    P = -P/3.0;
                } else {
                    P -= jobIn->bodies[b].particles.T(i,ZZ);
                    P = -P/2.0;
                }
                //if p is nan set to zero
                if (std::isnan(P) || std::isinf(P)){
                    P=0;
                }
                line << P << "\n";
                ffile << line.str();
            }
        }

        ffile << "TENSORS stress double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                //set nan stresses to 0
                double txx = jobIn->bodies[b].particles.T(i,XX);
                double txy = jobIn->bodies[b].particles.T(i,XY);
                double txz = jobIn->bodies[b].particles.T(i,XZ);
                double tyx = jobIn->bodies[b].particles.T(i,YX);
                double tyy = jobIn->bodies[b].particles.T(i,YY);
                double tyz = jobIn->bodies[b].particles.T(i,YZ);
                double tzx = jobIn->bodies[b].particles.T(i,ZX);
                double tzy = jobIn->bodies[b].particles.T(i,ZY);
                double tzz = jobIn->bodies[b].particles.T(i,ZZ);

                if (std::isnan(txx) || std::isinf(txx)){
                    txx = 0;
                }
                if (std::isnan(txy) || std::isinf(txy)){
                    txy = 0;
                }
                if (std::isnan(txz) || std::isinf(txz)){
                    txz = 0;
                }
                if (std::isnan(tyx) || std::isinf(tyx)){
                    tyx = 0;
                }
                if (std::isnan(tyy) || std::isinf(tyy)){
                    tyy = 0;
                }
                if (std::isnan(tyz) || std::isinf(tyz)){
                    tyz = 0;
                }
                if (std::isnan(tzx) || std::isinf(tzx)){
                    tzx = 0;
                }
                if (std::isnan(tzy) || std::isinf(tzy)){
                    tzy = 0;
                }
                if (std::isnan(tzz) || std::isinf(tzz)){
                    tzz = 0;
                }

                //line << jobIn->bodies[b].particles[i].T[XX] << " " << jobIn->bodies[b].particles[i].T[XY] << " " << jobIn->bodies[b].particles[i].T[XZ] << "\n";
                //line << jobIn->bodies[b].particles[i].T[YX] << " " << jobIn->bodies[b].particles[i].T[YY] << " " << jobIn->bodies[b].particles[i].T[YZ] << "\n";
                //line << jobIn->bodies[b].particles[i].T[ZX] << " " << jobIn->bodies[b].particles[i].T[ZY] << " " << jobIn->bodies[b].particles[i].T[ZZ] << "\n";

                line << txx << " " << txy << " " << txz << "\n";
                line << tyx << " " << tyy << " " << tyz << "\n";
                line << tzx << " " << tzy << " " << tzz << "\n";

                ffile << line.str() << std::endl;
            }
        }

        ffile << "TENSORS gradVelocity double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;

                //set nan gradients to 0
                double lxx = jobIn->bodies[b].particles.L(i,XX);
                double lxy = jobIn->bodies[b].particles.L(i,XY);
                double lxz = jobIn->bodies[b].particles.L(i,XZ);
                double lyx = jobIn->bodies[b].particles.L(i,YX);
                double lyy = jobIn->bodies[b].particles.L(i,YY);
                double lyz = jobIn->bodies[b].particles.L(i,YZ);
                double lzx = jobIn->bodies[b].particles.L(i,ZX);
                double lzy = jobIn->bodies[b].particles.L(i,ZY);
                double lzz = jobIn->bodies[b].particles.L(i,ZZ);

                if (std::isnan(lxx) || std::isinf(lxx)){
                    lxx = 0;
                }
                if (std::isnan(lxy) || std::isinf(lxy)){
                    lxy = 0;
                }
                if (std::isnan(lxz) || std::isinf(lxz)){
                    lxz = 0;
                }
                if (std::isnan(lyx) || std::isinf(lyx)){
                    lyx = 0;
                }
                if (std::isnan(lyy) || std::isinf(lyy)){
                    lyy = 0;
                }
                if (std::isnan(lyz) || std::isinf(lyz)){
                    lyz = 0;
                }
                if (std::isnan(lzx) || std::isinf(lzx)){
                    lzx = 0;
                }
                if (std::isnan(lzy) || std::isinf(lzy)){
                    lzy = 0;
                }
                if (std::isnan(lzz) || std::isinf(lzz)){
                    lzz = 0;
                }

                //line << jobIn->bodies[b].particles[i].L[XX] << " " << jobIn->bodies[b].particles[i].L[XY] << " " << jobIn->bodies[b].particles[i].L[XZ] << "\n";
                //line << jobIn->bodies[b].particles[i].L[YX] << " " << jobIn->bodies[b].particles[i].L[YY] << " " << jobIn->bodies[b].particles[i].L[YZ] << "\n";
                //line << jobIn->bodies[b].particles[i].L[ZX] << " " << jobIn->bodies[b].particles[i].L[ZY] << " " << jobIn->bodies[b].particles[i].L[ZZ] << "\n";

                line << lxx << " " << lxy << " " << lxz << "\n";
                line << lyx << " " << lyy << " " << lyz << "\n";
                line << lzx << " " << lzy << " " << lzz << "\n";
                ffile << line.str() << std::endl;
            }
        }

        ffile << "TENSORS F double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;

                //set nan gradients to 0
                double fxx = jobIn->bodies[b].particles.F(i,XX);
                double fxy = jobIn->bodies[b].particles.F(i,XY);
                double fxz = jobIn->bodies[b].particles.F(i,XZ);
                double fyx = jobIn->bodies[b].particles.F(i,YX);
                double fyy = jobIn->bodies[b].particles.F(i,YY);
                double fyz = jobIn->bodies[b].particles.F(i,YZ);
                double fzx = jobIn->bodies[b].particles.F(i,ZX);
                double fzy = jobIn->bodies[b].particles.F(i,ZY);
                double fzz = jobIn->bodies[b].particles.F(i,ZZ);

                if (std::isnan(fxx) || std::isinf(fxx)){
                    fxx = 0;
                }
                if (std::isnan(fxy) || std::isinf(fxy)){
                    fxy = 0;
                }
                if (std::isnan(fxz) || std::isinf(fxz)){
                    fxz = 0;
                }
                if (std::isnan(fyx) || std::isinf(fyx)){
                    fyx = 0;
                }
                if (std::isnan(fyy) || std::isinf(fyy)){
                    fyy = 0;
                }
                if (std::isnan(fyz) || std::isinf(fyz)){
                    fyz = 0;
                }
                if (std::isnan(fzx) || std::isinf(fzx)){
                    fzx = 0;
                }
                if (std::isnan(fzy) || std::isinf(fzy)){
                    fzy = 0;
                }
                if (std::isnan(fzz) || std::isinf(fzz)){
                    fzz = 0;
                }

                //line << jobIn->bodies[b].particles[i].F[XX] << " " << jobIn->bodies[b].particles[i].F[XY] << " " << jobIn->bodies[b].particles[i].F[XZ] << "\n";
                //line << jobIn->bodies[b].particles[i].F[YX] << " " << jobIn->bodies[b].particles[i].F[YY] << " " << jobIn->bodies[b].particles[i].F[YZ] << "\n";
                //line << jobIn->bodies[b].particles[i].F[ZX] << " " << jobIn->bodies[b].particles[i].F[ZY] << " " << jobIn->bodies[b].particles[i].F[ZZ] << "\n";

                line << fxx << " " << fxy << " " << fxz << "\n";
                line << fyx << " " << fyy << " " << fyz << "\n";
                line << fzx << " " << fzy << " " << fzz << "\n";
                ffile << line.str() << std::endl;
            }
        }

        ffile << "TENSORS Fp double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;

                //set nan gradients to 0
                double fxx = jobIn->bodies[b].particles.Fp(i,XX);
                double fxy = jobIn->bodies[b].particles.Fp(i,XY);
                double fxz = jobIn->bodies[b].particles.Fp(i,XZ);
                double fyx = jobIn->bodies[b].particles.Fp(i,YX);
                double fyy = jobIn->bodies[b].particles.Fp(i,YY);
                double fyz = jobIn->bodies[b].particles.Fp(i,YZ);
                double fzx = jobIn->bodies[b].particles.Fp(i,ZX);
                double fzy = jobIn->bodies[b].particles.Fp(i,ZY);
                double fzz = jobIn->bodies[b].particles.Fp(i,ZZ);

                if (std::isnan(fxx) || std::isinf(fxx)){
                    fxx = 0;
                }
                if (std::isnan(fxy) || std::isinf(fxy)){
                    fxy = 0;
                }
                if (std::isnan(fxz) || std::isinf(fxz)){
                    fxz = 0;
                }
                if (std::isnan(fyx) || std::isinf(fyx)){
                    fyx = 0;
                }
                if (std::isnan(fyy) || std::isinf(fyy)){
                    fyy = 0;
                }
                if (std::isnan(fyz) || std::isinf(fyz)){
                    fyz = 0;
                }
                if (std::isnan(fzx) || std::isinf(fzx)){
                    fzx = 0;
                }
                if (std::isnan(fzy) || std::isinf(fzy)){
                    fzy = 0;
                }
                if (std::isnan(fzz) || std::isinf(fzz)){
                    fzz = 0;
                }

                //line << jobIn->bodies[b].particles[i].Fp[XX] << " " << jobIn->bodies[b].particles[i].Fp[XY] << " " << jobIn->bodies[b].particles[i].Fp[XZ] << "\n";
                //line << jobIn->bodies[b].particles[i].Fp[YX] << " " << jobIn->bodies[b].particles[i].Fp[YY] << " " << jobIn->bodies[b].particles[i].Fp[YZ] << "\n";
                //line << jobIn->bodies[b].particles[i].Fp[ZX] << " " << jobIn->bodies[b].particles[i].Fp[ZY] << " " << jobIn->bodies[b].particles[i].Fp[ZZ] << "\n";

                line << fxx << " " << fxy << " " << fxz << "\n";
                line << fyx << " " << fyy << " " << fyz << "\n";
                line << fzx << " " << fzy << " " << fzz << "\n";
                ffile << line.str() << std::endl;
            }
        }

        ffile << "TENSORS Be double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;

                //set nan gradients to 0
                double fxx = jobIn->bodies[b].particles.Be(i,XX);
                double fxy = jobIn->bodies[b].particles.Be(i,XY);
                double fxz = jobIn->bodies[b].particles.Be(i,XZ);
                double fyx = jobIn->bodies[b].particles.Be(i,YX);
                double fyy = jobIn->bodies[b].particles.Be(i,YY);
                double fyz = jobIn->bodies[b].particles.Be(i,YZ);
                double fzx = jobIn->bodies[b].particles.Be(i,ZX);
                double fzy = jobIn->bodies[b].particles.Be(i,ZY);
                double fzz = jobIn->bodies[b].particles.Be(i,ZZ);

                if (std::isnan(fxx) || std::isinf(fxx)){
                    fxx = 0;
                }
                if (std::isnan(fxy) || std::isinf(fxy)){
                    fxy = 0;
                }
                if (std::isnan(fxz) || std::isinf(fxz)){
                    fxz = 0;
                }
                if (std::isnan(fyx) || std::isinf(fyx)){
                    fyx = 0;
                }
                if (std::isnan(fyy) || std::isinf(fyy)){
                    fyy = 0;
                }
                if (std::isnan(fyz) || std::isinf(fyz)){
                    fyz = 0;
                }
                if (std::isnan(fzx) || std::isinf(fzx)){
                    fzx = 0;
                }
                if (std::isnan(fzy) || std::isinf(fzy)){
                    fzy = 0;
                }
                if (std::isnan(fzz) || std::isinf(fzz)){
                    fzz = 0;
                }

                //line << jobIn->bodies[b].particles[i].Fp[XX] << " " << jobIn->bodies[b].particles[i].Fp[XY] << " " << jobIn->bodies[b].particles[i].Fp[XZ] << "\n";
                //line << jobIn->bodies[b].particles[i].Fp[YX] << " " << jobIn->bodies[b].particles[i].Fp[YY] << " " << jobIn->bodies[b].particles[i].Fp[YZ] << "\n";
                //line << jobIn->bodies[b].particles[i].Fp[ZX] << " " << jobIn->bodies[b].particles[i].Fp[ZY] << " " << jobIn->bodies[b].particles[i].Fp[ZZ] << "\n";

                line << fxx << " " << fxy << " " << fxz << "\n";
                line << fyx << " " << fyy << " " << fyz << "\n";
                line << fzx << " " << fzy << " " << fzz << "\n";
                ffile << line.str() << std::endl;
            }
        }

        ffile << "TENSORS state double\n";
        //ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;

                //set nan gradients to 0
                double fxx = jobIn->bodies[b].particles.state(i,XX);
                double fxy = jobIn->bodies[b].particles.state(i,XY);
                double fxz = jobIn->bodies[b].particles.state(i,XZ);
                double fyx = jobIn->bodies[b].particles.state(i,YX);
                double fyy = jobIn->bodies[b].particles.state(i,YY);
                double fyz = jobIn->bodies[b].particles.state(i,YZ);
                double fzx = jobIn->bodies[b].particles.state(i,ZX);
                double fzy = jobIn->bodies[b].particles.state(i,ZY);
                double fzz = jobIn->bodies[b].particles.state(i,ZZ);

                if (std::isnan(fxx) || std::isinf(fxx)){
                    fxx = 0;
                }
                if (std::isnan(fxy) || std::isinf(fxy)){
                    fxy = 0;
                }
                if (std::isnan(fxz) || std::isinf(fxz)){
                    fxz = 0;
                }
                if (std::isnan(fyx) || std::isinf(fyx)){
                    fyx = 0;
                }
                if (std::isnan(fyy) || std::isinf(fyy)){
                    fyy = 0;
                }
                if (std::isnan(fyz) || std::isinf(fyz)){
                    fyz = 0;
                }
                if (std::isnan(fzx) || std::isinf(fzx)){
                    fzx = 0;
                }
                if (std::isnan(fzy) || std::isinf(fzy)){
                    fzy = 0;
                }
                if (std::isnan(fzz) || std::isinf(fzz)){
                    fzz = 0;
                }

                //line << jobIn->bodies[b].particles[i].Fp[XX] << " " << jobIn->bodies[b].particles[i].Fp[XY] << " " << jobIn->bodies[b].particles[i].Fp[XZ] << "\n";
                //line << jobIn->bodies[b].particles[i].Fp[YX] << " " << jobIn->bodies[b].particles[i].Fp[YY] << " " << jobIn->bodies[b].particles[i].Fp[YZ] << "\n";
                //line << jobIn->bodies[b].particles[i].Fp[ZX] << " " << jobIn->bodies[b].particles[i].Fp[ZY] << " " << jobIn->bodies[b].particles[i].Fp[ZZ] << "\n";

                line << fxx << " " << fxy << " " << fxz << "\n";
                line << fyx << " " << fyy << " " << fyz << "\n";
                line << fzx << " " << fzy << " " << fzz << "\n";
                ffile << line.str() << std::endl;
            }
        }

        ffile << "SCALARS state9 double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                double P = jobIn->bodies[b].particles.state(i,9);
                //if p is nan set to zero
                if (std::isnan(P) || std::isinf(P)){
                    P=0;
                }
                line << P << "\n";
                ffile << line.str();
            }
        }

        ffile << "SCALARS state10 double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                double P = jobIn->bodies[b].particles.state(i,10);
                //if p is nan set to zero
                if (std::isnan(P) || std::isinf(P)){
                    P=0;
                }
                line << P << "\n";
                ffile << line.str();
            }
        }

        ffile.close();
    } else {
        std::cout << "Unable to open frame output file in MPMio object.\n";
    }
}

void MPMio::writeNodes(job_t* jobIn) {
    for (size_t b=0; b<jobIn->num_bodies;b++) {
        //open frame file to write
        std::ostringstream s;
        s << this->frameDirectory << "/n" << b << this->frameFile << "." << std::setw(10) << std::setfill('0') <<
        this->sampledFrames << ".vtk";
        std::ofstream ffile(s.str(), std::ios::trunc);

        //write frame data
        if (ffile.is_open()) {
            std::ostringstream fheader;
            fheader << "Frame: " << this->sampledFrames << ", Time: " << jobIn->t << "\n";
            std::ostringstream numPoints;
            numPoints << jobIn->num_nodes;
            std::ostringstream strSize;
            strSize << (jobIn->num_nodes * 2);
            ffile << "# vtk DataFile Version 3.0\n";
            ffile << fheader.str();
            ffile << "ASCII\n";
            ffile << "DATASET UNSTRUCTURED_GRID\n";

            ffile << "POINTS " << numPoints.str() << " double\n";
            for (size_t i = 0; i < jobIn->bodies[b].n; i++) {
                //position
                double x = jobIn->bodies[b].nodes.x[i];
                double y = jobIn->bodies[b].nodes.y[i];
                double z = jobIn->bodies[b].nodes.z[i];
                //ffile << jobIn->bodies[b].particles[i].x[0] << " ";
                //ffile << jobIn->bodies[b].particles[i].y[0] << " ";
                //ffile << jobIn->bodies[b].particles[i].z[0] << "\n";

                if (std::isnan(x) || std::isinf(x)) {
                    x = 0;
                }
                if (std::isnan(y) || std::isinf(y)) {
                    y = 0;
                }
                if (std::isnan(z) || std::isinf(z) || job->use_3d != 1) {
                    z = 0; // if 2d set z to 0
                }
                ffile << x << " " << y << " " << z << "\n";
            }

            ffile << "CELLS " << numPoints.str() << " " << strSize.str() << "\n";
            for (size_t i = 0; i < jobIn->bodies[b].n; i++) {
                std::ostringstream line;
                line << "1 " << i << "\n";
                ffile << line.str();
            }

            ffile << "CELL_TYPES " << numPoints.str() << "\n";
            for (size_t i = 0; i < jobIn->bodies[b].n; i++) {
                ffile << "1\n";
            }

            ffile << "POINT_DATA " << numPoints.str() << "\n";
            ffile << "SCALARS mass double 1\n";
            ffile << "LOOKUP_TABLE default\n";
            for (size_t i = 0; i < jobIn->bodies[b].n; i++) {
                std::ostringstream line;
                line << jobIn->bodies[b].nodes.m[i] << "\n";
                ffile << line.str();
            }

            ffile << "VECTORS velocity double\n";
            //ffile << "LOOKUP_TABLE default\n";
            for (size_t i = 0; i < jobIn->bodies[b].n; i++) {
                std::ostringstream line;
                //set velocity to zeros if nan
                double x_t = jobIn->bodies[b].nodes.contact_x_t[i];
                double y_t = jobIn->bodies[b].nodes.contact_y_t[i];
                double z_t = jobIn->bodies[b].nodes.contact_z_t[i];
                if (std::isnan(x_t) || std::isinf(x_t)) {
                    x_t = 0;
                }
                if (std::isnan(y_t) || std::isinf(y_t)) {
                    y_t = 0;
                }
                if (std::isnan(z_t) || std::isinf(z_t)) {
                    z_t = 0;
                }
                line << x_t << " " << y_t << " " << z_t << "\n";
                ffile << line.str();
            }

            ffile << "VECTORS force double\n";
            //ffile << "LOOKUP_TABLE default\n";
            for (size_t i = 0; i < jobIn->bodies[b].n; i++) {
                std::ostringstream line;
                //set body force to zero is nan
                double bx = jobIn->bodies[b].nodes.contact_fx[i];
                double by = jobIn->bodies[b].nodes.contact_fy[i];
                double bz = jobIn->bodies[b].nodes.contact_fz[i];
                if (std::isnan(bx) || std::isinf(bx)) {
                    bx = 0;
                }
                if (std::isnan(by) || std::isinf(by)) {
                    by = 0;
                }
                if (std::isnan(bz) || std::isinf(bz)) {
                    bz = 0;
                }
                line << bx << " " << by << " " << bz << "\n";
                ffile << line.str();
            }

            ffile.close();
        } else {
            std::cout << "Unable to open frame output file in MPMio object.\n";
        }
    }
}

void MPMio::writeCorner(job_t* jobIn, size_t cornerID){
    //open frame file to write
    std::ostringstream s;
    s << this->frameDirectory << "/" << this->frameFile << cornerID << "." << std::setw(10) << std::setfill('0') << this->sampledFrames << ".vtk";
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
                ffile << jobIn->bodies[b].particles.corner_x(i,cornerID) << " ";
                ffile << jobIn->bodies[b].particles.corner_y(i,cornerID) << " ";
                ffile << jobIn->bodies[b].particles.corner_z(i,cornerID) << "\n";
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
                line << jobIn->bodies[b].particles.m[i] << "\n";
                ffile << line.str();
            }
        }

        ffile << "SCALARS body double 1\n";
        ffile << "LOOKUP_TABLE default\n";
        for (size_t b=0; b<jobIn->num_bodies; b++){
            for (size_t i=0; i<jobIn->bodies[b].p; i++){
                std::ostringstream line;
                line << b << "\n";
                ffile << line.str();
            }
        }

        ffile.close();
    } else {
        std::cout << "Unable to open frame output file in MPMio object.\n";
    }
}
/************************************************/
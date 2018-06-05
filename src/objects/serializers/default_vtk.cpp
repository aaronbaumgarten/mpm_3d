//
// Created by aaron on 5/10/18.
// default_vtk.cpp
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <sys/stat.h>
#include <math.h>
#include "mpm_objects.hpp"
#include "job.hpp"
#include "serializers.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"

/*----------------------------------------------------------------------------*/

void DefaultVTK::init(Job* job){
    if (str_props.size() < 3 || fp64_props.size() < 1){
        std::cout << job->serializer->str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 4 properties defined ({sampleRate}, {frameDirectory, outputDirectory, outputName}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        frameDirectory = Parser::makeDirectory(str_props[0]);
        outputDirectory = Parser::makeDirectory(str_props[1]);
        outputName = str_props[2];

        sampleRate = fp64_props[0];
        sampledFrames = 0;
        t_last_frame = 0;

        printf("Serializer properties (frameDirectory = %s, outputDirectory = %s, outputName = %s, sampleRate = %g).\n",
               frameDirectory.c_str(), outputDirectory.c_str(), outputName.c_str(), sampleRate);
    }

    std::cout << "Serializer Initialized." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
void DefaultVTK::writeDefaultPointHeader(Job *job, Body *body, std::ofstream &pfile, int SPEC) {
    int plen = body->points->x.size();

    pfile << "ASCII\n";
    pfile << "DATASET UNSTRUCTURED_GRID\n";

    pfile << "POINTS " << plen << " double\n";
    for (int i=0;i<plen;i++){
        //vtk files require x,y,z
        for (int pos = 0; pos < 3; pos++){
            if (pos < body->points->x.DIM && (body->points->active(i) != 0) && std::isfinite(body->points->x(i,pos))){
                pfile << body->points->x(i,pos) << " ";
            } else {
                pfile << "0 ";
            }
        }
        pfile << "\n";
    }

    pfile << "CELLS " << plen << " " << 2*plen << "\n";
    for (int i=0;i<plen;i++){
        pfile << "1 " << i << "\n";
    }

    pfile << "CELL_TYPES " << plen << "\n";
    for (int i=0;i<plen;i++){
        pfile << "1\n";
    }

    pfile << "POINT_DATA " << plen << "\n";

    return;
}

void DefaultVTK::writeDefaultNodeHeader(Job *job, Body *body, std::ofstream &nfile, int SPEC) {
    int nlen = body->nodes->x.size();

    nfile << "ASCII\n";
    nfile << "DATASET UNSTRUCTURED_GRID\n";

    nfile << "POINTS " << nlen << " double\n";
    for (int i=0;i<nlen;i++){
        //vtk files require x,y,z
        for (int pos = 0; pos < 3; pos++){
            if (pos < body->nodes->x.DIM && (body->nodes->active(i) != 0) && std::isfinite(body->nodes->x(i,pos))){
                nfile << body->nodes->x(i,pos) << " ";
            } else {
                nfile << "0 ";
            }
        }
        nfile << "\n";
    }

    nfile << "CELLS " << nlen << " " << 2*nlen << "\n";
    for (int i=0;i<nlen;i++){
        nfile << "1 " << i << "\n";
    }

    nfile << "CELL_TYPES " << nlen << "\n";
    for (int i=0;i<nlen;i++){
        nfile << "1\n";
    }

    nfile << "POINT_DATA " << nlen << "\n";
}

int DefaultVTK::writeFrame(Job* job){
    if ((job->t - t_last_frame) >= (1.0/sampleRate) || sampledFrames == 0){ //job->t >= sampledFrames/sampleRate){
        t_last_frame = sampledFrames/sampleRate;
        sampledFrames += 1;
        //write frame for each body
        for (int b=0;b<job->bodies.size();b++) {
            currentBody = job->bodies[b].get();

            //open point file
            std::stringstream ss;
            ss << "fpd." << job->bodies[b]->id << "." << job->bodies[b]->name << "." << std::setw(10) << std::setfill('0') << (sampledFrames-1) << ".vtk";
            pfilename = ss.str();
            pfile = std::ofstream(frameDirectory+pfilename,std::ios::trunc);

            ss.str(std::string());
            ss.clear();

            //open node file
            ss << "fnd." << job->bodies[b]->id << "." << job->bodies[b]->name << "." << std::setw(10) << std::setfill('0') << (sampledFrames-1) << ".vtk";
            nfilename = ss.str();
            nfile = std::ofstream(frameDirectory+nfilename,std::ios::trunc);

            //set length of frame data
            plen = job->bodies[b]->points->x.size();
            nlen = job->bodies[b]->nodes->x.size();

            if (pfile.is_open()){
                pfile << "# vtk DataFile Version 3.0\n";
                pfile << "Frame: " << (sampledFrames-1) << ", Time: " << job->t << "\n";

                job->bodies[b]->points->writeHeader(job, currentBody, job->serializer.get(), pfile, VTK);

                //scalars, vectors and tensors here
            } else {
                std::cerr << "Could not open point frame: " << frameDirectory+pfilename << " !" << std::endl;
            }

            if (nfile.is_open()){
                nfile << "# vtk DataFile Version 3.0\n";
                nfile << "Frame: " << sampledFrames << ", Time: " << job->t << "\n";

                job->grid->writeHeader(job, currentBody, job->serializer.get(), nfile, VTK);

                //scalars, vectors and tensors here
            } else {
                std::cerr << "Could not open node frame: " << frameDirectory+nfilename << " !" << std::endl;
            }

            //call objects to write frame data
            if (job->activeBodies[b] != 0) {
                job->bodies[b]->points->writeFrame(job, job->bodies[b].get(), job->serializer.get());
                job->bodies[b]->nodes->writeFrame(job, job->bodies[b].get(), job->serializer.get());
                job->bodies[b]->material->writeFrame(job, job->bodies[b].get(), job->serializer.get());
                job->bodies[b]->boundary->writeFrame(job, job->bodies[b].get(), job->serializer.get());
                job->grid->writeFrame(job, job->serializer.get());
            }
            for (int c=0;c<job->contacts.size();c++) {
                job->contacts[c]->writeFrame(job, job->serializer.get());
            }

            pfile.close();
            nfile.close();
        }

        return 1;
    } else {
        //do not write frame
        return 0;
    }
}

/*----------------------------------------------------------------------------*/

void DefaultVTK::writeScalarArray(Eigen::VectorXd& scalarArray, std::string name){
    //check length of array vs. length of open files and write to correct file
    if (pfile.is_open() && scalarArray.rows() == plen){
        //write to point file
        pfile << "SCALARS " << name << " double 1\n";
        pfile << "LOOKUP_TABLE default\n";
        for (int i = 0; i < plen; i++){
            if (currentBody->points->active(i) == 1 && std::isfinite(scalarArray(i))) {
                pfile << scalarArray(i) << "\n";
            } else {
                pfile << "0" << "\n";
            }
        }
    }

    if (nfile.is_open() && scalarArray.rows() == nlen){
        //write to point file
        nfile << "SCALARS " << name << " double 1\n";
        nfile << "LOOKUP_TABLE default\n";
        for (int i = 0; i < nlen; i++){
            if (currentBody->nodes->active(i) == 1 && std::isfinite(scalarArray(i))) {
                nfile << scalarArray(i) << "\n";
            } else {
                nfile << "0" << "\n";
            }
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void DefaultVTK::writeVectorArray(MPMVectorArray& vectorArray, std::string name){
    //check length of array vs. length of open files and write to correct file
    if (pfile.is_open() && vectorArray.size() == plen){
        //write to point file
        pfile << "VECTORS " << name << " double\n";
        for (int i = 0; i < plen; i++){
            //vtk format requires x,y,z
            for (int pos = 0; pos < 3; pos++){
                if ((currentBody->points->active(i) == 1) && std::isfinite(vectorArray(i,pos))){
                    pfile << vectorArray(i,pos) << " ";
                } else {
                    pfile << "0 ";
                }
            }
            pfile << "\n";
        }
    }

    if (nfile.is_open() && vectorArray.size() == nlen){
        //write to point file
        nfile << "VECTORS " << name << " double\n";
        for (int i = 0; i < nlen; i++){
            //vtk format requires x,y,z
            for (int pos = 0; pos < 3; pos++){
                if ((currentBody->nodes->active(i) == 1) && std::isfinite(vectorArray(i,pos))){
                    nfile << vectorArray(i,pos) << " ";
                } else {
                    nfile << "0 ";
                }
            }
            nfile << "\n";
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

void DefaultVTK::writeTensorArray(MPMTensorArray& tensorArray, std::string name){
    //check length of array vs. length of open files and write to correct file
    if (pfile.is_open() && tensorArray.size() == plen){
        //write to point file
        pfile << "TENSORS " << name << " double\n";
        for (int i = 0; i < plen; i++){
            //vtk format requires x,y,z
            //brute force this one
            bool finite = true;
            for (int pos=0;pos<9;pos++){
                if (!std::isfinite(tensorArray(i,pos))){
                    finite = false;
                    break;
                }
            }
            if (currentBody->points->active(i) == 0 || !finite){
                pfile << "0 0 0\n";
                pfile << "0 0 0\n";
                pfile << "0 0 0\n";
                pfile << "\n";
            } else {
                pfile << tensorArray(i,0) << " " << tensorArray(i,1) <<  " " << tensorArray(i,2) << "\n";
                pfile << tensorArray(i,3) << " " << tensorArray(i,4) <<  " " << tensorArray(i,5) << "\n";
                pfile << tensorArray(i,6) << " " << tensorArray(i,7) <<  " " << tensorArray(i,8) << "\n";
                pfile << "\n";
            }
        }
    }

    if (nfile.is_open() && tensorArray.size() == nlen){
        //write to point file
        nfile << "TENSORS " << name << " double\n";
        for (int i = 0; i < nlen; i++){
            //vtk format requires x,y,z
            for (int pos = 0; pos < 3; pos++){
                //brute force this one
                bool finite = true;
                for (int pos=0;pos<9;pos++){
                    if (!std::isfinite(tensorArray(i,pos))){
                        finite = false;
                        break;
                    }
                }
                if (!finite){
                    nfile << "0 0 0\n";
                    nfile << "0 0 0\n";
                    nfile << "0 0 0\n";
                    nfile << "\n";
                } else {
                    nfile << tensorArray(i,0) << " " << tensorArray(i,1) <<  " " << tensorArray(i,2) << "\n";
                    nfile << tensorArray(i,3) << " " << tensorArray(i,4) <<  " " << tensorArray(i,5) << "\n";
                    nfile << tensorArray(i,6) << " " << tensorArray(i,7) <<  " " << tensorArray(i,8) << "\n";
                    nfile << "\n";
                }
            }
        }
    }
    return;
}


std::string DefaultVTK::saveState(Job* job){
    std::cerr << "NOT IMPLEMENTED! Oops." << std::endl;
    return "error";
}

/*----------------------------------------------------------------------------*/

int DefaultVTK::loadState(Job* job, std::string fullpath){
    std::cerr << "NOT IMPLEMENTED! Oops." << std::endl;
    return 0;
}

/*----------------------------------------------------------------------------*/

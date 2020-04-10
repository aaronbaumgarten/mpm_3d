//
// Created by aaron on 4/8/20.
// debug_vtk.cpp
//

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

void DebugVTK::init(Job* job){
    if (str_props.size() < 3 || fp64_props.size() < 3){
        std::cout << job->serializer->str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 6 properties defined ({sampleRate, debug_start, debug_end}, {frameDirectory, outputDirectory, outputName}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        frameDirectory = Parser::makeDirectory(str_props[0]);
        outputDirectory = Parser::makeDirectory(str_props[1]);
        outputName = str_props[2];

        sampleRate = fp64_props[0];
        t_debug_start = fp64_props[1];
        t_debug_end = fp64_props[2];
        sampledFrames = 0;
        t_last_frame = 0;

        //set buffer size to zero (see https://stackoverflow.com/questions/16605233/how-to-disable-buffering-on-a-stream)
        //pfile.rdbuf()->pubsetbuf(0, 0);
        //nfile.rdbuf()->pubsetbuf(0, 0);

        printf("Serializer properties (frameDirectory = %s, outputDirectory = %s, outputName = %s, sampleRate = %g, debug_start = %g, debug_end = %g).\n",
               frameDirectory.c_str(), outputDirectory.c_str(), outputName.c_str(), sampleRate, t_debug_start, t_debug_end);
    }

    std::cout << "Serializer Initialized." << std::endl;

    return;
}

int DebugVTK::writeFrame(Job* job){
    if ((job->t - job->t0) >= (sampledFrames/sampleRate) || sampledFrames == 0 || ((job->t > t_debug_start) && job->t < t_debug_end)){ //job->t >= sampledFrames/sampleRate){
        t_last_frame = job->t0 + sampledFrames/sampleRate;
        sampledFrames += 1;
        //write frame for each body
        for (int b=0;b<job->bodies.size();b++) {
            currentBody = job->bodies[b].get();

            //open point file
            std::stringstream ss;
            ss << "fpd." << job->bodies[b]->id << "." << job->bodies[b]->name << "." << std::setw(10) << std::setfill('0') << (sampledFrames-1) << ".vtk";
            pfilename = ss.str();
            //pfile = std::ofstream(frameDirectory+pfilename,std::ios::trunc);
            pfile.open(frameDirectory+pfilename,std::ios::trunc);

            ss.str(std::string());
            ss.clear();

            //open node file
            ss << "fnd." << job->bodies[b]->id << "." << job->bodies[b]->name << "." << std::setw(10) << std::setfill('0') << (sampledFrames-1) << ".vtk";
            nfilename = ss.str();
            //nfile = std::ofstream(frameDirectory+nfilename,std::ios::trunc);
            nfile.open(frameDirectory+nfilename,std::ios::trunc);

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

                if (job->bodies[b]->activeMaterial != 0){
                    job->bodies[b]->material->writeFrame(job, job->bodies[b].get(), job->serializer.get());
                }
                if (job->bodies[b]->activeBoundary != 0) {
                    job->bodies[b]->boundary->writeFrame(job, job->bodies[b].get(), job->serializer.get());
                }

                job->grid->writeFrame(job, job->serializer.get());
            }
            for (int c=0;c<job->contacts.size();c++) {
                job->contacts[c]->writeFrame(job, job->serializer.get());
            }

            pfile.close();
            pfile.clear();
            nfile.close();
            nfile.clear();
        }

        return 1;
    } else {
        //do not write frame
        return 0;
    }
}
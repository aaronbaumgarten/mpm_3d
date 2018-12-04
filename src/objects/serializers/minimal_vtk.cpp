//
// Created by aaron on 10/19/18.
// minimal_vtk.cpp
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

int MinimalVTK::writeFrame(Job* job){
    if ((job->t - job->t0) >= (sampledFrames/sampleRate) || sampledFrames == 0){ //job->t >= sampledFrames/sampleRate){
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

            //write out only position, velocity, pressure, and density
            //velocity
            writeVectorArray(job->bodies[b]->points->x_t, "velocity");
            //displacement
            writeVectorArray(job->bodies[b]->points->u, "displacement");
            //pressure
            Eigen::VectorXd tmpVec = Eigen::VectorXd(job->bodies[b]->points->T.size());
            for(int i=0;i<job->bodies[b]->points->T.size();i++){
                tmpVec(i) = -1.0/3.0 * job->bodies[b]->points->T[i].trace();
            }
            writeScalarArray(tmpVec,"pressure");
            //density
            tmpVec = job->bodies[b]->points->m.array() / job->bodies[b]->points->v.array();
            writeScalarArray(tmpVec,"density");

            //write nodal velocity and mass
            writeVectorArray(job->bodies[b]->nodes->x_t, "velocity");

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
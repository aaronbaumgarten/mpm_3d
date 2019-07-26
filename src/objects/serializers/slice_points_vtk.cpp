//
// Created by aaron on 7/2/19.
// slice_points_vtk.cpp
//

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

/*----------------------------------------------------------------------------*/

void SlicePointsVTK::init(Job* job){
    if (str_props.size() < 3 || fp64_props.size() < 1 || int_props.size() < 1){
        std::cout << job->serializer->str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 5 properties defined ({sampleRate}, {stride}, {frameDirectory, outputDirectory, outputName}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        frameDirectory = Parser::makeDirectory(str_props[0]);
        outputDirectory = Parser::makeDirectory(str_props[1]);
        outputName = str_props[2];

        sampleRate = fp64_props[0];
        sampledFrames = 0;
        t_last_frame = 0;

        stride = int_props[0];

        //set buffer size to zero (see https://stackoverflow.com/questions/16605233/how-to-disable-buffering-on-a-stream)
        //pfile.rdbuf()->pubsetbuf(0, 0);
        //nfile.rdbuf()->pubsetbuf(0, 0);

        printf("Serializer properties (frameDirectory = %s, outputDirectory = %s, outputName = %s, sampleRate = %g, stride = %i).\n",
               frameDirectory.c_str(), outputDirectory.c_str(), outputName.c_str(), sampleRate, stride);
    }

    std::cout << "Serializer Initialized." << std::endl;

    return;
}

int SlicePointsVTK::writeFrame(Job* job){
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
            //plen = job->bodies[b]->points->x.size();
            //use stride to get truncated length of point vector
            plen = (job->bodies[b]->points->x.size() / stride);
            nlen = job->bodies[b]->nodes->x.size();

            if (pfile.is_open()){
                pfile << "# vtk DataFile Version 3.0\n";
                pfile << "Frame: " << (sampledFrames-1) << ", Time: " << job->t << "\n";

                //job->bodies[b]->points->writeHeader(job, currentBody, job->serializer.get(), pfile, VTK);

                pfile << "ASCII\n";
                pfile << "DATASET UNSTRUCTURED_GRID\n";

                pfile << "POINTS " << plen << " double\n";
                int ii;
                for (int i=0;i<plen;i++){
                    //vtk files require x,y,z
                    ii = i*stride;
                    for (int pos = 0; pos < 3; pos++){
                        if (pos < job->bodies[b]->points->x.DIM && (job->bodies[b]->points->active(ii) != 0) && std::isfinite(job->bodies[b]->points->x(ii,pos))){
                            pfile << job->bodies[b]->points->x(ii,pos) << " ";
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


            //create temporary point velocity vector
            KinematicVectorArray tmpKVec = KinematicVectorArray(plen, job->bodies[b]->points->x_t.VECTOR_TYPE);

            //write out only position, velocity, pressure, tau_bar, gamma_dot, and density
            //velocity
            int ii;
            for (int i=0; i<plen; i++){
                ii = i*stride;
                tmpKVec(i) = job->bodies[b]->points->x_t(ii);
            }
            writeVectorArray(tmpKVec, "velocity");

            //displacement
            for (int i=0; i<plen; i++){
                ii = i*stride;
                tmpKVec(i) = job->bodies[b]->points->u(ii);
            }
            writeVectorArray(tmpKVec, "displacement");

            //pressure
            Eigen::VectorXd tmpVec = Eigen::VectorXd(plen);
            for(int i=0;i<plen;i++){
                ii = i*stride;
                tmpVec(i) = -1.0/3.0 * job->bodies[b]->points->T[ii].trace();
            }
            writeScalarArray(tmpVec,"pressure");

            //tau_bar
            for(int i=0;i<plen; i++){
                ii = i*stride;
                tmpVec(i) = (job->bodies[b]->points->T[ii] -
                             1.0/3.0*job->bodies[b]->points->T[ii].trace()*MaterialTensor::Identity()).norm() / std::sqrt(2.);
            }
            writeScalarArray(tmpVec,"tau_bar");

            //gamma_dot
            for(int i=0;i<plen; i++){
                ii = stride*i;
                tmpVec(i) = (0.5 * job->bodies[b]->points->L[ii] + 0.5 * job->bodies[b]->points->L[ii].transpose() -
                             1.0/3.0*job->bodies[b]->points->L[ii].trace()*MaterialTensor::Identity()).norm() * std::sqrt(2.);
            }
            writeScalarArray(tmpVec,"gamma_dot");

            //density
            for(int i=0; i<plen; i++){
                ii = stride*i;
                tmpVec(i) = job->bodies[b]->points->m[ii] / job->bodies[b]->points->v[ii];
            }
            writeScalarArray(tmpVec,"density");

            //write nodal velocity and density
            writeVectorArray(job->bodies[b]->nodes->x_t, "velocity");
            tmpVec = job->bodies[b]->nodes->m;
            for (int i=0; i < tmpVec.rows(); i++){
                tmpVec(i) /= job->grid->nodeVolume(job,i);
            }
            writeScalarArray(tmpVec, "density");

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

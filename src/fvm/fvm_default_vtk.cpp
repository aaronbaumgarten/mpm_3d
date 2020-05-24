//
// Created by aaron on 12/23/19.
// fvm_default_vtk.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sys/stat.h>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "mpm_objects.hpp"
#include "fvm_objects.hpp"
#include "fvm_serializers.hpp"

#include "job.hpp"

/*----------------------------------------------------------------------------*/
//initialization function
void FVMDefaultVTK::init(Job* job, FiniteVolumeDriver* driver){
    if (str_props.size() < 2){
        std::cout << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined (frameDirectory, outputName).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        frameDirectory = Parser::makeDirectory(str_props[0]);
        outputName = str_props[1];

        printf("FiniteVolumeSerializer properties (frameDirectory = %s, outputName = %s).\n",
               frameDirectory.c_str(), outputName.c_str());
    }

    std::cout << "FiniteVolumeSerializer Initialized." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//function to write FVM data to file
int FVMDefaultVTK::writeFrame(Job* job, FiniteVolumeDriver* driver){
    //always write frame when called
    sampledFrames += 1;

    //open file
    std::stringstream ss;
    ss << "fvd.0." << outputName << "." << std::setw(10) << std::setfill('0') << (sampledFrames-1) << ".vtk";
    filename = ss.str();
    //pfile = std::ofstream(frameDirectory+pfilename,std::ios::trunc);
    file.open(frameDirectory+filename,std::ios::trunc);

    ss.str(std::string());
    ss.clear();

    if (file.is_open()){
        file << "# vtk DataFile Version 3.0\n";
        file << "Frame: " << (sampledFrames-1) << ", Time: " << job->t << "\n";

        //ask grid to write VTK header for us
        driver->fluid_grid->writeHeader(file, VTK);

        //write data to file
        writeScalarArray(driver->fluid_body->rho, "density");
        writeVectorArray(driver->fluid_body->rho_x, "density_gradient");
        writeVectorArray(driver->fluid_body->p, "momentum");
        writeTensorArray(driver->fluid_body->p_x, "momentum_gradient");
        writeScalarArray(driver->fluid_body->rhoE, "energy");
        writeVectorArray(driver->fluid_body->rhoE_x, "enegry_gradient");
        writeScalarArray(driver->fluid_body->n_e, "porosity");
        writeVectorArray(driver->fluid_body->n_e_x, "porosity_gradient");

        driver->fluid_material->calculateElementTemperatures(job, driver);
        writeScalarArray(driver->fluid_body->theta, "temperature");

        driver->fluid_material->calculateElementPressures(job, driver);
        writeScalarArray(driver->fluid_body->P, "pressure");

        driver->fluid_material->calculateElementShearStresses(job, driver);
        writeTensorArray(driver->fluid_body->tau, "tau_f");

        //velocity and velocity gradients
        KinematicVectorArray u = KinematicVectorArray(driver->fluid_grid->element_count, job->JOB_TYPE);
        KinematicTensorArray u_x = driver->fluid_grid->getVelocityGradients(job, driver);
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            u[e] = driver->fluid_body->p[e]/driver->fluid_body->rho(e);
        }
        writeVectorArray(u, "velocity");
        writeTensorArray(u_x, "velocity_gradient");

        //pore density
        Eigen::VectorXd rho_f = driver->fluid_body->rho;
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            rho_f(e) /= driver->fluid_body->n_e(e);
        }
        writeScalarArray(rho_f, "true_density");

        //mach number
        Eigen::VectorXd mach = Eigen::VectorXd(driver->fluid_grid->element_count);
        for (int e=0; e<driver->fluid_grid->element_count; e++){
            mach(e) = u[e].norm()/driver->fluid_material->getSpeedOfSound(job,
                                                                          driver,
                                                                          driver->fluid_body->rho(e),
                                                                          driver->fluid_body->p[e],
                                                                          driver->fluid_body->rhoE(e),
                                                                          driver->fluid_body->n_e(e));
        }
        writeScalarArray(mach, "mach");

        //added this for debugging
        driver->solver->writeFrame(job, driver);

    } else {
        std::cerr << "Could not open frame: " << frameDirectory+filename << " !" << std::endl;
    }

    file.close();
    file.clear();

    return 1;
}

/*----------------------------------------------------------------------------*/
//functions for writing data into file
void FVMDefaultVTK::writeScalarArray(Eigen::VectorXd& scalarArray, std::string scalarName){
    //write to file
    file << "SCALARS " << scalarName << " double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < scalarArray.rows(); i++){
        if (std::isfinite(scalarArray(i))){
            file << scalarArray(i) << "\n";
        } else {
            file << "0\n";
        }
    }
    return;
}

void FVMDefaultVTK::writeVectorArray(MPMVectorArray& vectorArray, std::string vectorName){
    //write to file
    file << "VECTORS " << vectorName << " double\n";
    for (int i = 0; i < vectorArray.size(); i++){
        //vtk format requires x,y,z
        for (int pos = 0; pos < 3; pos++){
            if (std::isfinite(vectorArray(i,pos))){
                file << vectorArray(i,pos) << " ";
            } else {
                file << "0 ";
            }
        }
        file << "\n";
    }
    return;
}

void FVMDefaultVTK::writeTensorArray(MPMTensorArray& tensorArray, std::string tensorName){
    //write to point file
    file << "TENSORS " << tensorName << " double\n";
    for (int i = 0; i < tensorArray.size(); i++){
        //vtk format requires x,y,z
        //brute force this one
        bool finite = true;
        for (int pos=0;pos<9;pos++){
            if (!std::isfinite(tensorArray(i,pos))){
                finite = false;
                break;
            }
        }
        if (!finite){
            file << "0 0 0\n";
            file << "0 0 0\n";
            file << "0 0 0\n";
            file << "\n";
        } else {
            file << tensorArray(i,0) << " " << tensorArray(i,1) <<  " " << tensorArray(i,2) << "\n";
            file << tensorArray(i,3) << " " << tensorArray(i,4) <<  " " << tensorArray(i,5) << "\n";
            file << tensorArray(i,6) << " " << tensorArray(i,7) <<  " " << tensorArray(i,8) << "\n";
            file << "\n";
        }
    }
    return;
}
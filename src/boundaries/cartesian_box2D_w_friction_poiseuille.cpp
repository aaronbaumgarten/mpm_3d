//
// Created by aaron on 12/6/17.
// cartesian_box2D_w_friction_poiseuille .cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <grid.hpp>

#include "job.hpp"

#include "serializer.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "boundary.hpp"

double mu_f = 0;
double simulated_depth = 0;
double eta = 0;
Eigen::MatrixXi bcNodalMask; //store which dof to control
Eigen::VectorXd nodal_area;
Eigen::VectorXd pvec;

extern "C" void boundaryWriteFrame(Job* job, Body* body, Serializer* serializer); //write frame to serializer
extern "C" std::string boundarySaveState(Job* job, Body* body, Serializer* serializer, std::string filepath); //save state to serializer folder with returned filename
extern "C" int boundaryLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath); //read state from given full path

extern "C" void boundaryInit(Job* job, Body* body); //initialize boundary object
extern "C" void boundaryGenerateRules(Job* job, Body* body); //generate the rules given job and body state
extern "C" void boundaryApplyRules(Job* job, Body* body); //apply the rules given job and body state

/*----------------------------------------------------------------------------*/

void boundaryInit(Job* job, Body* body){
    if (job->DIM != 2){
        std::cerr << "\n FATAL ERROR: cartesian_box2D_w_fritcton.so requires DIM = 2. Given DIM = " << job->DIM << ".\nExiting..." << std::endl;
        exit(0);
    }

    if (job->grid.filename.compare("cartesian.so") != 0){
        std::cout << "\nBOUNDARY CONDITION WARNING!" << std::endl;
        std::cout << "\"cartesian_box.so\" boundary expects \"cartesian.so\" grid NOT \"" << job->grid.filename << "\"!\n" << std::endl;
    }

    //check that contact properties are set
    if ((body->boundary.fp64_props.size() < 3)){
        //need to coefficient of friction and bodies
        std::cout << body->boundary.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 3 property defined (mu_f, simulated_depth, eta).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        mu_f = body->boundary.fp64_props[0];
        simulated_depth = body->boundary.fp64_props[1];
        eta = body->boundary.fp64_props[2];
        std::cout << "Boundary properties ( mu_f = " << mu_f << ", simulated_depth = " << simulated_depth << ", eta = " << eta << ")." << std::endl;
    }

    //find bounds of box
    Eigen::VectorXd Lx = body->nodes.x.colwise().maxCoeff();

    //set bounding mask
    int len = body->nodes.x.rows();
    bcNodalMask = job->jobVectorArray<int>(len);
    bcNodalMask.setZero();
    pvec = Eigen::VectorXd(body->points.x.rows());
    nodal_area = Eigen::VectorXd(body->nodes.x.rows());

    bool is_edge = false;
    for (size_t i=0;i<len;i++){
        is_edge = false;
        for (size_t pos=0;pos<body->nodes.x.cols();pos++){
            if (body->nodes.x(i,pos) == 0){
                bcNodalMask(i,pos) = -1; //lower bound
                is_edge = true;
            } else if (body->nodes.x(i,pos) == Lx(pos)) {
                bcNodalMask(i,pos) = 1; //upper bound
                is_edge = true;
            }
        }
        for (size_t pos=0;pos<body->nodes.x.cols();pos++){
            if (is_edge && bcNodalMask(i,pos) != 1 && bcNodalMask(i,pos) != -1){
                bcNodalMask(i,pos) = 2;
            }
        }
    }

    std::cout << "Boundary Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void boundaryWriteFrame(Job* job, Body* body, Serializer* serializer){
    //write nodal mask to frame output
    Eigen::MatrixXd tmpMat = bcNodalMask.cast<double>();
    serializer->serializerWriteVectorArray(tmpMat,"bc_nodal_mask");
    return;
}

/*----------------------------------------------------------------------------*/

std::string boundarySaveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename
    std::stringstream s;
    s << "mpm_v2."  << body->name << "." << body->id << ".boundary." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    //write data
    if (ffile.is_open()) {
        ffile << "# mpm_v2 boundaries/cartesian_box.so\n"; //header
        ffile << mu_f << "\n";
        ffile << bcNodalMask.rows() << "\n"; //lines to read later
        job->jobVectorArrayToFile(bcNodalMask, ffile);

        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Boundary Saved: [" << body->name << "]." << std::endl;

    return filename;
}

/*----------------------------------------------------------------------------*/

int boundaryLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    std::string line; //read line

    std::ifstream fin(fullpath); //file to load from

    if (fin.is_open()) {
        std::getline(fin,line); //first line (header)
        std::getline(fin,line); //mu_f
        mu_f = std::stod(line);
        std::getline(fin,line); //length of file to be read
        bcNodalMask = job->jobVectorArray<int>(std::stoi(line)); //initialize vector array
        job->jobVectorArrayFromFile(bcNodalMask, fin); //read in from file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Boundary Loaded: [" << body->name << "]." << std::endl;

    return 1;
}

/*----------------------------------------------------------------------------*/

void boundaryGenerateRules(Job* job, Body* body){
    //nothing to do here
    //bc mask has already been set
    return;
}

/*----------------------------------------------------------------------------*/

void boundaryApplyRules(Job* job, Body* body){
    Eigen::VectorXd delta_momentum = job->jobVector<double>();
    double f_normal;

    for (size_t i=0;i<body->nodes.x_t.rows();i++){
        for (size_t pos=0;pos<body->nodes.x_t.cols();pos++){
            if (bcNodalMask(i,pos) == -1){
                //zero out velocity on lower boundary
                //only add friction for closing velocity
                delta_momentum(pos) = -std::min(0.0, body->nodes.mx_t(i,pos) + job->dt * body->nodes.f(i,pos));
                body->nodes.x_t(i,pos) = 0;
                body->nodes.mx_t(i,pos) = 0;
                body->nodes.f(i,pos) = 0;
            } else if (bcNodalMask(i,pos) == 1){
                //zero out velocity on upper boundary
                //only add friction for closing velocity
                delta_momentum(pos) = -std::max(0.0, body->nodes.mx_t(i,pos) + job->dt * body->nodes.f(i,pos));
                body->nodes.x_t(i,pos) = 0;
                body->nodes.mx_t(i,pos) = 0;
                body->nodes.f(i,pos) = 0;
            } else {
                delta_momentum(pos) = 0;
            }
        }
        f_normal = delta_momentum.norm() / job->dt;

        //apply friction to resulting motion
        delta_momentum = body->nodes.mx_t.row(i).transpose() + job->dt * body->nodes.f.row(i).transpose();
        if (delta_momentum.norm() > 0) {
            body->nodes.f.row(i) -=
                    std::min(delta_momentum.norm() / job->dt, f_normal * mu_f) * delta_momentum / delta_momentum.norm();
        }
    }

    //apply planar poiseuille flow rule to walls
    /*
    u_max = -dp/dx * H^2 / (8*mu)
    u_avg = 2/3 * u_max
    tau = H/2 (-dp/dx)

    -dp/dx = 12 * u_avg * mu / (H^2)
    tau = 6 * u_avg * mu / H

    f = 12 * u_avg * mu / H * da

    u_avg* = u_avg - f * dt / m
    u_avg* = f*H/(12*mu)*da
    f * H / (12 * mu * da) = u_avg - f*dt/m
    f ( H*da/(12*mu*da) + dt/m ) = u_avg
    f = m*u_avg / (m*H/(12*mu*da) + dt)
     */
    //2D sim, so we can cheat a little bit (expects depth = 1)
    body->bodyCalcNodalValues(job, nodal_area, body->points.v, Body::SET);

    if (eta != 0 && simulated_depth != 0) {
        for (size_t i=0; i<body->nodes.x_t.rows();i++){
            if (nodal_area(i) > 0) {
                delta_momentum = body->nodes.mx_t.row(i).transpose() + job->dt * body->nodes.f.row(i).transpose();
                body->nodes.f.row(i) -= delta_momentum.transpose() / (body->nodes.m(i) * simulated_depth / (12.0 * eta * nodal_area(i)) + job->dt);
            }
        }
    }

    return;
}

/*----------------------------------------------------------------------------*/


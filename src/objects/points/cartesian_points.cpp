//
// Created by aaron on 5/25/19.
// cartesian_points.cpp
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>
#include <sstream>
#include <Eigen/Core>
#include <ctime>

#include "mpm_objects.hpp"
#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "job.hpp"

#include "points.hpp"
#include "objects/bodies/bodies.hpp"

#include "registry.hpp"

/*----------------------------------------------------------------------------*/
//read point data from file
void CartesianPoints::readFromFile(Job *job, Body *body, std::string fileIN) {
    //for now, only implement 3D version
    /*
    if (job->JOB_TYPE != 3) {
        std::cerr << "ERROR: CartesianPoints is only implemented for 3D!" << std::endl;
        exit(0);
    }
     */

    file = fileIN;
    Lx = KinematicVector(job->JOB_TYPE);
    Nx = Eigen::VectorXi(job->grid->GRID_DIM);

    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::vector<std::string> params;
    std::string propName;
    std::string propValue;

    //section headers
    std::vector<std::string> headers = {"Lx","Nx","out_file","lmpp","rho","part"};

    //declare fstream from filename
    std::ifstream fin(file);

    //part registry
    Registry<Part> part_registry = Registry<Part>();

    //read config file and configure
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin, line)) {
            //remove spaces and comments from line
            line = Parser::removeComments(line);
            line = Parser::removeSpaces(line);

            //check if line matches header
            if (line.size() > 0) {
                //temporary objects for constructing job objects
                MPMObject tmp;

                line = Parser::removeBraces(line);
                line = Parser::removeQuotes(line);
                lvec = Parser::splitString(line, '=');

                //switch for headers
                switch (Parser::findStringID(headers, lvec[0])) {
                    case 0:
                        //Lx
                        lvec = Parser::splitString(lvec[1], ',');
                        for (int i = 0; i < lvec.size(); i++){
                            Lx[i] = std::stod(lvec[i]);
                        }
                        break;
                    case 1:
                        //Nx
                        lvec = Parser::splitString(lvec[1], ',');
                        for (int i = 0; i < lvec.size(); i++){
                            Nx(i) = std::stoi(lvec[i]);
                        }
                        break;
                    case 2:
                        //out_file
                        if (lvec.size() > 1){
                            out_file = lvec[1];
                        } else {
                            std::cerr << "ERROR: No output file passed to GmshPoints." << std::endl;
                        }
                        break;
                    case 3:
                        //lmpp
                        if (lvec.size() > 1){
                            lmpp = std::stoi(lvec[1]);
                        } else {
                            std::cerr << "ERROR: No lmpp value passed to GmshPoints." << std::endl;
                        }
                        break;
                    case 4:
                        //rho
                        if (lvec.size() > 1){
                            rho = std::stod(lvec[1]);
                        } else {
                            std::cerr << "ERROR: No rho value passed to GmshPoints." << std::endl;
                        }
                        break;
                    case 5:
                        //part
                        part_list.push_back(std::unique_ptr<Part>(nullptr));
                        int id = part_list.size() - 1;

                        params = {"class", "properties", "int-properties", "str-properties"};
                        std::getline(fin, line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin, line)) {
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line, '=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            part_list[id] = part_registry.get_object(propValue);
                                            break;
                                        case 1:
                                            //props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->serializer->fp64_props.push_back(std::stod(lvec[i]));
                                                tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 2:
                                            //int-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->serializer->int_props.push_back(std::stoi(lvec[i]));
                                                tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 3:
                                            //str-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->serializer->str_props.push_back(lvec[i]);
                                                tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized."
                                                      << std::endl;
                                    }
                                }
                            }
                            if (part_list[id]) {
                                //if serializer is valid
                                part_list[id]->fp64_props = tmp.fp64_props;
                                part_list[id]->int_props = tmp.int_props;
                                part_list[id]->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: Part object needs valid \'class\' defined. Exiting."
                                          << std::endl;
                                exit(0);
                            }
                            std::cout << "Part Configured: " << part_list[id]->object_name << std::endl;
                        }
                        break;
                }
            }
        }
    }

    //close file
    fin.close();

    //initialize parts
    for (int i=0; i<part_list.size(); i++){
        part_list[i]->init(job);
    }


    KinematicVector hx = KinematicVector(job->JOB_TYPE);
    KinematicVector l = KinematicVector(job->JOB_TYPE);
    //set hx
    for (int i=0; i<job->grid->GRID_DIM; i++){
        hx[i] = Lx[i]/Nx(i);
        l[i] = hx[i]/lmpp;
    }

    //calculate v
    double vTMP = 1;
    for (int i=0; i<job->grid->GRID_DIM; i++){
        vTMP *= l[i];
    }


    //create point initialization lists
    //m, v, x, x_t, a
    KinematicVectorArray xIN;
    KinematicVectorArray x_tIN;
    KinematicVector x_tTMP = KinematicVector(job->JOB_TYPE);
    KinematicVector X = KinematicVector(job->JOB_TYPE);
    x_tTMP.setZero();
    X.setZero();

    //create even distribution of points and check if in parts of body
    for (int i=0; i<Nx[0]*lmpp; i++){
        //check if 2D grid
        if (job->grid->GRID_DIM > 1){
            for (int j=0; j<Nx[1]*lmpp; j++){
                //check if 3D grid
                if (job->grid->GRID_DIM > 2){
                    for (int k=0; k<Nx[2]*lmpp; k++){
                        //sample point
                        X[0] = l[0] * (0.5 + i);
                        X[1] = l[1] * (0.5 + j);
                        X[2] = l[2] * (0.5 + k);

                        //check if any part
                        for (int p=0; p<part_list.size(); p++){
                            if (part_list[p]->encompasses(X)){
                                //add point to list
                                xIN.push_back(X);
                                x_tIN.push_back(x_tTMP);

                                break;
                            }
                        }
                    }
                } else {
                    //sample point
                    X[0] = l[0] * (0.5 + i);
                    X[1] = l[1] * (0.5 + j);

                    //check if any part
                    for (int p=0; p<part_list.size(); p++){
                        if (part_list[p]->encompasses(X)){
                            //add point to list
                            xIN.push_back(X);
                            x_tIN.push_back(x_tTMP);

                            break;
                        }
                    }
                }
            }
        } else {
            //sample point
            X[0] = l[0] * (0.5 + i);

            //check if any part
            for (int p=0; p<part_list.size(); p++){
                if (part_list[p]->encompasses(X)){
                    //add point to list
                    xIN.push_back(X);
                    x_tIN.push_back(x_tTMP);

                    break;
                }
            }
        }
    }

    //check if n-point grid quadrature points are in parts of body
    int len = xIN.size();

    //size KinematicVectors
    x = KinematicVectorArray(len, job->JOB_TYPE);
    u = KinematicVectorArray(len, job->JOB_TYPE);
    x_t = KinematicVectorArray(len, job->JOB_TYPE);
    mx_t = KinematicVectorArray(len, job->JOB_TYPE);
    b = KinematicVectorArray(len, job->JOB_TYPE);

    //size scalar vectors
    m.resize(len);
    v.resize(len);
    v0.resize(len);
    active.resize(len);
    extent.resize(len);

    //size tensor arrays
    T = MaterialTensorArray(len);
    L = KinematicTensorArray(len, job->JOB_TYPE);

    //zero out all entries to start
    x.setZero();
    u.setZero();
    x_t.setZero();
    m.setZero();
    v.setZero();
    v0.setZero();
    mx_t.setZero();
    b.setZero();
    T.setZero();
    L.setZero();
    active.setZero();
    extent.setZero();

    for (int i = 0; i < len; i++) {
        m[i] = rho*vTMP;
        v[i] = vTMP;
        x[i] = xIN[i];
        x_t[i] = x_tIN[i];
        active[i] = 1;

        //std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(x[i]).transpose() << std::endl;

        //correct volume and mass for axisymmetric simulation
        if (job->JOB_TYPE == job->JOB_AXISYM){
            m[i] *= x(i,0);
            v[i] *= x(i,0);
        }
    }

    return;
}
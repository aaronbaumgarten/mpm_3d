//
// Created by aaron on 11/30/18.
// gmsh_points.cpp
//

//
// Created by aaron on 5/14/18.
// default_points.cpp
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
//initialize point state (assumes that readFromFile has been called)
//no safety check on this, so be careful please
void GmshPoints::init(Job* job, Body* body){
    //v0 initialization
    v0 = v;

    //extent initialization (for grid integration)
    if(job->grid->GRID_DIM == 1){
        for (int i=0;i<v.rows();i++){
            extent[i] = 0.5 * v[i];
        }
    } else if (job->JOB_TYPE == job->JOB_AXISYM){
        for (int i=0;i<v.rows();i++){
            extent[i] = 0.5 * std::sqrt(v[i]/x(i,0));
        }
    } else if (job->grid->GRID_DIM == 2){
        for (int i=0;i<v.rows();i++){
            extent[i] = 0.5 * std::sqrt(v[i]);
        }
    } else if (job->grid->GRID_DIM == 3){
        for (int i = 0; i < v.rows(); i++) {
            extent[i] = 0.5 * std::cbrt(v[i]);
        }
    }

    //A matrix initialization
    int cp = 1;
    for (int i=0;i<job->grid->GRID_DIM;i++){
        cp *= 2; //square or cube
    }
    //initialize A matrix
    //for mapping position in cube to id
    //0 -> -1,-1,-1
    //1 -> +1,-1,-1
    //...
    //8 -> +1,+1,+1
    Eigen::VectorXi onoff = -1*Eigen::VectorXi::Ones(job->grid->GRID_DIM);
    A = Eigen::MatrixXi(cp, job->grid->GRID_DIM);
    for (int c=0; c<cp;c++){
        for (int i=0;i<onoff.rows();i++){
            A(c,i) = onoff(i);
        }
        for (int i=0;i<onoff.rows();i++) {
            if (onoff(i) == -1){
                onoff(i) = 1;
                break;
            } else {
                onoff(i) = -1;
            }
        }
    }

    std::cout << "Points Initialized: [" << file << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//read point data from file
void GmshPoints::readFromFile(Job *job, Body *body, std::string fileIN) {
    //for now, only implement 3D version
    if (job->JOB_TYPE != 3) {
        std::cerr << "ERROR: GmshPoints is only implemented for 3D!" << std::endl;
        exit(0);
    }

    file = fileIN;

    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::vector<std::string> params;
    std::string propName;
    std::string propValue;

    //section headers
    std::vector<std::string> headers = {"msh_file","out_file","part"};

    //declare fstream from filename
    std::ifstream fin(file);

    //list size
    int part_list_size = 0;
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
                MPMObject tmp, point_tmp, node_tmp, material_tmp, boundary_tmp;

                line = Parser::removeBraces(line);
                line = Parser::removeQuotes(line);
                lvec = Parser::splitString(line, '=');

                //switch for headers
                switch (Parser::findStringID(headers, lvec[0])) {
                    case 0:
                        //msh_file
                        if (lvec.size() > 1){
                            msh_file = lvec[1];
                        } else {
                            std::cerr << "ERROR: No .msh file passed to GmshPoints. Exiting." << std::endl;
                        }
                        break;
                    case 1:
                        //out_file
                        if (lvec.size() > 1){
                            out_file = lvec[1];
                        } else {
                            std::cerr << "ERROR: No output file passed to GmshPoints. Exiting." << std::endl;
                        }
                        break;
                    case 2:
                        //part
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
                                            part_list.push_back(
                                            job->serializer = serializer_registry.get_object(propValue);
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
                            if (job->serializer) {
                                //if serializer is valid
                                job->serializer->fp64_props = tmp.fp64_props;
                                job->serializer->int_props = tmp.int_props;
                                job->serializer->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: Serializer object needs valid \'class\' defined. Exiting."
                                          << std::endl;
                                exit(0);
                            }
                            std::cout << "Serializer Configured: " << job->serializer->object_name << std::endl;
                        }
                        break;
                }
            }
        }
    }





    int node_count, element_count;
    KinematicVectorArray x_n;
    Eigen::MatrixXi nodeIDs;
    int npe = 4;

    int len;
    std::string line; //read line
    std::vector<std::string> lvec;
    std::ifstream fin(msh_file); //file to load from
    if (fin.is_open()) {
        //if open, read lines
        std::vector<std::string> headers= {"$Nodes","$Elements"};
        while (std::getline(fin,line)) {
            switch(Parser::findStringID(headers,line)){
                case 0:
                    //nodes
                    std::getline(fin,line);
                    len = std::stoi(line); //number of nodes;

                    //nodeTags.resize(len);
                    x_n = KinematicVectorArray(len,job->JOB_TYPE);
                    node_count = len;

                    //loop to read in all nodes
                    for (int i=0;i<len;i++){
                        std::getline(fin,line);
                        lvec = Parser::splitString(line,' ');
                        //lvec[0] gives gmsh id (1-indexed)
                        x_n(i,0) = std::stod(lvec[1]); //x-coord
                        x_n(i,1) = std::stod(lvec[2]); //y-coord
                        x_n(i,2) = std::stod(lvec[3]); //z-coord
                    }
                    break;
                case 1:
                    //elements
                    if (true) {
                        std::getline(fin, line);
                        len = std::stoi(line); //number of potential elements

                        nodeIDs.resize(len, npe); //element to node map

                        int i = -1;
                        while (std::getline(fin, line)) {
                            if (line.compare("$EndElements") == 0) {
                                break;
                            }

                            lvec = Parser::splitString(line, ' ');
                            //check that element type is tetrahedron
                            if (std::stoi(lvec[1]) == 4) {
                                i++; //increment i
                                len = std::stoi(lvec[2]); //num tags
                                nodeIDs(i, 0) = std::stoi(lvec[3 + len]) - 1;
                                nodeIDs(i, 1) = std::stoi(lvec[4 + len]) - 1;
                                nodeIDs(i, 2) = std::stoi(lvec[5 + len]) - 1;
                                nodeIDs(i, 3) = std::stoi(lvec[6 + len]) - 1;
                            }
                        }

                        element_count = i + 1;
                        nodeIDs.conservativeResize(i + 1, npe);
                    }

                    break;
                default:
                    //do nothing
                    break;
            }
        }
    } else {
        std::cerr << "ERROR: Cannot open " << msh_file << "! Exiting." << std::endl;
        exit(0);
    }

    //create point initialization lists
    //m, v, x, x_t, a
    std::vector<double> vIN;
    std::vector<KinematicVector> xIN;
    std::vector<KinematicVector> x_tIN;



    std::string line;
    std::ifstream fin(file);
    std::stringstream ss;
    std::vector<std::string> s_vec;

    if (fin.is_open()) {
        std::getline(fin, line);
        int len = std::stoi(line);

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
            std::getline(fin, line);
            s_vec = Parser::splitString(line, ' ');
            if (s_vec.size() == (1 + 1 + job->DIM + job->DIM + 1)) {
                m[i] = std::stod(s_vec[0]);                 //first column is mass
                v[i] = std::stod(s_vec[1]);                 //second column is volume
                for (int pos = 0; pos < job->DIM; pos++) {
                    x(i, pos) = std::stod(s_vec[2 + pos]);    //following cols are position
                }
                for (int pos = 0; pos < job->DIM; pos++) {
                    x_t(i, pos) = std::stod(s_vec[2 + job->DIM + pos]);     //following cols are velocity
                }
                active[i] = std::stod(s_vec[2 + 2 * job->DIM]);

            } else if (s_vec.size() == (1 + 1 + job->grid->GRID_DIM + job->grid->GRID_DIM + 1)) {
                m[i] = std::stod(s_vec[0]);                 //first column is mass
                v[i] = std::stod(s_vec[1]);                 //second column is volume
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++) {
                    x(i, pos) = std::stod(s_vec[2 + pos]);    //following cols are position
                }
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++) {
                    x_t(i, pos) = std::stod(s_vec[2 + job->grid->GRID_DIM + pos]);     //following cols are velocity
                }
                active[i] = std::stod(s_vec[2 + 2 * job->grid->GRID_DIM]);

            } else {
                std::cerr << "ERROR: Unable to read line: " << file << std::endl;
                return;
            }

            //correct volume and mass for axisymmetric simulation
            if (job->JOB_TYPE == job->JOB_AXISYM){
                m[i] *= x(i,0);
                v[i] *= x(i,0);
            }
        }

    } else {
        std::cerr << "ERROR: Unable to open file: " << file << std::endl;
        return;
    }

    return;
}
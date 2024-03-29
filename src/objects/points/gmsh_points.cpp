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
//read point data from file
void GmshPoints::readFromFile(Job *job, Body *body, std::string fileIN) {
    //for now, only implement 3D version
    // May 18, 2022: attempting to implement 2D point generation
    if (job->JOB_TYPE != 3 && job->JOB_TYPE != 2) {
        std::cerr << "ERROR: GmshPoints is only implemented for 2D/3D!" << std::endl;
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
    std::vector<std::string> headers = {"msh_file","out_file","lmpp","rho","part"};

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
                        //msh_file
                        if (lvec.size() > 1){
                            msh_file = lvec[1];
                        } else {
                            std::cerr << "ERROR: No .msh file passed to GmshPoints." << std::endl;
                        }
                        break;
                    case 1:
                        //out_file
                        if (lvec.size() > 1){
                            out_file = lvec[1];
                        } else {
                            std::cerr << "ERROR: No output file passed to GmshPoints." << std::endl;
                        }
                        break;
                    case 2:
                        //lmpp
                        if (lvec.size() > 1){
                            lmpp = std::stoi(lvec[1]);
                        } else {
                            std::cerr << "ERROR: No lmpp value passed to GmshPoints." << std::endl;
                        }
                        break;
                    case 3:
                        //rho
                        if (lvec.size() > 1){
                            rho = std::stod(lvec[1]);
                        } else {
                            std::cerr << "ERROR: No rho value passed to GmshPoints." << std::endl;
                        }
                        break;
                    case 4:
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

    //read gmsh file
    int node_count, element_count;
    KinematicVectorArray x_n;
    Eigen::MatrixXi nodeIDs;
    int npe = 4;
    
    //assign npe for 2D
    if (job->JOB_TYPE == 2){
        npe = 3;
    }

    int len;
    //std::string line; //read line
    //std::vector<std::string> lvec;
    fin = std::ifstream(msh_file); //file to load from
    if (fin.is_open()) {
        //if open, read lines
        std::vector<std::string> headers= {"$MeshFormat","$Nodes","$Elements"};
        while (std::getline(fin,line)) {
            switch(Parser::findStringID(headers,line)){
                case 0:
                    //mesh format
                    std::getline(fin,line);
                    lvec = Parser::splitString(line,' ');
                    //mesh file version
                    msh_version = std::stod(lvec[0]);
                    break;
                case 1:
                    //nodes
                    std::getline(fin,line);

                    if (floor(msh_version) == 2) {
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
                    } else if (floor(msh_version) == 4){
                        lvec = Parser::splitString(line,' '); //entity_blocks num_nodes
                        len = std::stoi(lvec[1]);

                        //nodeTags.resize(len);
                        x_n = KinematicVectorArray(len,job->JOB_TYPE);
                        node_count = len;
                        //version 4 nodes are broken into blocks
                        int k = 0;
                        int block_length;
                        while (k < len){
                            //read line
                            std::getline(fin,line);
                            lvec = Parser::splitString(line, ' ');
                            //length of block
                            block_length = std::stoi(lvec[3]);
                            //skip tags
                            for (int i=0;i<block_length;i++){
                                std::getline(fin,line);
                            }
                            //read nodes in block
                            for (int i=0;i<block_length;i++){
                                std::getline(fin,line);
                                lvec = Parser::splitString(line,' ');
                                //lvec[0] gives gmsh id (1-indexed)
                                x_n(k,0) = std::stod(lvec[0]); //x-coord
                                x_n(k,1) = std::stod(lvec[1]); //y-coord
                                x_n(k,2) = std::stod(lvec[2]); //z-coord
                                //increment counter
                                k++;
                            }
                        }
                    } else {
                        std::cerr << "Unrecognized Gmsh version: " << msh_version << ". Exiting." << std::endl;
                        exit(0);
                    }
                    break;
                case 2:
                    //elements
                    if (floor(msh_version) == 2) {
                        std::getline(fin, line);
                        len = std::stoi(line); //number of potential elements

                        nodeIDs.resize(len, npe); //element to node map

                        int i = -1;
                        while (std::getline(fin, line)) {
                            if (line.compare("$EndElements") == 0) {
                                break;
                            }

                            lvec = Parser::splitString(line, ' ');
                            //check that element type is tetrahedron (or triangle in 2D)
                            if (std::stoi(lvec[1]) == 4) {
                                i++; //increment i
                                len = std::stoi(lvec[2]); //num tags
                                nodeIDs(i, 0) = std::stoi(lvec[3 + len]) - 1;
                                nodeIDs(i, 1) = std::stoi(lvec[4 + len]) - 1;
                                nodeIDs(i, 2) = std::stoi(lvec[5 + len]) - 1;
                                nodeIDs(i, 3) = std::stoi(lvec[6 + len]) - 1;
                            } else if (std::stoi(lvec[1]) == 2 && job->JOB_TYPE == 2){
                                i++; //increment i
                                len = std::stoi(lvec[2]); //num tags
                                nodeIDs(i, 0) = std::stoi(lvec[3 + len]) - 1;
                                nodeIDs(i, 1) = std::stoi(lvec[4 + len]) - 1;
                                nodeIDs(i, 2) = std::stoi(lvec[5 + len]) - 1;
                            }
                        }

                        element_count = i + 1;
                        nodeIDs.conservativeResize(i + 1, npe);

                    } else if (floor(msh_version) == 4){
                        std::getline(fin, line);
                        lvec = Parser::splitString(line,' ');
                        len = std::stoi(lvec[1]); //number of potential entities

                        nodeIDs.resize(len, npe); //element to node map

                        int i = -1;
                        int block_length, block_type;
                        while (std::getline(fin, line)){
                            if (line.compare("$EndElements") == 0){
                                break;
                            }

                            lvec = Parser::splitString(line, ' ');
                            //block info
                            block_type = std::stoi(lvec[2]);
                            block_length = std::stoi(lvec[3]);
                            if (block_type == 4) {
                                for (int k = 0; k < block_length; k++) {
                                    std::getline(fin, line);
                                    lvec = Parser::splitString(line, ' ');

                                    i++; //increment i
                                    nodeIDs(i, 0) = std::stoi(lvec[1]) - 1;
                                    nodeIDs(i, 1) = std::stoi(lvec[2]) - 1;
                                    nodeIDs(i, 2) = std::stoi(lvec[3]) - 1;
                                    nodeIDs(i, 3) = std::stoi(lvec[4]) - 1;
                                }
                            } else if (block_type == 2 && job->JOB_TYPE == 2){
                                for (int k = 0; k < block_length; k++) {
                                    std::getline(fin, line);
                                    lvec = Parser::splitString(line, ' ');

                                    i++; //increment i
                                    nodeIDs(i, 0) = std::stoi(lvec[1]) - 1;
                                    nodeIDs(i, 1) = std::stoi(lvec[2]) - 1;
                                    nodeIDs(i, 2) = std::stoi(lvec[3]) - 1;
                                }
                            } else {
                                //read through block and continue
                                for (int k=0; k<block_length; k++){
                                    std::getline(fin, line);
                                }
                            }
                        }

                        element_count = i + 1;
                        nodeIDs.conservativeResize(i + 1, npe);

                    } else {
                        std::cerr << "Unrecognized Gmsh version: " << msh_version << ". Exiting." << std::endl;
                        exit(0);
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
    KinematicVectorArray xIN = KinematicVectorArray(0, job->JOB_TYPE);
    KinematicVectorArray x_tIN = KinematicVectorArray(0, job->JOB_TYPE);
    KinematicTensor A = KinematicTensor(job->JOB_TYPE);
    KinematicTensor Ainv = KinematicTensor(job->JOB_TYPE);
    KinematicVector avec, bvec, cvec;
    avec = KinematicVector(job->JOB_TYPE);
    bvec = KinematicVector(job->JOB_TYPE);
    cvec = KinematicVector(job->JOB_TYPE);
    double v_e, vTMP;
    KinematicVector x_tTMP = KinematicVector(job->JOB_TYPE);
    x_tTMP.setZero();

    KinematicVector Xi = KinematicVector(job->JOB_TYPE);
    KinematicVector X = KinematicVector(job->JOB_TYPE);

    //check if n-point grid quadrature points are in parts of body
    for (int e = 0; e < element_count; e++){
        if (job->JOB_TYPE == 3){
            avec = x_n(nodeIDs(e,1)) - x_n(nodeIDs(e,0));
            bvec = x_n(nodeIDs(e,2)) - x_n(nodeIDs(e,0));
            cvec = x_n(nodeIDs(e,3)) - x_n(nodeIDs(e,0));

            //volume of element
            v_e = (avec.dot((bvec.cross(cvec))))/6.0;

            //mapping from xi to x
            A(0,0) = avec(0); A(1,0) = avec(1); A(2,0) = avec(2);
            A(0,1) = bvec(0); A(1,1) = bvec(1); A(2,1) = bvec(2);
            A(0,2) = cvec(0); A(1,2) = cvec(1); A(2,2) = cvec(2);

            //mapping from x to xi
            Ainv = A.inverse();

            //check middle tet
            for (int i=1; i<=lmpp; i++){
                for (int j=1; j<=lmpp; j++){
                    for (int k=1; k<=lmpp; k++){
                        if (i + j + k <= (lmpp+1)){
                            //n-point quadrature location
                            Xi(0) = i - 0.5;
                            Xi(1) = j - 0.5;
                            Xi(2) = k - 0.5;
                            Xi /= lmpp;

                            //map to x
                            X = x_n(nodeIDs(e,0)) + A*Xi;

                            //check if any part
                            for (int p=0; p<part_list.size(); p++){
                                if (part_list[p]->encompasses(X)){
                                    //add point to list
                                    xIN.push_back(X);
                                    x_tIN.push_back(x_tTMP);

                                    vTMP = 2.0*v_e/(lmpp*lmpp*lmpp);
                                    vIN.push_back(vTMP);
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            //check edge tets
            for (int i=0; i<=lmpp; i++){
                for (int j=0; j<=lmpp; j++){
                    for (int k=0; k<lmpp; k++){
                        if ((lmpp - i - j - k)%2 == 1){
                            for (int ii=-1; ii<2; ii+=2){
                                if (i == 0){
                                    //skip -1
                                    ii = 1;
                                }
                                for (int jj=-1; jj<2; jj+=2){
                                    if (j == 0){
                                        //skip -1
                                        jj = 1;
                                    }
                                    for (int kk=-1; kk<2; kk+=2){
                                        if (k == 0){
                                            //skip -1
                                            kk = 1;
                                        }
                                        //add points in octagon
                                        //n-point quadrature location
                                        Xi(0) = i + ii/4.0;
                                        Xi(1) = j + jj/4.0;
                                        Xi(2) = k + kk/4.0;
                                        Xi /= lmpp;

                                        //map to x
                                        X = x_n(nodeIDs(e,0)) + A*Xi;

                                        //check if any part
                                        for (int p=0; p<part_list.size(); p++){
                                            if (part_list[p]->encompasses(X)){
                                                //add point to list
                                                xIN.push_back(X);
                                                x_tIN.push_back(x_tTMP);

                                                vTMP = v_e/(lmpp*lmpp*lmpp);
                                                vIN.push_back(vTMP);
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else if (job->JOB_TYPE == 2){
            avec = x_n(nodeIDs(e,1)) - x_n(nodeIDs(e,0));
            bvec = x_n(nodeIDs(e,2)) - x_n(nodeIDs(e,0));

            //volume of element
            v_e = (avec.cross(bvec)).norm()/2.0;

            //mapping from xi to x
            A(0,0) = avec(0); A(1,0) = avec(1);
            A(0,1) = bvec(0); A(1,1) = bvec(1);

            //mapping from x to xi
            Ainv = A.inverse();

            //check left-handed triangles
            for (int i=1; i<=lmpp; i++){
                for (int j=1; j<=lmpp; j++){
                    if (i + j <= (lmpp+1)){
                        //centroid position
                        Xi(0) = i - 2.0/3.0;
                        Xi(1) = j - 2.0/3.0;
                        Xi /= lmpp;
                        
                        //map to x
                        X = x_n(nodeIDs(e,0)) + A*Xi;

                        //check if any part
                        for (int p=0; p<part_list.size(); p++){
                            if (part_list[p]->encompasses(X)){
                                //add point to list
                                xIN.push_back(X);
                                x_tIN.push_back(x_tTMP);

                                vTMP = v_e/(lmpp*lmpp);
                                vIN.push_back(vTMP);
                                break;
                            }
                        }
                    }
                }
            }
            
            //check right-handed triangles
            for (int i=1; i<=lmpp; i++){
                for (int j=1; j<=lmpp; j++){
                    if (i + j <= lmpp){
                        //centroid position
                        Xi(0) = i - 1.0/3.0;
                        Xi(1) = j - 1.0/3.0;
                        Xi /= lmpp;
                        
                        //map to x
                        X = x_n(nodeIDs(e,0)) + A*Xi;

                        //check if any part
                        for (int p=0; p<part_list.size(); p++){
                            if (part_list[p]->encompasses(X)){
                                //add point to list
                                xIN.push_back(X);
                                x_tIN.push_back(x_tTMP);

                                vTMP = v_e/(lmpp*lmpp);
                                vIN.push_back(vTMP);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    len = xIN.size();

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
        m[i] = rho*vIN[i];
        v[i] = vIN[i];
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

    //load output file and write ASCII data
    std::ofstream fout(out_file, std::ios::trunc);
    if (fout.is_open()){
        //write first line with number of points
        fout << x.size() << "\n";

        //write one line of data per point
        for (int i=0; i<x.size(); i++){
            fout << m[i] << " ";
            fout << v[i] << " ";
            for (int pos=0; pos<x.DIM; pos++){
                fout << x(i,pos) << " ";
            }
            for (int pos=0; pos<x_t.DIM; pos++){
                fout << x_t(i,pos) << " ";
            }
            fout << "1\n";
        }

        fout.close();
        std::cout << "Point File Written: " << out_file << "." << std::endl;
    } else {
        std::cerr << "Could not open point data file: " << out_file << " !" << std::endl;
    }

    return;
}

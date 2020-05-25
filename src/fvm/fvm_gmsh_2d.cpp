//
// Created by aaron on 12/26/19.
// fvm_gmsh_2d.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <fstream>
#include <job.hpp>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"
#include "mpm_objects.hpp"
#include "fvm_objects.hpp"
#include "fvm_grids.hpp"
#include "threadpool.hpp"

static const bool USE_NODE_NEIGHBORS = false;

void FVMGmsh2D::init(Job* job, FiniteVolumeDriver* driver){
    //assign grid dimension
    GRID_DIM = 2;
    //check that MPM grid is also 2D
    if (job->grid->GRID_DIM != 2){
        std::cerr << "ERROR: Dimension mismatch between " << job->grid->object_name << " and " << object_name << "! Exiting." << std::endl;
        exit(0);
    }

    //check size of properties passed to driver object
    if (str_props.size() < 1) {
        std::cout << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 1 property defined (filename).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        filename = str_props[0];
    }

    //call initializer for base class
    FVMGridBase::init(job, driver);

    //initialize tags and values vectors;
    std::vector<FiniteVolumeMethod::BCContainer> tmp_bc_info;


    //if integers are passed, they define boundary conditions
    if (int_props.size() == 0){
        //no bc tags given
        tmp_bc_info = std::vector<FiniteVolumeMethod::BCContainer>(0);
    } else if (int_props.size() >= 1) {
        //counter to avoid overflow error
        int fp64_iterator = 0; //already read in GRID_DIM properties

        //initialize tmp container
        tmp_bc_info = std::vector<FiniteVolumeMethod::BCContainer>(int_props.size());

        //bc_tags given
        for (int i=0; i<tmp_bc_info.size(); i++){
            tmp_bc_info[i].tag = int_props[i];
            tmp_bc_info[i].values = std::array<double, 2>();
            tmp_bc_info[i].values[0] = 0; tmp_bc_info[i].values[1] = 0;
            tmp_bc_info[i].vector = KinematicVector(job->JOB_TYPE);
            //check that bc tag is in coded list
            if (tmp_bc_info[i].tag != VELOCITY_INLET &&
                tmp_bc_info[i].tag != VELOCITY_TEMP_INLET &&
                tmp_bc_info[i].tag != VELOCITY_DENSITY_INLET &&
                tmp_bc_info[i].tag != PRESSURE_INLET &&
                tmp_bc_info[i].tag != PRESSURE_OUTLET &&
                tmp_bc_info[i].tag != DAMPED_OUTLET &&
                tmp_bc_info[i].tag != ADIABATIC_WALL &&
                tmp_bc_info[i].tag != THERMAL_WALL &&
                tmp_bc_info[i].tag != SYMMETRIC_WALL &&
                tmp_bc_info[i].tag != SUPERSONIC_INLET &&
                tmp_bc_info[i].tag != SUPERSONIC_OUTLET &&
                tmp_bc_info[i].tag != PERIODIC &&
                tmp_bc_info[i].tag != DAMPED_WALL &&
                tmp_bc_info[i].tag != STAGNATION_INLET){
                std::cerr << "ERROR: Boundary tag " << tmp_bc_info[i].tag << " not defined for FVMCartesian grid object! Exiting." << std::endl;
                exit(0);
            }

            //add required input parameters for boundary type
            switch (tmp_bc_info[i].tag){
                case VELOCITY_INLET:
                    //check that boundary condition is valid for this simulation type
                    if (driver->TYPE != FiniteVolumeDriver::ISOTHERMAL){
                        std::cerr << "ERROR: VELOCITY_INLET boundary condition requires ISOTHERMAL simulation type for stability. Exiting.";
                    } else {
                        //first GRID_DIM properties are the components of velocity at boundary
                        for (int pos=0; pos<GRID_DIM; pos++){
                            if (fp64_props.size() > fp64_iterator) {
                                tmp_bc_info[i].vector[pos] = fp64_props[fp64_iterator];
                                fp64_iterator++;
                            } else {
                                std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                                exit(0);
                            }
                        }
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : VELOCITY_INLET : u = ";
                    std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(tmp_bc_info[i].vector).transpose() << std::endl;
                    break;

                case VELOCITY_TEMP_INLET:
                    //check that boundary condition is valid for this simulation type
                    if (driver->TYPE != FiniteVolumeDriver::THERMAL){
                        std::cerr << "ERROR: VELOCITY_TEMP_INLET boundary condition requires THERMAL simulation type. Exiting.";
                    } else {
                        //first GRID_DIM properties are the components of velocity at boundary
                        for (int pos=0; pos<GRID_DIM; pos++){
                            if (fp64_props.size() > fp64_iterator) {
                                tmp_bc_info[i].vector[pos] = fp64_props[fp64_iterator];
                                fp64_iterator++;
                            } else {
                                std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                                exit(0);
                            }
                        }
                        //next property is temperature
                        if (fp64_props.size() > fp64_iterator) {
                            tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                            fp64_iterator++;
                        } else {
                            std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                            exit(0);
                        }
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : VELOCITY_TEMP_INLET : u = ";
                    std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(tmp_bc_info[i].vector).transpose();
                    std::cout << ", T = " << tmp_bc_info[i].values[0] << std::endl;
                    break;

                case VELOCITY_DENSITY_INLET:
                    //first GRID_DIM properties are the components of velocity at boundary
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (fp64_props.size() > fp64_iterator) {
                            tmp_bc_info[i].vector[pos] = fp64_props[fp64_iterator];
                            fp64_iterator++;
                        } else {
                            std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                            exit(0);
                        }
                    }
                    //next property is density
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : VELOCITY_DENSITY_INLET : u = ";
                    std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(tmp_bc_info[i].vector).transpose();
                    std::cout << ", rho = " << tmp_bc_info[i].values[0] << std::endl;
                    break;

                case PRESSURE_INLET:
                    //first property is density
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }
                    //second property is temperature
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[1] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : PRESSURE_INLET : P = " << tmp_bc_info[i].values[0];
                    std::cout << ", T = " << tmp_bc_info[i].values[1] << std::endl;
                    break;

                case PRESSURE_OUTLET:
                    //first property is density
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }
                    //second property is temperature
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[1] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : PRESSURE_OUTLET : P = " << tmp_bc_info[i].values[0];
                    std::cout << ", T* = " << tmp_bc_info[i].values[1] << std::endl;
                    break;

                case DAMPED_OUTLET:
                    //first property is density
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }
                    //second property is temperature
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[1] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : DAMPED_OUTLET : P = " << tmp_bc_info[i].values[0];
                    std::cout << ", T* = " << tmp_bc_info[i].values[1] << std::endl;
                    break;

                case ADIABATIC_WALL:
                    //for readability, take in one unused property
                    fp64_iterator++;

                    //print boundary condition info
                    std::cout << " - " << i << " : ADIABATIC_WALL" << std::endl;
                    break;

                case THERMAL_WALL:
                    //one property, temperature
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : THERMAL_WALL : T = " << tmp_bc_info[i].values[0] << std::endl;
                    break;

                case SYMMETRIC_WALL:
                    //for readability, take in one unused property
                    fp64_iterator++;

                    //print boundary condition info
                    std::cout << " - " << i << " : SYMMETRIC_WALL" << std::endl;
                    break;

                case SUPERSONIC_INLET:
                    //first GRID_DIM properties are the components of velocity at boundary
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (fp64_props.size() > fp64_iterator) {
                            tmp_bc_info[i].vector[pos] = fp64_props[fp64_iterator];
                            fp64_iterator++;
                        } else {
                            std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                            exit(0);
                        }
                    }
                    //next property is density
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }
                    //final property is temperature
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[1] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : SUPERSONIC_INLET : u = ";
                    std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(tmp_bc_info[i].vector).transpose();
                    std::cout << ", rho = " << tmp_bc_info[i].values[0];
                    std::cout << ", T = " << tmp_bc_info[i].values[1] << std::endl;
                    break;

                case SUPERSONIC_OUTLET:
                    //for readability, take in one unused property
                    fp64_iterator++;

                    //print boundary condition info
                    std::cout << " - " << i << " : SUPERSONIC_OUTLET" << std::endl;
                    break;

                case PERIODIC:
                    //for readability, take in one unused property
                    fp64_iterator++;

                    //print boundary condition info
                    std::cout << " - " << i << " : PERIODIC" << std::endl;
                    std::cerr << "ERROR: FVMGmsh2D does not allow PERIODIC boundaries! Exiting." << std::endl;
                    exit(0);
                    break;

                case DAMPED_WALL:
                    //one property, damping coefficient
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : DAMPED_WALL : nu = " << tmp_bc_info[i].values[0] << std::endl;
                    break;

                case STAGNATION_INLET:
                    //two properties: P^t and T^t
                    if (fp64_props.size() > fp64_iterator+1){
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                        tmp_bc_info[i].values[1] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : STAGNATION_INLET : P^t = " << tmp_bc_info[i].values[0];
                    std::cout << ", T^t = " << tmp_bc_info[i].values[1] << std::endl;
                    break;

                default:
                    //do nothing
                    break;
            }

        }
    }

    if (job->JOB_TYPE != 2){
        std::cerr << "WARNING: Cannot use FVMGmsh2D with JOB_TYPE of " << job->JOB_TYPE << "!" << std::endl;
    }

    //temporary vector listing faces on boundary
    std::vector<std::array<int,3>> bounding_faces; //id, node1, node2

    double msh_version;
    int len;
    std::string line; //read line
    std::vector<std::string> lvec;
    std::ifstream fin(filename); //file to load from
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
                        x_n = KinematicVectorArray(len, job->JOB_TYPE);
                        node_count = len;

                        //loop to read in all nodes
                        for (int i = 0; i < len; i++) {
                            std::getline(fin, line);
                            lvec = Parser::splitString(line, ' ');
                            //lvec[0] gives gmsh id (1-indexed)
                            x_n(i, 0) = std::stod(lvec[1]); //x-coord
                            x_n(i, 1) = std::stod(lvec[2]); //y-coord
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
                        v_e.resize(len); //element volume

                        int i = -1;
                        while (std::getline(fin, line)) {
                            if (line.compare("$EndElements") == 0) {
                                break;
                            }

                            lvec = Parser::splitString(line, ' ');
                            //check that element type is triangle
                            if (std::stoi(lvec[1]) == 2) {
                                i++; //increment i
                                len = std::stoi(lvec[2]); //num tags
                                nodeIDs(i, 0) = std::stoi(lvec[3 + len]) - 1;
                                nodeIDs(i, 1) = std::stoi(lvec[4 + len]) - 1;
                                nodeIDs(i, 2) = std::stoi(lvec[5 + len]) - 1;
                            } else if (std::stoi(lvec[1]) == 1) {
                                //line element
                                len = std::stoi(lvec[2]); //num tags
                                std::array<int,3> tmp_array = std::array<int,3>();
                                tmp_array[0] = std::stoi(lvec[3]);
                                tmp_array[1] = std::stoi(lvec[3 + len]) - 1;
                                tmp_array[2] = std::stoi(lvec[4 + len]) - 1;
                                //add to list of faces on boundary
                                bounding_faces.push_back(tmp_array);
                            }
                        }

                        element_count = i + 1;
                        nodeIDs.conservativeResize(i + 1, 3);
                    }  else if (floor(msh_version) == 4){
                        std::getline(fin, line);
                        lvec = Parser::splitString(line,' ');
                        len = std::stoi(lvec[1]); //number of potential entities

                        nodeIDs.resize(len, npe); //element to node map

                        int i = -1;
                        int block_length, block_type, tag;
                        while (std::getline(fin, line)){
                            if (line.compare("$EndElements") == 0){
                                break;
                            }

                            lvec = Parser::splitString(line, ' ');
                            //block info
                            tag = std::stoi(lvec[1]);
                            block_type = std::stoi(lvec[2]);
                            block_length = std::stoi(lvec[3]);
                            //check that entity is triangle
                            if (block_type == 2) {
                                for (int k = 0; k < block_length; k++) {
                                    std::getline(fin, line);
                                    lvec = Parser::splitString(line, ' ');

                                    i++; //increment i
                                    nodeIDs(i, 0) = std::stoi(lvec[1]) - 1;
                                    nodeIDs(i, 1) = std::stoi(lvec[2]) - 1;
                                    nodeIDs(i, 2) = std::stoi(lvec[3]) - 1;
                                }
                            } else if (block_type == 1) {
                                //if entity is a line, it defines a set of physical tags for nodes
                                for (int k=0; k<block_length; k++) {
                                    std::array<int, 3> tmp_array = std::array<int, 3>();
                                    tmp_array[0] = tag;
                                    tmp_array[1] = std::stoi(lvec[1]) - 1;
                                    tmp_array[2] = std::stoi(lvec[2]) - 1;
                                    bounding_faces.push_back(tmp_array);
                                }
                            } else {
                                //read through block and continue
                                for (int k = 0; k < block_length; k++) {
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
        std::cerr << "ERROR: Cannot open " << filename << "! Exiting." << std::endl;
        exit(0);
    }

    //element volume
    v_e.resize(nodeIDs.rows());
    x_e = KinematicVectorArray(element_count, job->JOB_TYPE);
    //herons formula
    double a, b, c, s;
    for (int e=0; e<nodeIDs.rows(); e++){
        a = (x_n(nodeIDs(e,0)) - x_n(nodeIDs(e,1))).norm();
        b = (x_n(nodeIDs(e,1)) - x_n(nodeIDs(e,2))).norm();
        c = (x_n(nodeIDs(e,2)) - x_n(nodeIDs(e,0))).norm();
        s = (a+b+c)/2.0;

        v_e(e) = std::sqrt(s*(s-a)*(s-b)*(s-c));

        //element centroid
        x_e[e] = (x_n(nodeIDs(e,0)) + x_n(nodeIDs(e,1)) + x_n(nodeIDs(e,2)))/3.0;
    }

    //define faces from mesh (for now, just face->node map, element->face map, and face->element map)
    element_faces = std::vector<std::vector<int>>(element_count);
    std::vector<std::vector<std::array<int,2>>> tmp_node_to_face_map = std::vector<std::vector<std::array<int,2>>>(node_count);
    //[node id]: [other node, face_id]
    int f_tmp=0;
    int n0, n1, n2;
    bool face_defined;
    for (int e=0; e<element_count; e++){
        //loop over elements
        n0 = nodeIDs(e,0);
        n1 = nodeIDs(e,1);
        n2 = nodeIDs(e,2);
        //FACE 1 <0, 1>
        //check if already defined
        face_defined = false;
        for (int i=0; i<tmp_node_to_face_map[n0].size(); i++){
            if (n1 == tmp_node_to_face_map[n0][i][0]){
                //face is defined
                face_defined = true;
                f_tmp = tmp_node_to_face_map[n0][i][1]; //store face id
                //exit loop
                break;
            }
        }
        //if face is not defined, then define it
        if (!face_defined){
            std::array<int,2> tmp_array = {n0, n1};
            face_nodes.push_back(tmp_array);                //add nodes to define face
            f_tmp = face_nodes.size() - 1;
            tmp_array[1] = f_tmp;
            tmp_node_to_face_map[n1].push_back(tmp_array);  //add to temporary container so we don't duplicate faces
            tmp_array[0] = n1;
            tmp_node_to_face_map[n0].push_back(tmp_array);  //make sure container is symmetric

            //add element to face map
            tmp_array[0] = e; tmp_array[1] = -1;
            face_elements.push_back(tmp_array);
        } else {
            //add element to face map
            face_elements[f_tmp][1] = e;
        }

        //add face to element list
        element_faces[e].push_back(f_tmp);

        //FACE 2 <1, 2>
        //check if already defined
        face_defined = false;
        for (int i=0; i<tmp_node_to_face_map[n1].size(); i++){
            if (n2 == tmp_node_to_face_map[n1][i][0]){
                //face is defined
                face_defined = true;
                f_tmp = tmp_node_to_face_map[n1][i][1]; //store face id
                //exit loop
                break;
            }
        }
        //if face is not defined, then define it
        if (!face_defined){
            std::array<int,2> tmp_array = {n1, n2};
            face_nodes.push_back(tmp_array);                //add nodes to define face
            f_tmp = face_nodes.size() - 1;
            tmp_array[1] = f_tmp;
            tmp_node_to_face_map[n2].push_back(tmp_array);  //add to temporary container so we don't duplicate faces
            tmp_array[0] = n2;
            tmp_node_to_face_map[n1].push_back(tmp_array);  //make sure container is symmetric

            //add element to face map
            tmp_array[0] = e; tmp_array[1] = -1;
            face_elements.push_back(tmp_array);
        } else {
            //add element to face map
            face_elements[f_tmp][1] = e;
        }

        //add face to element list
        element_faces[e].push_back(f_tmp);

        //FACE 3 <2, 0>
        //check if already defined
        face_defined = false;
        for (int i=0; i<tmp_node_to_face_map[n2].size(); i++){
            if (n0 == tmp_node_to_face_map[n2][i][0]){
                //face is defined
                face_defined = true;
                f_tmp = tmp_node_to_face_map[n2][i][1]; //store face id
                //exit loop
                break;
            }
        }
        //if face is not defined, then define it
        if (!face_defined){
            std::array<int,2> tmp_array = {n2, n0};
            face_nodes.push_back(tmp_array);                //add nodes to define face
            f_tmp = face_nodes.size() - 1;
            tmp_array[1] = f_tmp;
            tmp_node_to_face_map[n0].push_back(tmp_array);  //add to temporary container so we don't duplicate faces
            tmp_array[0] = n0;
            tmp_node_to_face_map[n2].push_back(tmp_array);  //make sure container is symmetric

            //add element to face map
            tmp_array[0] = e; tmp_array[1] = -1;
            face_elements.push_back(tmp_array);
        } else {
            //add element to face map
            face_elements[f_tmp][1] = e;
        }

        //add face to element list
        element_faces[e].push_back(f_tmp);
    }

    //consistency check
    for (int e=0; e<element_count; e++){
        if (element_faces[e].size() != 3){
            std::cout << "ERROR! Element " << e << " has incorrect number of faces: " << element_faces[e].size() << " != 3." << std::endl;
        }
    }
    if (face_nodes.size() != face_elements.size()){
        //something is wrong
        std::cerr << "ERROR! Somehow face_nodes and face_elements maps are different sizes!" << std::endl;
        exit(0);
    }

    //finish creating face definitions (face areas, normals, and centroids)
    face_count = face_nodes.size();
    face_normals = KinematicVectorArray(face_count, job->JOB_TYPE);
    face_areas = Eigen::VectorXd(face_count);
    x_f = KinematicVectorArray(face_count, job->JOB_TYPE);
    for (int f=0; f<face_count; f++){
        //node ids
        n0 = face_nodes[f][0];
        n1 = face_nodes[f][1];

        //face area = length of face
        face_areas(f) = (x_n[n0] - x_n[n1]).norm();

        //face centroid
        x_f[f] = (x_n[n0] + x_n[n1])/2.0;

        //face normal
        face_normals(f, 0) = -x_n(n0, 1) + x_n(n1, 1);
        face_normals(f, 1) = x_n(n0, 0) - x_n(n1, 0);
        face_normals[f] /= face_normals[f].norm();

        //ensure face normal consistent with element definitions
        if ((x_e[face_elements[f][0]] - x_f[f]).dot(face_normals[f]) > 0){
            //normal is pointing into element A, so it must be switched: A -|-> B
            face_normals[f] *= -1.0;
        }
    }

    //define element neighbors map
    element_neighbors = std::vector<std::vector<int>>(element_count);
    if (USE_NODE_NEIGHBORS) {
        //use all connected nodes
        std::vector<std::vector<int>> tmp_node_to_element_map = std::vector<std::vector<int>>(node_count);
        for (int e = 0; e < element_count; e++) {
            //fill in list of elements associated with each node
            n0 = nodeIDs(e, 0);
            n1 = nodeIDs(e, 1);
            n2 = nodeIDs(e, 2);

            tmp_node_to_element_map[n0].push_back(e);
            tmp_node_to_element_map[n1].push_back(e);
            tmp_node_to_element_map[n2].push_back(e);
        }
        int e_tmp;
        bool in_list;
        for (int e = 0; e < element_count; e++) {
            //loop over elements
            for (int j = 0; j < npe; j++) {
                //for each node
                n0 = nodeIDs(e, j);
                //check if elements associated with that node are already in the neighbor list
                for (int i = 0; i < tmp_node_to_element_map[n0].size(); i++) {
                    e_tmp = tmp_node_to_element_map[n0][i];
                    in_list = false;
                    for (int ii = 0; i < element_neighbors[e].size(); i++) {
                        if (element_neighbors[e][ii] == e_tmp) {
                            in_list = true;
                            break;
                        }
                    }
                    if (!in_list) {
                        //if element is not in neighbor list, add it
                        element_neighbors[e].push_back(e_tmp);
                    }
                }
            }
        }
    } else {
        //use face neighbors only
        int e_tmp;
        bool in_list;
        for (int f = 0; f < face_count; f++) {
            //loop over faces
            if (face_elements[f][0] > -1 && face_elements[f][1] > -1){
                //these can only show up once
                element_neighbors[face_elements[f][0]].push_back(face_elements[f][1]);
                element_neighbors[face_elements[f][1]].push_back(face_elements[f][0]);
            }
        }
    }

    //define which boundary is associated with each face
    bc_info = std::vector<FiniteVolumeMethod::BCContainer>(face_count);
    //initalize bc_info
    for (int f=0; f<face_count; f++){
        bc_info[f].tag = -1;
        bc_info[f].values = std::array<double,2>();
        bc_info[f].values[0] = 0;
        bc_info[f].values[1] = 0;
        bc_info[f].vector = KinematicVector(job->JOB_TYPE);
        bc_info[f].vector.setZero();
    }

    /*
    std::cout << "faces:" << std::endl;
    for (int i=0; i<bounding_faces.size(); i++){
        std::cout << "[" << i << "]: " << bounding_faces[i][0] << std::endl;
    }
     */

    //loop over stored boundary faces
    int which_boundary;
    for (int i=0; i<bounding_faces.size(); i++){
        which_boundary = bounding_faces[i][0];
        n0 = bounding_faces[i][1];
        n1 = bounding_faces[i][2];
        //check faces associated with n0 for n1 from node to face map create earlier
        for (int ii=0; ii<tmp_node_to_face_map[n0].size(); ii++){
            if (n1 == tmp_node_to_face_map[n0][ii][0]){
                //found face
                f_tmp = tmp_node_to_face_map[n0][ii][1];
                break;
            }
        }
        if (which_boundary >= tmp_bc_info.size()){
            //if no bc defined by input file, assume dirichlet
            //no bc tags given
            if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
                bc_info[f_tmp].tag = FiniteVolumeGrid::VELOCITY_INLET;
                bc_info[f_tmp].values = std::array<double, 2>();
                bc_info[f_tmp].vector = KinematicVector(job->JOB_TYPE);
                bc_info[f_tmp].vector.setZero();
            } else if (driver->TYPE == FiniteVolumeDriver::THERMAL) {
                bc_info[f_tmp].tag = FiniteVolumeGrid::ADIABATIC_WALL;
                bc_info[f_tmp].values = std::array<double, 2>();
                bc_info[f_tmp].vector = KinematicVector(job->JOB_TYPE);
                bc_info[f_tmp].vector.setZero();
            } else {
                std::cerr << "ERROR: FVMGmsh2D is not implemented for FiniteVolumeDriver TYPE " << driver->TYPE << "! Exiting.";
                exit(0);
            }
        } else {
            //copy associated input bc info to face vector
            bc_info[f_tmp] = tmp_bc_info[which_boundary];
        }
    }

    //penultimately, need to form least squares system of equations
    //construct A and b for least squares calculations
    A_e = std::vector<Eigen::MatrixXd>(element_count);
    A_inv = std::vector<Eigen::MatrixXd>(element_count);
    b_e = std::vector<Eigen::VectorXd>(element_count);
    int length_of_A; //including element neighbors and dirichlet faces
    //initialize matrix A and vector B
    KinematicVector x_0, x;
    double tmp_dif;
    for (int e = 0; e<element_count; e++) {
        //size A and b for element based on number of neighboring cells and dirichlet faces
        length_of_A = element_neighbors[e].size();

        //get element centroid
        x_0 = getElementCentroid(job, e);
        A_e[e] = Eigen::MatrixXd(length_of_A,GRID_DIM);
        b_e[e] = Eigen::VectorXd(length_of_A);
        //create system of equations
        for (int ii = 0; ii < element_neighbors[e].size(); ii++) {
            x = getElementCentroid(job, element_neighbors[e][ii]);
            for (int pos = 0; pos < GRID_DIM; pos++) {
                tmp_dif = x[pos] - x_0[pos];
                A_e[e](ii, pos) = x[pos] - x_0[pos];
            }
        }

        //psuedo inverse
        A_inv[e] = Eigen::MatrixXd(GRID_DIM,A_e[e].rows());
        Eigen::MatrixXd AtA_inv = Eigen::MatrixXd(GRID_DIM,GRID_DIM);
        AtA_inv.setZero();
        Eigen::MatrixXd AtA = A_e[e].transpose()*A_e[e];
        Eigen::VectorXd a_tmp;
        Eigen::VectorXd e_vec = Eigen::VectorXd(AtA.rows());
        for (int ii=0; ii<AtA.rows(); ii++){
            e_vec.setZero();
            e_vec(ii) = 1;
            //a_tmp = AtA.householderQr().solve(e_vec);
            a_tmp = AtA.fullPivHouseholderQr().solve(e_vec);
            for (int jj=0; jj<GRID_DIM; jj++){
                AtA_inv(jj,ii) = a_tmp(jj);
            }
        }
        A_inv[e] = AtA_inv*A_e[e].transpose();
    }

    //lastly, quadrature rule
    if (driver->ORDER == 1 || USE_REDUCED_QUADRATURE){
        qpe = 1; //quad points per element
        qpf = 1; //quad points per face
    } else if (driver->ORDER == 2){
        qpe = 3;
        qpf = 2;
    } else {
        std::cerr << "ERROR: FVMCartesian not defined for simulation ORDER > 2." << std::endl;
        qpe = 3;
        qpf = 2;
    }
    ext_quad_count = face_count*qpf;
    int_quad_count = element_count*qpe;

    //define face, element, and quadrature centroids
    x_q = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE); //volume integrals first, then surface
    w_q = Eigen::VectorXd(int_quad_count + ext_quad_count);
    q_b = std::vector<bool>(int_quad_count + ext_quad_count);

    //initialize porosity and solid velocity fields
    n_q = Eigen::VectorXd(int_quad_count + ext_quad_count);
    n_q.setConstant(1.0);
    gradn_q = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE);
    gradn_q.setZero();
    v_sq = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE);
    v_sq.setZero();

    //quadrature points are default NOT on boundary
    for (int i=0; i<int_quad_count+ext_quad_count; i++){
        q_b[i] = false;
    }

    //quadrature points
    double offset_factor = 1.0/sqrt(3.0);
    KinematicVector x1, x2, x3;

    //internal quadrature rule
    for (int e=0; e<element_count; e++){
        if (qpe == 1){
            //quad point collocated with element centroid
            x_q[e*qpe] = x_e[e];
            w_q(e*qpe) = v_e[e];
        } else if (qpe == 3){
            //quad points offset from center
            x_q[e*qpe + 0] = 0.5*(x_e[e] + x_n[nodeIDs(e,0)]);
            w_q(e*qpe + 0) = v_e[e]/qpe;
            x_q[e*qpe + 1] = 0.5*(x_e[e] + x_n[nodeIDs(e,1)]);
            w_q(e*qpe + 1) = v_e[e]/qpe;
            x_q[e*qpe + 2] = 0.5*(x_e[e] + x_n[nodeIDs(e,2)]);
            w_q(e*qpe + 2) = v_e[e]/qpe;
        }
    }

    //external quadrature rule
    for (int f=0; f<face_count; f++){
        for (int q=0; q<qpf; q++){
            if (bc_info[f].tag > -1 && bc_info[f].tag != PERIODIC) {
                //face is bounding face and quadrature point should be included in boundary integrals
                q_b[int_quad_count + f * qpf + q] = true;
            }
        }

        if (qpf == 1){
            //quad point collocated with face centroid
            x_q[int_quad_count + f * qpf] = x_f[f];
            w_q(int_quad_count + f * qpf) = face_areas(f);
        } else if (qpf == 2){
            //quad points offset from center
            x_q[int_quad_count + f*qpf + 0] = x_f[f] + offset_factor*(x_n[face_nodes[f][0]] - x_f[f]);
            w_q(int_quad_count + f*qpf + 0) = getFaceArea(f)/qpf;
            x_q[int_quad_count + f*qpf + 1] = x_f[f] + offset_factor*(x_n[face_nodes[f][1]] - x_f[f]);
            w_q(int_quad_count + f*qpf + 1) = getFaceArea(f)/qpf;
        }
    }

    //initialize grid mappings
    generateMappings(job, driver);

    std::cout << "FiniteVolumeGrid initialized." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//function to write header of VTK file
void FVMGmsh2D::writeHeader(std::ofstream& file, int SPEC){
    if (SPEC != Serializer::VTK){
        std::cerr << "ERROR: Unknown file SPEC in writeHeader: " << SPEC  << "! Exiting." << std::endl;
        exit(0);
    }

    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "POINTS " << node_count << " double\n";
    for (int i=0;i<node_count;i++){
        //vtk requires 3D position
        file << x_n(i,0) << " " << x_n(i,1) << " " << 0 << "\n";
    }
    file << "CELLS " << element_count << " " << 4*element_count << "\n";
    for (int e=0;e<element_count;e++){
        file << "3 " << nodeIDs(e,0) << " " << nodeIDs(e,1) << " " << nodeIDs(e,2) << "\n";
    }

    file << "CELL_TYPES " << element_count << "\n";
    for (int e=0;e<element_count;e++){
        file << "5\n";
    }
    file << "CELL_DATA " << element_count << "\n";
    return;
}

/*----------------------------------------------------------------------------*/
//reconstruct momentum field
void FVMGmsh2D::constructMomentumField(Job* job, FiniteVolumeDriver* driver){
    if (driver->ORDER == 1){
        driver->fluid_body->p_x.setZero(); //let momentum be constant within an element
    } else if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMGmsh2D does not currently implement higher ORDER momentum reconstruction."
                      << std::endl;
        }
        //initialize holder matrices
        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        KinematicVector p_0, p, p_max, p_min, min_dif, tmp_val;
        double rho_0;
        for (int e = 0; e<element_count; e++){
            //least squares fit of u_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            p_0 = driver->fluid_body->p[e] / driver->fluid_body->n_e(e);
            rho_0 = driver->fluid_body->rho(e) / driver->fluid_body->n_e(e);
            p_max = p_0; p_min = p_0;
            for (int mom_index = 0; mom_index<GRID_DIM; mom_index++){
                //create system of equations
                for (int ii=0; ii<element_neighbors[e].size(); ii++){
                    //ill in b vector
                    p = driver->fluid_body->p[element_neighbors[e][ii]] / driver->fluid_body->n_e(element_neighbors[e][ii]);
                    b_e[e](ii) = p[mom_index] - p_0[mom_index];

                    //update maximum and minimum velocities
                    if (p[mom_index] > p_max[mom_index]){
                        p_max[mom_index] = p[mom_index];
                    } else if (p[mom_index] < p_min[mom_index]){
                        p_min[mom_index] = p[mom_index];
                    }
                }

                //solve for component of gradient
                sol = A_inv[e]*b_e[e];
                for (int pos = 0; pos<GRID_DIM; pos++){
                    driver->fluid_body->p_x(e, mom_index, pos) = sol(pos);
                }
            }

            //limit gradient to ensure monotonicity
            //check max, min at each node
            for (int n=0; n<npe; n++){
                //p(x) = grad(p) * (x - x_0)
                tmp_val = driver->fluid_body->p_x[e]*(x_n[nodeIDs(e,n)] - x_e[e]);
                //check for both overshoot and undershoot in each direction
                for (int mom_index=0; mom_index < GRID_DIM; mom_index++) {
                    if (tmp_val[mom_index] > (p_max[mom_index] - p_0[mom_index])){
                        //limit gradient
                        for (int pos=0; pos<GRID_DIM; pos++) {
                            driver->fluid_body->p_x(e, mom_index, pos) *= (p_max[mom_index] - p_0[mom_index])/tmp_val[mom_index];
                        }
                    } else if (tmp_val[mom_index] < (p_min[mom_index] - p_0[mom_index])){
                        //limit gradient
                        for (int pos=0; pos<GRID_DIM; pos++) {
                            driver->fluid_body->p_x(e, mom_index, pos) *= (p_min[mom_index] - p_0[mom_index])/tmp_val[mom_index];
                        }
                    }
                }
            }

            //transform true momentum gradient into effective momentum gradient
            driver->fluid_body->p_x[e] *= driver->fluid_body->n_e(e);
        }
    }
    return;
}

void FVMGmsh2D::constructDensityField(Job* job, FiniteVolumeDriver* driver){
    if (driver->ORDER == 1){
        driver->fluid_body->rho_x.setZero(); //let density be constant within an element
    } if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMCartesian does not currently implement higher ORDER momentum reconstruction."
                      << std::endl;
        }
        //initialize least squares fields
        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        double rho_0, rho, rho_max, rho_min, min_dif;
        double tmp_val;
        for (int e = 0; e<element_count; e++){
            //least squares fit of rho_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            rho_0 = driver->fluid_body->rho(e) / driver->fluid_body->n_e(e);
            rho_max = rho_0;
            rho_min = rho_0;

            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                //fill in b vector
                rho = driver->fluid_body->rho(element_neighbors[e][ii]) / driver->fluid_body->n_e(element_neighbors[e][ii]);
                b_e[e](ii) = rho - rho_0;

                //update maximum and minimum velocities
                if (rho > rho_max){
                    rho_max = rho;
                } else if (rho < rho_min){
                    rho_min = rho;
                }
            }

            //solve for components of gradient
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                driver->fluid_body->rho_x(e, pos) = sol(pos);
            }

            //limit gradient to ensure monotonicity
            //check max, min at each node
            for (int n=0; n<npe; n++){
                //p(x) = grad(p) * (x - x_0)
                tmp_val = driver->fluid_body->rho_x[e].dot(x_n[nodeIDs(e,n)] - x_e[e]);
                //check for both overshoot and undershoot
                if (tmp_val > (rho_max - rho_0)){
                    //limit gradient
                    for (int pos=0; pos<GRID_DIM; pos++) {
                        driver->fluid_body->rho_x[e] *= (rho_max - rho_0)/tmp_val;
                    }
                } else if (tmp_val < (rho_min - rho_0)){
                    //limit gradient
                    for (int pos=0; pos<GRID_DIM; pos++) {
                        driver->fluid_body->rho_x[e] *= (rho_min - rho_0)/tmp_val;
                    }
                }
            }

            //transform true density gradient into effective density gradient
            driver->fluid_body->rho_x[e] *= driver->fluid_body->n_e(e);
        }
    }
    return;
}


void FVMGmsh2D::constructEnergyField(Job* job, FiniteVolumeDriver* driver){
    if (driver->ORDER == 1){
        driver->fluid_body->rhoE_x.setZero(); //let density be constant within an element
    } if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMCartesian does not currently implement higher ORDER momentum reconstruction."
                      << std::endl;
        }
        //initialize least squares fields
        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        double rhoE_0, rhoE, rhoE_max, rhoE_min, min_dif;
        double tmp_val;
        for (int e = 0; e<element_count; e++){
            //least squares fit of rho_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            rhoE_0 = driver->fluid_body->rhoE(e) / driver->fluid_body->n_e(e);
            rhoE_max = rhoE_0;
            rhoE_min = rhoE_0;

            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                //fill in b vector
                rhoE = driver->fluid_body->rhoE(element_neighbors[e][ii]) / driver->fluid_body->n_e(element_neighbors[e][ii]);
                b_e[e](ii) = rhoE - rhoE_0;

                //update maximum and minimum velocities
                if (rhoE > rhoE_max){
                    rhoE_max = rhoE;
                } else if (rhoE < rhoE_min){
                    rhoE_min = rhoE;
                }
            }

            //solve for components of gradient
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                driver->fluid_body->rhoE_x(e, pos) = sol(pos);
            }

            //limit gradient to ensure monotonicity
            //check max, min at each node
            for (int n=0; n<npe; n++){
                //p(x) = grad(p) * (x - x_0)
                tmp_val = driver->fluid_body->rhoE_x[e].dot(x_n[nodeIDs(e,n)] - x_e[e]);
                //check for both overshoot and undershoot
                if (tmp_val > (rhoE_max - rhoE_0)){
                    //limit gradient
                    for (int pos=0; pos<GRID_DIM; pos++) {
                        driver->fluid_body->rhoE_x[e] *= (rhoE_max - rhoE_0)/tmp_val;
                    }
                } else if (tmp_val < (rhoE_min - rhoE_0)){
                    //limit gradient
                    for (int pos=0; pos<GRID_DIM; pos++) {
                        driver->fluid_body->rhoE_x[e] *= (rhoE_min - rhoE_0)/tmp_val;
                    }
                }
            }

            //transform true energy gradient to effectiv energy gradient
            driver->fluid_body->rhoE_x[e] *= driver->fluid_body->n_e(e);
        }
    }
    return;
}

void FVMGmsh2D::constructPorosityField(Job* job, FiniteVolumeDriver* driver){
    //use porosity field gradients to adjust field gradients
    if (USE_LOCAL_POROSITY_CORRECTION) {
        Eigen::VectorXd gradn_star = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x, u, p, tmp_gradn_star;
        KinematicVector rho_x, rhoE_x;
        KinematicTensor p_x;
        double rho, rhoE, M, c;
        double n_0, n, n_max, n_min, min_dif;
        double tmp_dif, tmp_val;
        double dn_ds, P;
        KinematicVector u_s = KinematicVector(job->JOB_TYPE);

        //loop over elements
        for (int e = 0; e < element_count; e++) {
            //estimate gradn_star (reconstructed porosity gradient)
            if (driver->ORDER == 1) {
                //gradients are initially zero
                gradn_star.setZero();

                tmp_gradn_star = KinematicVector(job->JOB_TYPE);
                tmp_gradn_star.setZero();
            }
            if (driver->ORDER >= 2) {
                if (driver->ORDER > 2) {
                    std::cout << "ERROR: FVMCartesian does not currently implement higher ORDER porosity reconstruction." << std::endl;
                }

                //least squares fit of n to neighbors of element e
                x_0 = getElementCentroid(job, e);
                n_0 = driver->fluid_body->n_e(e);
                n_max = n_0;
                n_min = n_0;

                //create system of equations
                for (int ii = 0; ii < element_neighbors[e].size(); ii++) {
                    //fill in b vector
                    n = driver->fluid_body->n_e(element_neighbors[e][ii]);
                    b_e[e](ii) = n - n_0;

                    //update maximum and minimum velocities
                    if (n > n_max) {
                        n_max = n;
                    } else if (n < n_min) {
                        n_min = n;
                    }
                }

                //solve for components of gradient
                gradn_star = A_inv[e] * b_e[e];

                tmp_gradn_star = KinematicVector(job->JOB_TYPE);
                for (int pos = 0; pos < GRID_DIM; pos++) {
                    tmp_gradn_star[pos] = gradn_star(pos);
                }

                //limit gradient to ensure monotonicity
                //check max, min at each node
                for (int i=0; i<npe; i++){
                    //p(x) = grad(p) * (x - x_0)
                    tmp_val = tmp_gradn_star.dot(x_n[nodeIDs(e,i)] - x_e[e]);
                    //check for both overshoot and undershoot
                    if (tmp_val > (n_max - n_0)){
                        //limit gradient
                        for (int pos=0; pos<GRID_DIM; pos++) {
                            tmp_gradn_star *= (n_max - n_0)/tmp_val;
                        }
                    } else if (tmp_val < (n_min - n_0)){
                        //limit gradient
                        for (int pos=0; pos<GRID_DIM; pos++) {
                            tmp_gradn_star *= (n_min - n_0)/tmp_val;
                        }
                    }
                }
            }

            //update gradient estimates based on difference between n_e_x and gradn_star
            //pg 37-38 nb #7
            //get cell velocity
            p = driver->fluid_body->p(e);
            rho = driver->fluid_body->rho(e);
            rhoE = driver->fluid_body->rhoE(e);
            u = p/rho;
            n = driver->fluid_body->n_e(e);

            //get average solid velocity in cell
            u_s.setZero();
            for (int q=0; q<qpe; q++){
                u_s += v_sq[e*qpe+q] * w_q(e*qpe+q) / getElementVolume(e);
            }

            //convert p, rho, rhoE to solid reference frame
            p -= rho*u_s;
            rhoE += -p.dot(u_s) + 0.5*rho*u_s.dot(u_s);
            u -= u_s;

            //estimate Mach number
            c = driver->fluid_material->getSpeedOfSound(job, driver, rho, p, rhoE, n);
            M = u.norm()/c;

            if (u.norm() > 1e-10 && n > 1e-10 && M < (1.0 - delta)) {
                //estimate dn/ds correction (amount not 'seen' by standard reconstruction operation)
                dn_ds = (driver->fluid_body->n_e_x[e] - tmp_gradn_star).dot(u) / u.norm();

                //calculate gradients of properties in solid frame
                rho_x = rho/n * dn_ds * u/u.norm() * (M*M)/(1.0 - M*M);
                p_x = -p.tensor(u)/(u.norm()*n) * dn_ds;

                if (driver->TYPE == driver->THERMAL) {
                    //get pressure
                    P = driver->fluid_material->getPressure(job, driver, rho, p, rhoE, n);
                    rhoE_x = -(rho*c*c - rhoE - P) / n * dn_ds * u / u.norm() * (M * M) / (1.0 - M * M);
                }

                //convert gradients back to spatial frame and update gradients
                driver->fluid_body->rho_x[e] += rho_x;
                driver->fluid_body->p_x[e] += p_x + u_s.tensor(rho_x);

                if (driver->TYPE == driver->THERMAL) {
                    driver->fluid_body->rhoE_x[e] += rhoE_x + p_x.transpose()*u_s + 0.5*u_s.dot(u_s)*rho_x;
                }
            } else {
                //flow too slow, too low porosity, or trans/supersonic
                //do nothing
            }
        }
    } else {
        //do nothing
    }

    return;
}
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

    //if integers are passed, they define boundary conditions
    std::vector<int> tmp_bc_tags = std::vector<int>(0);
    std::vector<KinematicVector> tmp_bc_values = std::vector<KinematicVector>(0);
    std::vector<double> tmp_bc_density = std::vector<double>(0);                    //for supersonic inlet
    int supersonic_props = 0;
    if (int_props.size() >= 1){
        tmp_bc_tags.resize(int_props.size());
        tmp_bc_values.resize(int_props.size());
        tmp_bc_density.resize(int_props.size());

        if (fp64_props.size() < int_props.size()*GRID_DIM){
            std::cerr << "ERROR: FVMGmsh2D passed mismatched integer (" << int_props.size() << ") and double (" << fp64_props.size() << ") properties!" << std::endl;
            exit(0);
        }

        for (int i = 0; i<int_props.size(); i++){
            tmp_bc_tags[i] = int_props[i];

            //check for supersonic inlet
            if (tmp_bc_tags[i] == SUPERSONIC_INLET){
                //need to have another floating point inlet property
                supersonic_props++;
                if (fp64_props.size() < int_props.size()*GRID_DIM + supersonic_props){
                    std::cerr << "ERROR: FVMGmsh2D passed mismatched integer ("
                              << int_props.size() << ") and double ("
                              << fp64_props.size() << ") properties with "
                              << supersonic_props << " additional properties required!" << std::endl;
                    exit(0);
                }
                tmp_bc_density[i] = fp64_props[2*i + supersonic_props - 1];
            } else {
                tmp_bc_density[i] = -1; //this should never actually be used
            }

            tmp_bc_values[i][0] = fp64_props[2*i + supersonic_props];
            tmp_bc_values[i][1] = fp64_props[2*i + 1 + supersonic_props];

            //check that bc tag is in coded list
            if (tmp_bc_tags[i] != DIRICHLET
                && tmp_bc_tags[i] != NEUMANN
                && tmp_bc_tags[i] != PERIODIC
                && tmp_bc_tags[i] != NEUMANN_DAMPING
                && tmp_bc_tags[i] != SUPERSONIC_INLET){
                std::cerr << "ERROR: Boundary tag " << tmp_bc_tags[i] << " not defined for FVMGmsh2D grid object! Exiting." << std::endl;
                exit(0);
            }
        }
    }


    if (job->JOB_TYPE != 2){
        std::cerr << "WARNING: Cannot use FVMGmsh2D with JOB_TYPE of " << job->JOB_TYPE << "!" << std::endl;
    }

    std::cout << "FiniteVolumeGrid properties: (filename = " << filename << ")" << std::endl;
    for (int i=0; i<tmp_bc_tags.size(); i++){
        std::cout << "    " << i << " : " << tmp_bc_tags[i] << " : " << tmp_bc_density[i] << " - <" << tmp_bc_values[i][0] << ", " << tmp_bc_values[i][1] << ">" << std::endl;
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
    std::vector<std::vector<int>> tmp_node_to_element_map = std::vector<std::vector<int>>(node_count);
    for (int e=0; e<element_count; e++){
        //fill in list of elements associated with each node
        n0 = nodeIDs(e,0);
        n1 = nodeIDs(e,1);
        n2 = nodeIDs(e,2);

        tmp_node_to_element_map[n0].push_back(e);
        tmp_node_to_element_map[n1].push_back(e);
        tmp_node_to_element_map[n2].push_back(e);
    }
    int e_tmp;
    bool in_list;
    for (int e=0; e<element_count; e++){
        //loop over elements
        for (int j=0; j<npe; j++) {
            //for each node
            n0 = nodeIDs(e,j);
            //check if elements associated with that node are already in the neighbor list
            for (int i = 0; i < tmp_node_to_element_map[n0].size(); i++) {
                e_tmp = tmp_node_to_element_map[n0][i];
                in_list = false;
                for (int ii = 0; i < element_neighbors[e].size(); i++) {
                    if (element_neighbors[e][ii] == e_tmp){
                        in_list = true;
                        break;
                    }
                }
                if (!in_list){
                    //if element is not in neighbor list, add it
                    element_neighbors[e].push_back(e_tmp);
                }
            }
        }
    }

    //define boundary conditions
    bc_tags = Eigen::VectorXi(face_count);
    bc_tags.setConstant(-1);                                        //by default, interior faces are set to -1
    bc_values = KinematicVectorArray(face_count, job->JOB_TYPE);
    bc_values.setZero();                                            //by default, boundary conditions are zero
    bc_density = Eigen::VectorXd(face_count);
    bc_density.setZero();

    //loop over stored boundary faces
    int which_boundary;
    for (int i=0; i<bounding_faces.size(); i++){
        which_boundary = bounding_faces[i][0];
        n0 = bounding_faces[i][1];
        n1 = bounding_faces[i][2];
        //check faces associated with n0 for n1 from node to face map create earlier
        for (int i=0; i<tmp_node_to_face_map[n0].size(); i++){
            if (n1 == tmp_node_to_face_map[n0][i][0]){
                //found face
                f_tmp = tmp_node_to_face_map[n0][i][1];
                break;
            }
        }
        if (which_boundary >= tmp_bc_tags.size()){
            //if no bc defined by input file, assume dirichlet
            bc_tags[f_tmp] = DIRICHLET;
        } else {
            bc_tags[f_tmp] = tmp_bc_tags[which_boundary];
            bc_values[f_tmp] = tmp_bc_values[which_boundary];
            bc_density[f_tmp] = tmp_bc_density[which_boundary];
        }
    }

    //lastly, need to form least squares system of equations
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
        for (int j = 0; j < element_faces[e].size(); j++) {
            //only add BCs for dirichlet conditions
            int f = element_faces[e][j];
            if (bc_tags[f] == DIRICHLET || bc_tags[f] == SUPERSONIC_INLET) {
                length_of_A++; //add dirichlet conditions to A
            }
        }

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

        //add dirichlet boundary conditions
        int i = element_neighbors[e].size();
        for (int j = 0; j < element_faces[e].size(); j++) {
            //only add BCs for dirichlet conditions
            int f = element_faces[e][j];
            if (bc_tags[f] == DIRICHLET || bc_tags[f] == SUPERSONIC_INLET) {
                x = getFaceCentroid(job, f);
                for (int pos = 0; pos < GRID_DIM; pos++) {
                    tmp_dif = x[pos] - x_0[pos];
                    A_e[e](i, pos) = x[pos] - x_0[pos];
                }
                //increment counter
                i++;
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
            a_tmp = AtA.householderQr().solve(e_vec);
            for (int jj=0; jj<GRID_DIM; jj++){
                AtA_inv(jj,ii) = a_tmp(jj);
            }
        }
        A_inv[e] = AtA_inv*A_e[e].transpose();
    }

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
//helper functions
double FVMGmsh2D::getElementVolume(int e){
    return v_e(e);
}

int FVMGmsh2D::getElementTag(int e){
    //elements are not tagged
    return -1;
}

double FVMGmsh2D::getFaceArea(int f){
    return face_areas(f);
}

int FVMGmsh2D::getFaceTag(int f){
    return bc_tags(f);
}

KinematicVector FVMGmsh2D::getFaceNormal(Job* job, int f){
    return face_normals[f];
}

std::vector<int> FVMGmsh2D::getElementFaces(int e){
    return element_faces[e];
}

std::array<int,2> FVMGmsh2D::getOrientedElementsByFace(int f){
    return face_elements[f];
}

std::vector<int> FVMGmsh2D::getElementNeighbors(int e){
    return element_neighbors[e];
}

KinematicVector FVMGmsh2D::getElementCentroid(Job* job, int e){
    return x_e[e];
}

KinematicVector FVMGmsh2D::getFaceCentroid(Job* job, int f){
    return x_f[f];
}

/*----------------------------------------------------------------------------*/
//generate mapping matrix between MPM grid and FVM grid
void FVMGmsh2D::generateMappings(Job* job, FiniteVolumeDriver* driver){
    //for now do nothing
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
            p_0 = driver->fluid_body->p[e];
            rho_0 = driver->fluid_body->rho(e);
            p_max = p_0; p_min = p_0;
            for (int mom_index = 0; mom_index<GRID_DIM; mom_index++){
                //create system of equations
                for (int ii=0; ii<element_neighbors[e].size(); ii++){
                    //ill in b vector
                    p = driver->fluid_body->p[element_neighbors[e][ii]];
                    b_e[e](ii) = p[mom_index] - p_0[mom_index];

                    //update maximum and minimum velocities
                    if (p[mom_index] > p_max[mom_index]){
                        p_max[mom_index] = p[mom_index];
                    } else if (p[mom_index] < p_min[mom_index]){
                        p_min[mom_index] = p[mom_index];
                    }
                }

                //add boundary conditions where applicable
                int i=element_neighbors[e].size();
                for (int j = 0; j<element_faces[e].size(); j++){
                    //only add BCs for dirichlet conditions
                    int f = element_faces[e][j];
                    if (bc_tags[f] == FiniteVolumeGrid::DIRICHLET) {
                        //fill in b vector
                        p = rho_0 * bc_values[f];
                        b_e[e](i) = p[mom_index] - p_0[mom_index];

                        //update maximum and minimum velocities
                        if (p[mom_index] > p_max[mom_index]) {
                            p_max[mom_index] = p[mom_index];
                        } else if (p[mom_index] < p_min[mom_index]) {
                            p_min[mom_index] = p[mom_index];
                        }
                        //increment i
                        i++;
                    } else if (bc_tags[f] == SUPERSONIC_INLET) {
                        //fill in b vector
                        p = bc_density[f] * bc_values[f];
                        b_e[e](i) = p[mom_index] - p_0[mom_index];

                        //update maximum and minimum velocities
                        if (p[mom_index] > p_max[mom_index]) {
                            p_max[mom_index] = p[mom_index];
                        } else if (p[mom_index] < p_min[mom_index]) {
                            p_min[mom_index] = p[mom_index];
                        }
                        //increment i
                        i++;
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
            rho_0 = driver->fluid_body->rho(e);
            rho_max = rho_0;
            rho_min = rho_0;

            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                //fill in b vector
                rho = driver->fluid_body->rho(element_neighbors[e][ii]);
                b_e[e](ii) = rho - rho_0;

                //update maximum and minimum velocities
                if (rho > rho_max){
                    rho_max = rho;
                } else if (rho < rho_min){
                    rho_min = rho;
                }
            }

            //add boundary conditions where applicable (total number of eq'ns still less than rows of A)
            int i=element_neighbors[e].size();
            for (int j = 0; j<element_faces[e].size(); j++){
                //only add BCs for dirichlet conditions
                int f = element_faces[e][j];
                if (bc_tags[f] == FiniteVolumeGrid::DIRICHLET) {
                    b_e[e](i) = 0;

                    //increment i
                    i++;
                } else if (bc_tags[f] == SUPERSONIC_INLET) {
                    b_e[e](i) = bc_density[f] - rho_0;

                    //increment i
                    i++;
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
        }
    }
    return;
}

KinematicTensorArray FVMGmsh2D::getVelocityGradients(Job* job, FiniteVolumeDriver* driver){
    //reconstruct velocity field
    KinematicTensorArray u_x = KinematicTensorArray(element_count, job->JOB_TYPE);
    KinematicVector u = KinematicVector(job->JOB_TYPE);
    KinematicVector p = KinematicVector(job->JOB_TYPE);
    double rho, rho_0, tmp_dif;

    //initialize least squares fitting variables
    Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
    KinematicVector x_0, x;
    KinematicVector u_0;

    for (int e=0; e<element_count; e++){
        p = driver->fluid_body->p[e];
        rho_0 = driver->fluid_body->rho(e);
        u_0 = p/rho_0;

        //least squares fit of u_x to neighbors of element e
        x_0 = getElementCentroid(job, e);
        for (int dir = 0; dir<GRID_DIM; dir++){
            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                rho = driver->fluid_body->rho(element_neighbors[e][ii]);
                p = driver->fluid_body->p[element_neighbors[e][ii]];
                b_e[e](ii) = p[dir]/rho - u_0[dir];
            }

            //add boundary conditions where applicable (total number of eq'ns still less than rows of A)
            int i=element_neighbors[e].size();
            for (int j = 0; j<element_faces[e].size(); j++){
                //only add BCs for dirichlet conditions
                int f = element_faces[e][j];
                if (bc_tags[f] == DIRICHLET || bc_tags[f] == SUPERSONIC_INLET) {
                    u = bc_values[f];
                    b_e[e](i) = u[dir] - u_0[dir];

                    //increment i
                    i++;
                }
            }

            //solve for u_pos component of gradient
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                u_x(e, dir, pos) = sol(pos);
            }
        }

    }
    return u_x;
}

/*----------------------------------------------------------------------------*/
//functions to compute element-wise fluxes of field variables using reconstructed velocity field
Eigen::VectorXd FVMGmsh2D::calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, Eigen::VectorXd& phi){
    if (phi.rows() != element_count){
        std::cerr << "ERROR! Length of vector passed to calculateElementFluxIntegrals is wrong size! ";
        std::cerr << phi.rows() << " != " << element_count << std::endl;
        exit(0);
    }

    //flux rates
    Eigen::VectorXd result = Eigen::VectorXd(element_count);
    result.setZero();

    double v_plus, v_minus; //value
    KinematicVector u_plus, u_minus, normal, u_bar; //velocity
    int e_plus, e_minus;    //elements
    double rho_plus, rho_minus;
    double flux, area;
    double a_1, a_2, a_6, c, rho_bar, phi_bar, theta_bar;
    double lambda_1, lambda_2;

    //number of quadrature points per face
    int num_quad = 2;

    //vector of quad point positions relative to face center
    std::array<double,2> quad_points = {-1.0/std::sqrt(3), 1.0/std::sqrt(3)};

    //different flux calculation for different ORDER approximations
    if (driver->ORDER == 1){
        //rho, u, value constant within element
        //loop over faces
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);
        KinematicVector x;

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (bc_tags[f] == -1){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];
                v_minus = phi(e_minus);
                v_plus  = phi(e_plus);
                rho_minus = driver->fluid_body->rho(e_minus);
                rho_plus = driver->fluid_body->rho(e_plus);
                u_minus = driver->fluid_body->p[e_minus]/rho_minus; //u = p/rho
                u_plus  = driver->fluid_body->p[e_plus]/rho_plus;

                //approximate Roe advective rate
                u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                rho_bar = std::sqrt(rho_plus*rho_minus);
                phi_bar = (v_plus*std::sqrt(rho_minus) + v_minus*std::sqrt(rho_plus))/(std::sqrt(rho_minus) + std::sqrt(rho_plus));
                theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                c = driver->fluid_material->getSpeedOfSound(job, driver, x, rho_bar, theta_bar);

                //roe eigenvalues
                lambda_1 = std::abs(u_bar.dot(normal) - c);
                lambda_2 = std::abs(u_bar.dot(normal) + c);
                if (lambda_1 < delta*c){
                    lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                }
                if (lambda_2 < delta*c){
                    lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                }

                //calculate Roe eigenvector coefficients
                a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                a_2 = (rho_plus - rho_minus) - a_1;
                a_6 = v_plus - v_minus - (rho_plus - rho_minus)*phi_bar/rho_bar;

                //flux in n direction
                flux = area*0.5*(v_plus*u_plus.dot(normal) + v_minus*u_minus.dot(normal)
                                 - a_1*lambda_1*phi_bar/rho_bar
                                 - a_2*lambda_2*phi_bar/rho_bar
                                 - a_6*std::abs(u_bar.dot(normal)));

                //add flux to element integrals
                result(e_minus) -= flux;
                result(e_plus) += flux;
            } else if (bc_tags[f] == DIRICHLET || bc_tags[f] == SUPERSONIC_INLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and assigned velocity
                    v_minus = phi(e_minus);
                    flux = v_minus*area*bc_values[f].dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and assigned velocity
                    v_plus = phi(e_plus);
                    flux = v_plus*area*bc_values[f].dot(normal);
                    result(e_plus) += flux;
                }
            } else if (bc_tags[f] == NEUMANN || bc_tags[f] == NEUMANN_DAMPING){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and cell velocity
                    v_minus = phi(e_minus);
                    u_minus = driver->fluid_body->p[e_minus]/driver->fluid_body->rho(e_minus);
                    flux = v_minus*area*u_minus.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and cell velocity
                    v_plus = phi(e_plus);
                    u_plus = driver->fluid_body->p[e_plus]/driver->fluid_body->rho(e_plus);
                    flux = v_plus*area*u_plus.dot(normal);
                    result(e_plus) += flux;
                }
            }
        }
    } else if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMGmsh2D does not currently implement higher ORDER flux reconstruction."
                      << std::endl;
        }

        //first reconstruct input field
        KinematicVectorArray phi_x = KinematicVectorArray(element_count, job->JOB_TYPE);

        //least squares fit
        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x, x_face, x_quad;
        x_quad = KinematicVector(job->JOB_TYPE);
        double value_0, value, value_max, value_min, min_dif;
        double tmp_dif, tmp_val;
        for (int e = 0; e<element_count; e++){
            //least squares fit of phi to neighbors of element e
            x_0 = getElementCentroid(job, e);
            value_0 = phi(e);
            min_dif = 0;

            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                value = phi(element_neighbors[e][ii]);
                b_e[e](ii) = value - value_0;

                //update maximum and minimum velocities
                if (value > value_max){
                    value_max = value;
                } else if (value < value_min){
                    value_min = value;
                }
            }

            //zero remainder of system of eq'ns
            for (int ii=element_neighbors[e].size(); ii<b_e[e].rows(); ii++){
                for (int pos=0; pos<GRID_DIM; pos++){
                    b_e[e](ii) = 0;
                }
            }

            //solve for u_pos component of gradient
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                phi_x(e, pos) = sol(pos);
            }

            //limit gradient to ensure monotonicity
            //check max, min at each node
            for (int n=0; n<npe; n++){
                //p(x) = grad(p) * (x - x_0)
                tmp_val = phi_x[e].dot(x_n[nodeIDs(e,n)] - x_e[e]);
                //check for both overshoot and undershoot
                if (tmp_val > (value_max - value_0)){
                    //limit gradient
                    for (int pos=0; pos<GRID_DIM; pos++) {
                        phi_x[e] *= (value_max - value_0)/tmp_val;
                    }
                } else if (tmp_val < (value_min - value_0)){
                    //limit gradient
                    for (int pos=0; pos<GRID_DIM; pos++) {
                        phi_x[e] *= (value_min - value_0)/tmp_val;
                    }
                }
            }
        }

        //loop over faces and use quadrature to reconstruct flux integral
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x_face = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (bc_tags[f] == -1){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    x = x_f[f] - x_e[e_minus];
                    x[0] += -quad_points[q]*normal[1]*area;
                    x[1] += quad_points[q]*normal[0]*area;

                    //calculate A properties
                    v_minus = phi(e_minus) + phi_x(e_minus).dot(x);
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    u_minus = (driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x)/rho_minus;

                    //relative position to centroid of B
                    x = x_f[f] - x_e[e_plus];
                    x[0] += -quad_points[q]*normal[1]*area;
                    x[1] += quad_points[q]*normal[0]*area;

                    //calculate B properties
                    v_plus = phi(e_plus) + phi_x(e_plus).dot(x);
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x)/rho_plus;

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = std::sqrt(rho_plus*rho_minus);
                    phi_bar = (v_plus*std::sqrt(rho_minus) + v_minus*std::sqrt(rho_plus))/(std::sqrt(rho_minus) + std::sqrt(rho_plus));
                    theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                    c = driver->fluid_material->getSpeedOfSound(job, driver, x_quad, rho_bar, theta_bar);

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    if (lambda_1 < delta*c){
                        lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                    }
                    if (lambda_2 < delta*c){
                        lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                    }

                    //calculate Roe eigenvector coefficients
                    a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                    a_2 = (rho_plus - rho_minus) - a_1;
                    a_6 = v_plus - v_minus - (rho_plus - rho_minus)*phi_bar/rho_bar;

                    //flux in n direction
                    flux = area/num_quad * 0.5 * (v_plus*u_plus.dot(normal) + v_minus*u_minus.dot(normal)
                                                  - a_1*lambda_1*phi_bar/rho_bar
                                                  - a_2*lambda_2*phi_bar/rho_bar
                                                  - a_6*std::abs(u_bar.dot(normal)));

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else if (bc_tags[f] == DIRICHLET || bc_tags[f] == SUPERSONIC_INLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {

                    if (e_minus > -1) {
                        //relative position to centroid of A
                        x = x_f[f] - x_e[e_minus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;

                        //calculate A properties
                        v_minus = phi(e_minus) + phi_x(e_minus).dot(x);
                        flux = v_minus*area/num_quad*bc_values[f].dot(normal);
                        result(e_minus) -= flux;
                    }

                    if (e_plus > -1) {
                        //relative position to centroid of B
                        x = x_f[f] - x_e[e_plus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;

                        //calculate B properties
                        v_plus = phi(e_plus) + phi_x(e_plus).dot(x);
                        flux = v_plus*area/num_quad*bc_values[f].dot(normal);
                        result(e_plus) += flux;
                    }
                }
            } else if (bc_tags[f] == NEUMANN || bc_tags[f] == NEUMANN_DAMPING){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {

                    if (e_minus > -1) {
                        //relative position to centroid of A
                        x = x_f[f] - x_e[e_minus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;

                        //calculate A properties
                        v_minus = phi(e_minus) + phi_x(e_minus).dot(x);
                        rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        u_minus = (driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x)/rho_minus;
                        flux = v_minus*area/num_quad*u_minus.dot(normal);
                        result(e_minus) -= flux;
                    }

                    if (e_plus > -1) {
                        //relative position to centroid of B
                        x = x_f[f] - x_e[e_plus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;

                        //calculate B properties
                        v_plus = phi(e_plus) + phi_x(e_plus).dot(x);
                        rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x)/rho_plus;
                        flux = v_plus*area/num_quad*u_plus.dot(normal);
                        result(e_plus) += flux;
                    }
                }
            }
        }
    }
    return result;
}

//functions to compute element mass flux
Eigen::VectorXd FVMGmsh2D::calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver){
    //flux rates
    Eigen::VectorXd result = Eigen::VectorXd(element_count);
    result.setZero();

    KinematicVector u_plus, u_minus, normal, u_bar; //velocity
    int e_plus, e_minus;    //elements
    double rho_plus, rho_minus;
    double flux, area;
    double lambda_1, lambda_2, a_1, a_2, c, rho_bar, theta_bar;

    //number of quadrature points per face
    int num_quad = 2;

    //vector of quad point positions relative to face center
    std::array<double,2> quad_points = {-1.0/std::sqrt(3), 1.0/std::sqrt(3)};

    //different flux calculation for different ORDER approximations
    if (driver->ORDER == 1){
        //rho, u, value constant within element
        //loop over faces
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);
        KinematicVector x;

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (bc_tags[f] == -1){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];
                rho_minus = driver->fluid_body->rho(e_minus);
                rho_plus = driver->fluid_body->rho(e_plus);
                u_minus = driver->fluid_body->p[e_minus]/rho_minus; //u = p/rho
                u_plus  = driver->fluid_body->p[e_plus]/rho_plus;

                //approximate Roe advective rate
                u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                rho_bar = std::sqrt(rho_plus*rho_minus);
                theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                c = driver->fluid_material->getSpeedOfSound(job, driver, x, rho_bar, theta_bar);

                //roe eigenvalues
                lambda_1 = std::abs(u_bar.dot(normal) - c);
                lambda_2 = std::abs(u_bar.dot(normal) + c);
                if (lambda_1 < delta*c){
                    lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                }
                if (lambda_2 < delta*c){
                    lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                }

                //calculate Roe eigenvector coefficients
                a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                a_2 = (rho_plus - rho_minus) - a_1;

                //flux in n direction
                flux = area * 0.5 * (rho_plus*u_plus.dot(normal) + rho_minus*u_minus.dot(normal)
                                     - a_1*lambda_1 - a_2*lambda_2);

                //add flux to element integrals
                result(e_minus) -= flux;
                result(e_plus) += flux;
            } else if (bc_tags[f] == DIRICHLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and assigned velocity
                    rho_minus = driver->fluid_body->rho(e_minus);
                    flux = rho_minus*area*bc_values[f].dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and assigned velocity
                    rho_plus = driver->fluid_body->rho(e_plus);
                    flux = rho_plus*area*bc_values[f].dot(normal);
                    result(e_plus) += flux;
                }
            } else if (bc_tags[f] == SUPERSONIC_INLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and assigned velocity
                    rho_minus = bc_density[f];
                    flux = rho_minus*area*bc_values[f].dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and assigned velocity
                    rho_plus = bc_density[f];
                    flux = rho_plus*area*bc_values[f].dot(normal);
                    result(e_plus) += flux;
                }
            } else if (bc_tags[f] == NEUMANN || bc_tags[f] == NEUMANN_DAMPING){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and cell velocity
                    rho_minus = driver->fluid_body->rho(e_minus);
                    u_minus = driver->fluid_body->p[e_minus]/driver->fluid_body->rho(e_minus);
                    flux = rho_minus*area*u_minus.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and cell velocity
                    rho_plus = driver->fluid_body->rho(e_plus);
                    u_plus = driver->fluid_body->p[e_plus]/driver->fluid_body->rho(e_plus);
                    flux = rho_plus*area*u_plus.dot(normal);
                    result(e_plus) += flux;
                }
            }
        }
    } else if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMGmsh2D does not currently implement higher ORDER flux reconstruction."
                      << std::endl;
        }

        KinematicVector x = KinematicVector(job->JOB_TYPE);
        KinematicVector x_face = KinematicVector(job->JOB_TYPE);
        KinematicVector x_quad = KinematicVector(job->JOB_TYPE);

        //loop over faces and use quadrature to reconstruct flux integral
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x_face = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (bc_tags[f] == -1){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    x = x_f[f] - x_e[e_minus];
                    x[0] += -quad_points[q]*normal[1]*area;
                    x[1] += quad_points[q]*normal[0]*area;

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    u_minus = (driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x)/rho_minus;

                    //relative position to centroid of B
                    x = x_f[f] - x_e[e_plus];
                    x[0] += -quad_points[q]*normal[1]*area;
                    x[1] += quad_points[q]*normal[0]*area;

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x)/rho_plus;

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = std::sqrt(rho_plus*rho_minus);
                    theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                    c = driver->fluid_material->getSpeedOfSound(job, driver, x_quad, rho_bar, theta_bar);

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    if (lambda_1 < delta*c){
                        lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                    }
                    if (lambda_2 < delta*c){
                        lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                    }

                    //calculate Roe eigenvector coefficients
                    a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                    a_2 = (rho_plus - rho_minus) - a_1;

                    //flux in n direction
                    flux = area/num_quad * 0.5 * (rho_plus*u_plus.dot(normal) + rho_minus*u_minus.dot(normal)
                                                  - a_1*lambda_1 - a_2*lambda_2);

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else if (bc_tags[f] == DIRICHLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {

                    if (e_minus > -1) {
                        //relative position to centroid of A
                        x = x_f[f] - x_e[e_minus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        //calculate A properties
                        rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        flux = rho_minus*area/num_quad*bc_values[f].dot(normal);
                        result(e_minus) -= flux;
                    }

                    if (e_plus > -1) {
                        //relative position to centroid of B
                        x = x_f[f] - x_e[e_plus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        //calculate B properties
                        rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        flux = rho_plus*area/num_quad*bc_values[f].dot(normal);
                        result(e_plus) += flux;
                    }
                }
            } else if (bc_tags[f] == SUPERSONIC_INLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {

                    if (e_minus > -1) {
                        //relative position to centroid of A
                        x = x_f[f] - x_e[e_minus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        //calculate A properties
                        rho_minus = bc_density[f];
                        flux = rho_minus*area/num_quad*bc_values[f].dot(normal);
                        result(e_minus) -= flux;
                    }

                    if (e_plus > -1) {
                        //relative position to centroid of B
                        x = x_f[f] - x_e[e_plus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        //calculate B properties
                        rho_plus = bc_density[f];
                        flux = rho_plus*area/num_quad*bc_values[f].dot(normal);
                        result(e_plus) += flux;
                    }
                }
            } else if (bc_tags[f] == NEUMANN || bc_tags[f] == NEUMANN_DAMPING){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {

                    if (e_minus > -1) {
                        //relative position to centroid of A
                        x = x_f[f] - x_e[e_minus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;

                        //calculate A properties
                        rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        u_minus = (driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x)/rho_minus;
                        flux = rho_minus*area/num_quad*u_minus.dot(normal);
                        result(e_minus) -= flux;
                    }


                    if (e_plus > -1) {
                        //relative position to centroid of B
                        x = x_f[f] - x_e[e_plus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;

                        //calculate B properties
                        rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x)/rho_plus;
                        flux = rho_plus*area/num_quad*u_plus.dot(normal);
                        result(e_plus) += flux;
                    }
                }
            }
        }
    }
    return result;
}

//functions to compute element momentum fluxes
KinematicVectorArray FVMGmsh2D::calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver){

    //flux rates
    KinematicVectorArray result = KinematicVectorArray(element_count, job->JOB_TYPE);
    result.setZero();

    //intermediate variables
    KinematicVector u_plus, u_minus, normal, u_bar, flux; //velocity
    KinematicVector p_plus, p_minus;
    int e_plus, e_minus;    //elements
    double rho_plus, rho_minus;
    double area;
    double rho_bar, c, a_1, a_2, lambda_1, lambda_2;
    KinematicVector x_0 = KinematicVector(job->JOB_TYPE);
    KinematicVector x_face = KinematicVector(job->JOB_TYPE);
    KinematicVector x_quad = KinematicVector(job->JOB_TYPE);
    KinematicVector a_3 = KinematicVector(job->JOB_TYPE);

    //traction calculatons
    KinematicTensorArray L = getVelocityGradients(job, driver);
    KinematicTensor L_tmp = KinematicTensor(job->JOB_TYPE);
    MaterialTensor tau_plus, tau_minus;
    double P_plus, P_minus, theta_bar;

    //number of quadrature points per face
    int num_quad = 2;

    //vector of quad point positions relative to face center
    std::array<double,2> quad_points = {-1.0/std::sqrt(3), 1.0/std::sqrt(3)};

    //different flux calculation for different ORDER approximations
    if (driver->ORDER == 1){
        //rho, u, value constant within element
        //loop over faces
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x_face = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (bc_tags[f] == -1){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];
                rho_minus = driver->fluid_body->rho(e_minus);
                rho_plus = driver->fluid_body->rho(e_plus);
                p_minus = driver->fluid_body->p[e_minus];
                p_plus = driver->fluid_body->p[e_plus];
                u_minus = p_minus/rho_minus; //u = p/rho
                u_plus  = p_plus/rho_plus;

                //approximate Roe advective rate
                u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                rho_bar = std::sqrt(rho_plus*rho_minus);
                theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                c = driver->fluid_material->getSpeedOfSound(job, driver, x_face, rho_bar, theta_bar);

                //roe eigenvalues
                lambda_1 = std::abs(u_bar.dot(normal) - c);
                lambda_2 = std::abs(u_bar.dot(normal) + c);
                if (lambda_1 < delta*c){
                    lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                }
                if (lambda_2 < delta*c){
                    lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                }

                //calculate Roe eigenvector coefficients
                a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                a_2 = (rho_plus - rho_minus) - a_1;
                a_3 = p_plus - p_minus - (rho_plus-rho_minus)*u_bar;
                a_3 = a_3 - a_3.dot(normal)*normal; //remove normal component of a_3 vector

                //flux in n direction
                flux = area*0.5*(p_plus*u_plus.dot(normal) + p_minus*u_minus.dot(normal)
                                 - a_1*lambda_1*(u_bar - c*normal)
                                 - a_2*lambda_2*(u_bar + c*normal)
                                 - a_3*std::abs(u_bar.dot(normal)));

                //add tractions to momentum flux
                //for now use simple reconstruction of theta
                tau_minus = driver->fluid_material->getShearStress(job, driver, x_face, L[e_minus], rho_bar, theta_bar);
                tau_plus  = driver->fluid_material->getShearStress(job, driver, x_face, L[e_plus], rho_bar, theta_bar);
                P_plus = c*c*(rho_plus - rho_bar) + driver->fluid_material->getPressure(job, driver, x_face, rho_bar, theta_bar);
                P_minus = c*c*(rho_minus - rho_plus) + P_plus;

                flux += area * 0.5 * ((P_plus + P_minus)*normal - KinematicVector((tau_plus + tau_minus)*normal, job->JOB_TYPE));

                //add flux to element integrals
                result(e_minus) -= flux;
                result(e_plus) += flux;
            } else if (bc_tags[f] == DIRICHLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and assigned velocity
                    p_minus = driver->fluid_body->p[e_minus];

                    //tractions
                    rho_bar = driver->fluid_body->rho(e_minus);
                    theta_bar = driver->fluid_body->theta(e_minus);

                    //estimate L
                    x_0 = getElementCentroid(job, e_minus);
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_values[f][ii] - p_minus[ii]/rho_bar)/(x_face - x_0).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, x_face, L_tmp, rho_bar, theta_bar);
                    P_minus = driver->fluid_material->getPressure(job, driver, x_face, rho_bar, theta_bar);

                    flux = p_minus*area*bc_values[f].dot(normal)
                           + area*P_minus*normal
                           - area*KinematicVector(tau_minus*normal, job->JOB_TYPE);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and assigned velocity
                    p_plus = driver->fluid_body->p[e_plus];

                    //tractions
                    rho_bar = driver->fluid_body->rho(e_plus);
                    theta_bar = driver->fluid_body->theta(e_plus);

                    //estimate L
                    x_0 = getElementCentroid(job, e_plus);
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_values[f][ii] - p_plus[ii]/rho_bar)/(x_face - x_0).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, x_face, L_tmp, rho_bar, theta_bar);
                    P_plus = driver->fluid_material->getPressure(job, driver, x_face, rho_bar, theta_bar);

                    flux = p_plus*area*bc_values[f].dot(normal)
                           + area*P_plus*normal
                           - area*KinematicVector(tau_plus*normal, job->JOB_TYPE);

                    result(e_plus) += flux;
                }
            } else if (bc_tags[f] == SUPERSONIC_INLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and assigned velocity
                    p_minus = driver->fluid_body->p[e_minus];

                    //tractions
                    rho_bar = driver->fluid_body->rho(e_minus);
                    theta_bar = driver->fluid_body->theta(e_minus);

                    //estimate L
                    x_0 = getElementCentroid(job, e_minus);
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_values[f][ii] - p_minus[ii]/rho_bar)/(x_face - x_0).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, x_face, L_tmp, bc_density[f], theta_bar);
                    P_minus = driver->fluid_material->getPressure(job, driver, x_face, bc_density[f], theta_bar);

                    flux = bc_density[f]*bc_values[f]*area*bc_values[f].dot(normal)
                           + area*P_minus*normal
                           - area*KinematicVector(tau_minus*normal, job->JOB_TYPE);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and assigned velocity
                    p_plus = driver->fluid_body->p[e_plus];

                    //tractions
                    rho_bar = driver->fluid_body->rho(e_plus);
                    theta_bar = driver->fluid_body->theta(e_plus);

                    //estimate L
                    x_0 = getElementCentroid(job, e_plus);
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_values[f][ii] - p_plus[ii]/rho_bar)/(x_face - x_0).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, x_face, L_tmp, bc_density[f], theta_bar);
                    P_plus = driver->fluid_material->getPressure(job, driver, x_face, bc_density[f], theta_bar);

                    flux = bc_density[f]*bc_values[f]*area*bc_values[f].dot(normal)
                           + area*P_plus*normal
                           - area*KinematicVector(tau_plus*normal, job->JOB_TYPE);

                    result(e_plus) += flux;
                }
            } else if (bc_tags[f] == NEUMANN){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and cell velocity
                    rho_minus = driver->fluid_body->rho(e_minus);
                    p_minus = driver->fluid_body->p[e_minus];
                    u_minus = p_minus/rho_minus;
                    flux = p_minus*area*u_minus.dot(normal);

                    //add traction directly to flux integral
                    result(e_minus) -= flux;
                    result(e_minus) += area*bc_values[f];
                }

                if (e_plus > -1){
                    //flux defined by cell value and cell velocity
                    rho_plus = driver->fluid_body->rho(e_plus);
                    p_plus = driver->fluid_body->p[e_plus];
                    u_plus = p_plus/rho_plus;
                    flux = p_plus*area*u_plus.dot(normal);

                    //add traction directly to flux integral
                    result(e_plus) += flux;
                    result(e_minus) += area*bc_values[f];
                }
            } else if (bc_tags[f] == NEUMANN_DAMPING){
                //face has prescribed traction for weighted average
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and cell velocity
                    rho_minus = driver->fluid_body->rho(e_minus);
                    p_minus = driver->fluid_body->p[e_minus];
                    u_minus = p_minus/rho_minus;
                    flux = p_minus*area*u_minus.dot(normal);

                    theta_bar = driver->fluid_body->theta(e_minus);
                    P_minus = driver->fluid_material->getPressure(job, driver, x_face, rho_minus, theta_bar);

                    //add traction directly to flux integral
                    result(e_minus) -= flux + lambda*area*P_minus*normal;
                    result(e_minus) += (1.0-lambda)*area*bc_values[f];
                }

                if (e_plus > -1){
                    //flux defined by cell value and cell velocity
                    rho_plus = driver->fluid_body->rho(e_plus);
                    p_plus = driver->fluid_body->p[e_plus];
                    u_plus = p_plus/rho_plus;
                    flux = p_plus*area*u_plus.dot(normal);

                    theta_bar = driver->fluid_body->theta(e_plus);
                    P_plus = driver->fluid_material->getPressure(job, driver, x_face, rho_plus, theta_bar);

                    //add traction directly to flux integral
                    result(e_plus) += flux + lambda*area*P_plus*normal;
                    result(e_plus) += (1.0 - lambda)*area*bc_values[f];
                }
            }
        }
    } else if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMCartesian does not currently implement higher ORDER flux reconstruction."
                      << std::endl;
        }

        KinematicVector x = KinematicVector(job->JOB_TYPE);

        //loop over faces and use quadrature to reconstruct flux integral
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x_face = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (bc_tags[f] == -1){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    x = x_f[f] - x_e[e_minus];
                    x[0] += -quad_points[q]*normal[1]*area;
                    x[1] += quad_points[q]*normal[0]*area;
                    x_quad = x_e[e_minus] + x;                    //also need exact quadrature point location

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;
                    u_minus = p_minus/rho_minus;

                    //relative position to centroid of B
                    x = x_f[f] - x_e[e_plus];
                    x[0] += -quad_points[q]*normal[1]*area;
                    x[1] += quad_points[q]*normal[0]*area;

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;
                    u_plus = p_plus/rho_plus;

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = std::sqrt(rho_plus*rho_minus);
                    theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                    c = driver->fluid_material->getSpeedOfSound(job, driver, x_quad, rho_bar, theta_bar);

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    if (lambda_1 < delta*c){
                        lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                    }
                    if (lambda_2 < delta*c){
                        lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                    }

                    //calculate Roe eigenvector coefficients
                    a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                    a_2 = (rho_plus - rho_minus) - a_1;
                    a_3 = p_plus - p_minus - (rho_plus-rho_minus)*u_bar;
                    a_3 = a_3 - a_3.dot(normal)*normal; //remove normal component of a_3 vector

                    /*
                    std::cout << "alpha_int: " << a_1 << " " << a_2 << std::endl;
                    std::cout << "quadrature_location: " << x[0] << ", " << x[1] << std::endl;
                    std::cout << "rho: " << rho_plus << ", " << rho_minus << std::endl;
                    std::cout << ", u: " << u_plus[0] << ", " << u_minus[0] << std::endl;
                    std::cout << "u_bar: " << u_bar[0] << ", c_bar: " << c << std::endl;
                     */

                    //flux in n direction
                    flux = area/num_quad * 0.5 * (p_plus*u_plus.dot(normal) + p_minus*u_minus.dot(normal)
                                                  - a_1*lambda_1*(u_bar - c*normal)
                                                  - a_2*lambda_2*(u_bar + c*normal)
                                                  - a_3*std::abs(u_bar.dot(normal)));

                    //add tractions to momentum flux
                    //for now use simple reconstruction of theta
                    tau_minus = driver->fluid_material->getShearStress(job, driver, x_quad, L[e_minus], rho_bar, theta_bar);
                    tau_plus  = driver->fluid_material->getShearStress(job, driver, x_quad, L[e_plus], rho_bar, theta_bar);
                    P_plus = c*c*(rho_plus - rho_bar) + driver->fluid_material->getPressure(job, driver, x_quad, rho_bar, theta_bar);
                    P_minus = c*c*(rho_minus - rho_plus) + P_plus;

                    flux += area/num_quad * 0.5 * ((P_plus + P_minus)*normal - KinematicVector((tau_plus + tau_minus)*normal, job->JOB_TYPE));

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else if (bc_tags[f] == DIRICHLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    if (e_minus > -1) {
                        //relative position to centroid of A
                        x = x_f[f] - x_e[e_minus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        x_quad = x_e[e_minus] + x;

                        //calculate A properties
                        p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;

                        //tractions
                        rho_bar = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        theta_bar = driver->fluid_body->theta(e_minus);

                        //estimate L
                        x_0 = getElementCentroid(job, e_minus);
                        for (int ii=0; ii<GRID_DIM; ii++){
                            for (int jj=0; jj<GRID_DIM; jj++){
                                L_tmp(ii,jj) = (bc_values[f][ii] - driver->fluid_body->p[e_minus][ii]/rho_bar)
                                               /(x_face - x_0).dot(normal)*normal[jj];
                            }
                        }

                        tau_minus = driver->fluid_material->getShearStress(job, driver, x_quad, L_tmp, rho_bar, theta_bar);
                        P_minus = driver->fluid_material->getPressure(job, driver, x_quad, rho_bar, theta_bar);

                        flux = p_minus*area/num_quad*bc_values[f].dot(normal)
                               + area/num_quad*P_minus*normal
                               - area/num_quad*KinematicVector(tau_minus*normal, job->JOB_TYPE);

                        result(e_minus) -= flux;
                    }

                    if (e_plus > -1) {
                        //relative position to centroid of B
                        x = x_f[f] - x_e[e_plus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        x_quad = x_e[e_plus] + x;

                        //calculate B properties
                        p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;

                        //tractions
                        rho_bar = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        theta_bar = driver->fluid_body->theta(e_plus);

                        //estimate L
                        x_0 = getElementCentroid(job, e_plus);
                        for (int ii=0; ii<GRID_DIM; ii++){
                            for (int jj=0; jj<GRID_DIM; jj++){
                                L_tmp(ii,jj) = (bc_values[f][ii] - driver->fluid_body->p[e_plus][ii]/rho_bar)
                                               /(x_face - x_0).dot(normal)*normal[jj];
                            }
                        }

                        tau_plus = driver->fluid_material->getShearStress(job, driver, x_quad, L_tmp, rho_bar, theta_bar);
                        P_plus = driver->fluid_material->getPressure(job, driver, x_quad, rho_bar, theta_bar);

                        flux = p_plus*area/num_quad*bc_values[f].dot(normal)
                               + area/num_quad*P_plus*normal
                               - area/num_quad*KinematicVector(tau_plus*normal, job->JOB_TYPE);

                        result(e_plus) += flux;
                    }
                }
            } else if (bc_tags[f] == SUPERSONIC_INLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    if (e_minus > -1) {
                        //relative position to centroid of A
                        x = x_f[f] - x_e[e_minus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        x_quad = x_e[e_minus] + x;

                        //calculate A properties
                        p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;

                        //tractions
                        rho_bar = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        theta_bar = driver->fluid_body->theta(e_minus);

                        //estimate L
                        x_0 = getElementCentroid(job, e_minus);
                        for (int ii=0; ii<GRID_DIM; ii++){
                            for (int jj=0; jj<GRID_DIM; jj++){
                                L_tmp(ii,jj) = (bc_values[f][ii] - driver->fluid_body->p[e_minus][ii]/rho_bar)
                                               /(x_face - x_0).dot(normal)*normal[jj];
                            }
                        }

                        tau_minus = driver->fluid_material->getShearStress(job, driver, x_quad, L_tmp, bc_density[f], theta_bar);
                        P_minus = driver->fluid_material->getPressure(job, driver, x_quad, bc_density[f], theta_bar);

                        flux = bc_density[f]*bc_values[f]*area/num_quad*bc_values[f].dot(normal)
                               + area/num_quad*P_minus*normal
                               - area/num_quad*KinematicVector(tau_minus*normal, job->JOB_TYPE);

                        result(e_minus) -= flux;
                    }

                    if (e_plus > -1) {
                        //relative position to centroid of B
                        x = x_f[f] - x_e[e_plus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        x_quad = x_e[e_plus] + x;

                        //calculate B properties
                        p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;

                        //tractions
                        rho_bar = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        theta_bar = driver->fluid_body->theta(e_plus);

                        //estimate L
                        x_0 = getElementCentroid(job, e_plus);
                        for (int ii=0; ii<GRID_DIM; ii++){
                            for (int jj=0; jj<GRID_DIM; jj++){
                                L_tmp(ii,jj) = (bc_values[f][ii] - driver->fluid_body->p[e_plus][ii]/rho_bar)
                                               /(x_face - x_0).dot(normal)*normal[jj];
                            }
                        }

                        tau_plus = driver->fluid_material->getShearStress(job, driver, x_quad, L_tmp, bc_density[f], theta_bar);
                        P_plus = driver->fluid_material->getPressure(job, driver, x_quad, bc_density[f], theta_bar);

                        flux = bc_density[f]*bc_values[f]*area/num_quad*bc_values[f].dot(normal)
                               + area/num_quad*P_plus*normal
                               - area/num_quad*KinematicVector(tau_plus*normal, job->JOB_TYPE);

                        result(e_plus) += flux;
                    }
                }
            } else if (bc_tags[f] == NEUMANN){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    if (e_minus > -1) {
                        //relative position to centroid of A
                        x = x_f[f] - x_e[e_minus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        x_quad = x_e[e_minus] + x;

                        //calculate A properties
                        rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;
                        u_minus = p_minus/rho_minus;
                        flux = p_minus*area/num_quad*u_minus.dot(normal);

                        //add traction directly
                        result(e_minus) -= flux;
                        result(e_minus) += area/num_quad*bc_values[f];
                    }

                    if (e_plus > -1) {
                        //relative position to centroid of B
                        x = x_f[f] - x_e[e_plus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        x_quad = x_e[e_plus] + x;

                        //calculate B properties
                        rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;
                        u_plus = p_plus/rho_plus;
                        flux = p_plus*area/num_quad*u_plus.dot(normal);

                        //add traction directly
                        result(e_plus) += flux;
                        result(e_minus) += area/num_quad*bc_values[f];
                    }
                }
            } else if (bc_tags[f] == NEUMANN_DAMPING){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    if (e_minus > -1) {
                        //relative position to centroid of A
                        x = x_f[f] - x_e[e_minus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        x_quad = x_e[e_minus] + x;

                        //calculate A properties
                        rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;
                        u_minus = p_minus/rho_minus;
                        flux = p_minus*area/num_quad*u_minus.dot(normal);

                        theta_bar = driver->fluid_body->theta(e_minus);
                        P_minus = driver->fluid_material->getPressure(job, driver, x_quad, driver->fluid_body->rho(e_minus), theta_bar);

                        //add traction directly
                        result(e_minus) -= flux + lambda*area/num_quad*P_minus*normal;
                        result(e_minus) += (1-lambda)*area/num_quad*bc_values[f];
                    }

                    if (e_plus > -1) {
                        //relative position to centroid of B
                        x = x_f[f] - x_e[e_plus];
                        x[0] += -quad_points[q]*normal[1]*area;
                        x[1] += quad_points[q]*normal[0]*area;
                        x_quad = x_e[e_plus] + x;

                        //calculate B properties
                        rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;
                        u_plus = p_plus/rho_plus;
                        flux = p_plus*area/num_quad*u_plus.dot(normal);

                        theta_bar = driver->fluid_body->theta(e_plus);
                        P_plus = driver->fluid_material->getPressure(job, driver, x_quad, driver->fluid_body->rho(e_plus), theta_bar);

                        //add traction directly
                        result(e_plus) += flux + lambda*area/num_quad*P_plus*normal;
                        result(e_minus) += (1-lambda)*area/num_quad*bc_values[f];
                    }
                }
            }
        }
    }
    return result;
}

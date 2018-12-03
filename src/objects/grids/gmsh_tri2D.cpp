//
// Created by aaron on 5/25/18.
// gmsh_tri2D.cpp
//


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <math.h>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "grids.hpp"

/*----------------------------------------------------------------------------*/
//
int TriangularGridLinear::ijk_to_n(Job* job, const Eigen::VectorXi& ijk){
    int idOUT = 0;
    for (int pos=0;pos<ijk.size();pos++){
        if (pos == 2){
            idOUT += ijk(pos) * Nx(pos-1) * Nx(pos-2);
        } else if (pos == 1){
            idOUT += ijk(pos) * Nx(pos-1);
        } else {
            idOUT += ijk(pos);
        }
    }
    return idOUT;
}

Eigen::VectorXi TriangularGridLinear::n_to_ijk(Job* job, int n) {
    Eigen::VectorXi ijk = Eigen::VectorXi(GRID_DIM);//Eigen::VectorXi(job->DIM);
    int tmp = n;
    //find i,j,k representation of node id
    for (int i = 0; i < ijk.rows(); i++) {
        ijk(i) = tmp % (Nx(i));
        tmp = tmp / (Nx(i));
    }

    return ijk;
}

/*----------------------------------------------------------------------------*/
//
bool TriangularGridLinear::line_segment_intersect(Eigen::VectorXd s0_p0, Eigen::VectorXd s0_p1, Eigen::VectorXd s1_p0, Eigen::VectorXd s1_p1){
    //return true if line segments intersect (in 2-space)

    //solve for ax + by = c
    Eigen::Matrix2d ab_mat;
    Eigen::Vector2d c_vec, xy_vec;
    if (s0_p0(0) == s0_p1(0)){
        ab_mat(0,1) = 0;
        ab_mat(0,0) = 1;
        c_vec(0) = s0_p0(0);
    } else if (s0_p0(1) == s0_p1(1)){
        ab_mat(0,1) = 1;
        ab_mat(0,0) = 0;
        c_vec(0) = s0_p0(1);
    } else {
        ab_mat(0,1) = 1;
        ab_mat(0,0) = (s0_p1(1) - s0_p0(1))/(s0_p0(0) - s0_p1(0));
        c_vec(0) = ab_mat(0,0)*s0_p0(0) + s0_p0(1);
    }

    if (s1_p0(0) == s1_p1(0)){
        ab_mat(1,1) = 0;
        ab_mat(1,0) = 1;
        c_vec(1) = s1_p0(0);
    } else if (s1_p0(1) == s1_p1(1)){
        ab_mat(1,1) = 1;
        ab_mat(1,0) = 0;
        c_vec(1) = s1_p0(1);
    } else {
        ab_mat(1,1) = 1;
        ab_mat(1,0) = (s1_p1(1) - s1_p0(1))/(s1_p0(0) - s1_p1(0));
        c_vec(1) = ab_mat(1,0)*s1_p0(0) + s1_p0(1);
    }

    //solve for x,y intercept
    if (ab_mat(0,0) == ab_mat(1,0)){
        //parallel lines
        return false;
    }

    xy_vec = ab_mat.inverse() * c_vec;

    //check that x,y within bounds
    double xx00, xx01, xx10, xx11;
    xx00 = std::sqrt((xy_vec(0) - s0_p0(0))*(xy_vec(0) - s0_p0(0)) + (xy_vec(1) - s0_p0(1))*(xy_vec(1) - s0_p0(1))); //(xy_vec - s0_p0).norm();
    xx01 = std::sqrt((xy_vec(0) - s0_p1(0))*(xy_vec(0) - s0_p1(0)) + (xy_vec(1) - s0_p1(1))*(xy_vec(1) - s0_p1(1))); //(xy_vec - s0_p1).norm();
    xx10 = std::sqrt((xy_vec(0) - s1_p0(0))*(xy_vec(0) - s1_p0(0)) + (xy_vec(1) - s1_p0(1))*(xy_vec(1) - s1_p0(1))); //(xy_vec - s1_p0).norm();
    xx11 = std::sqrt((xy_vec(0) - s1_p1(0))*(xy_vec(0) - s1_p1(0)) + (xy_vec(1) - s1_p1(1))*(xy_vec(1) - s1_p1(1))); //(xy_vec - s1_p1).norm();

    double sl0, sl1;
    sl0 = std::sqrt((s0_p1(0) - s0_p0(0))*(s0_p1(0) - s0_p0(0)) + (s0_p1(1) - s0_p0(1))*(s0_p1(1) - s0_p0(1)));
    sl1 = std::sqrt((s1_p1(0) - s1_p0(0))*(s1_p1(0) - s1_p0(0)) + (s1_p1(1) - s1_p0(1))*(s1_p1(1) - s1_p0(1)));

    return ( (std::max(xx00,xx01) <= sl0) && (std::max(xx10,xx11) <= sl1) );
}


/*----------------------------------------------------------------------------*/
//
int TriangularGridLinear::whichSearchCell(const KinematicVector &xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "whichSearchCell failed");
    bool inDomain;
    int elementID = floor(xIN[0] / hx[0]); //id in x-dimension
    for (int i = 0; i < GRID_DIM; i++) {
        inDomain = (xIN[i] < Lx[i] && xIN[i] >= 0);
        if (!inDomain) { //if xIn is outside domain, return -1
            return -1;
        }
        if (i == 1) {
            //add number of elements in next layer of id in higher dimensions
            elementID += floor(xIN[i] / hx[i]) * (Nx[i - 1]);
        } else if (i == 2) {
            elementID += floor(xIN[i] / hx[i]) * (Nx[i - 1]) * (Nx[i - 2]);
        }
    }
    return elementID;
}

bool TriangularGridLinear::inElement(Job* job, const KinematicVector& xIN, int idIN){
    //check if point is interior to element
    KinematicVector xi = KinematicVector(job->JOB_TYPE);
    KinematicTensor tmpMat = KinematicTensor(job->JOB_TYPE);

    //xi = Ainv * (x - x0)
    xi = Ainv(idIN)*(xIN - x_n(nodeIDs(idIN,0)));

    return (xi(0) >= 0 && xi(1) >= 0 && (xi(0) + xi(1)) <= 1);
}


/*----------------------------------------------------------------------------*/
//
void TriangularGridLinear::init(Job* job){
    if (GRID_DIM != 2){
        std::cerr << "TriangularGridLinear requires GRID_DIM = 2, got: " << GRID_DIM << std::endl;
        exit(0);
    }

    if (fp64_props.size() < 1|| str_props.size() < 1){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need 1 filename and 1 length scale defined.\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        msh_filename = str_props[0];
        lc = fp64_props[0];

        std::cout << "Grid properties (filename = " << msh_filename << ", lc = " << lc << ")." << std::endl;
    }

    int len;
    std::string line; //read line
    std::vector<std::string> lvec;
    std::ifstream fin(msh_filename); //file to load from
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
                        v_n.resize(len);
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
                            //read nodes in block
                            for (int i=0;i<block_length;i++){
                                std::getline(fin,line);
                                lvec = Parser::splitString(line,' ');
                                //lvec[0] gives gmsh id (1-indexed)
                                x_n(i,0) = std::stod(lvec[1]); //x-coord
                                x_n(i,1) = std::stod(lvec[2]); //y-coord
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
                        int block_length, block_type;
                        while (std::getline(fin, line)){
                            if (line.compare("$EndElements") == 0){
                                break;
                            }

                            lvec = Parser::splitString(line, ' ');
                            //block info
                            block_type = std::stoi(lvec[3]);
                            block_length = std::stoi(lvec[4]);
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
                            } else {
                                //read through block and continue
                                for (int k=0; k<block_length; k++){
                                    std::getline(fin, line);
                                }
                            }

                            element_count = i + 1;
                            nodeIDs.conservativeResize(i + 1, npe);
                        }

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
        std::cerr << "ERROR: Cannot open " << msh_filename << "! Exiting." << std::endl;
        exit(0);
    }

    //call initializing function
    hiddenInit(job);

    std::cout << "Grid Initialized." << std::endl;

    return;
}
void TriangularGridLinear::hiddenInit(Job* job){

    //called from loading or initializing functions

    v_e.resize(nodeIDs.rows());
    v_n.resize(x_n.size());
    v_n.setZero();

    //element volume
    //herons formula
    double a, b, c, s;
    for (int e=0; e<nodeIDs.rows(); e++){
        a = (x_n(nodeIDs(e,0)) - x_n(nodeIDs(e,1))).norm();
        b = (x_n(nodeIDs(e,1)) - x_n(nodeIDs(e,2))).norm();
        c = (x_n(nodeIDs(e,2)) - x_n(nodeIDs(e,0))).norm();
        s = (a+b+c)/2.0;

        v_e(e) = std::sqrt(s*(s-a)*(s-b)*(s-c));

        //nodal volume
        for (int c=0;c<nodeIDs.cols();c++) {
            v_n(nodeIDs(e, c)) += v_e(e) / nodeIDs.cols();
        }
    }

    //setup A matrices
    A = KinematicTensorArray(nodeIDs.rows(), job->JOB_TYPE);
    Ainv = KinematicTensorArray(nodeIDs.rows(), job->JOB_TYPE);


    for (int e=0;e<nodeIDs.rows();e++){
        //A = [-x0+x1, -x0+x2; -y0+y1, -y0+y2]
        A(e,0,0) = -x_n(nodeIDs(e,0),0) + x_n(nodeIDs(e,1),0);
        A(e,0,1) = -x_n(nodeIDs(e,0),0) + x_n(nodeIDs(e,2),0);
        A(e,1,0) = -x_n(nodeIDs(e,0),1) + x_n(nodeIDs(e,1),1);
        A(e,1,1) = -x_n(nodeIDs(e,0),1) + x_n(nodeIDs(e,2),1);

        //Ainv = 1 / (ad - bc) [d, -b; -c, a]
        s = 1 / (A(e,0,0)*A(e,1,1) - A(e,0,1)*A(e,1,0));
        Ainv(e,0,0) = s*A(e,1,1);
        Ainv(e,0,1) = -s*A(e,0,1);
        Ainv(e,1,0) = -s*A(e,1,0);
        Ainv(e,1,1) = s*A(e,0,0);
    }


    //setup search grid
    //store length, number of linear nodes, and deltas
    Lx = KinematicVector(job->JOB_TYPE);
    Lx.setZero();
    for (int i=0;i<x_n.size();i++){
        for (int pos=0;pos<GRID_DIM;pos++){
            if (x_n(i,pos) > Lx(pos)){
                Lx(pos) = x_n(i,pos);
            }
        }
    }
    Nx = Eigen::VectorXi(job->DIM);;
    for (size_t pos=0;pos<GRID_DIM;pos++){
        Nx(pos) = std::floor(Lx(pos)/lc);
    }

    hx = KinematicVector(job->JOB_TYPE);
    for (int pos=0;pos<GRID_DIM;pos++){
        hx(pos) = Lx(pos) / Nx(pos);
    }

    for (int pos=GRID_DIM;pos<hx.size();pos++){
        Lx(pos) = 0;
        Nx(pos) = 0;
        hx(pos) = 0;
    }

    int len = 1; //number of cells in search grid
    for (int pos=0;pos<GRID_DIM;pos++){
        len *= Nx(pos);
    }

    Eigen::VectorXd ijk;
    Eigen::Vector2d s1_p0, s1_p1;
    std::vector<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d>> list_p0;
    std::vector<Eigen::Vector2d,Eigen::aligned_allocator<Eigen::Vector2d>> list_p1;
    Eigen::VectorXd xe0, xe1, xe2;

    s1_p0 << 0,0; s1_p1 << 1,0;
    list_p0.push_back(s1_p0);
    list_p1.push_back(s1_p1);

    s1_p0 << 1,0; s1_p1 << 1,1;
    list_p0.push_back(s1_p0);
    list_p1.push_back(s1_p1);

    s1_p0 << 1,1; s1_p1 << 0,1;
    list_p0.push_back(s1_p0);
    list_p1.push_back(s1_p1);

    s1_p0 << 0,1; s1_p1 << 0,0;
    list_p0.push_back(s1_p0);
    list_p1.push_back(s1_p1);

    for (int i=0;i<len;i++){
        search_offsets.push_back(search_cells.size());
        //fill search_cells
        //check every side for intersect with every edge
        ijk = (n_to_ijk(job,i)).cast<double>();
        for (int side = 0; side < 4; side++){
            s1_p0(0) = (ijk(0) + list_p0[side](0));
            s1_p0(1) = (ijk(1) + list_p0[side](1));
            s1_p1(0) = (ijk(0) + list_p1[side](0));
            s1_p1(1) = (ijk(1) + list_p1[side](1));
            for (int pos = 0; pos < s1_p0.size(); pos++){
                s1_p0[pos] *= hx[pos];
                s1_p1[pos] *= hx[pos];
            }

            for (int e = 0; e<element_count; e++){
                //check corners within cell first
                if (i == whichSearchCell(x_n(nodeIDs(e, 0)))) {
                    search_cells.push_back(e);
                    continue;
                } else if (i == whichSearchCell(x_n(nodeIDs(e, 1)))) {
                    search_cells.push_back(e);
                    continue;
                } else if (i == whichSearchCell(x_n(nodeIDs(e, 2)))) {
                    search_cells.push_back(e);
                    continue;
                }

                //check 3-edges for intersect
                xe0 = EIGEN_MAP_OF_KINEMATIC_VECTOR(x_n(nodeIDs(e,0)));
                xe1 = EIGEN_MAP_OF_KINEMATIC_VECTOR(x_n(nodeIDs(e,1)));
                xe2 = EIGEN_MAP_OF_KINEMATIC_VECTOR(x_n(nodeIDs(e,2)));


                if (line_segment_intersect(xe0, xe1, s1_p0, s1_p1)){
                    search_cells.push_back(e);
                    continue;
                } else if (line_segment_intersect(xe1, xe2, s1_p0, s1_p1)){
                    search_cells.push_back(e);
                    continue;
                } else if (line_segment_intersect(xe2, xe0, s1_p0, s1_p1)){
                    search_cells.push_back(e);
                    continue;
                }
            }
        }
        std::cout << "Initializing... " << i << "/" << len << "\r";
    }
    search_offsets.push_back(search_cells.size());

    return;
}


/*----------------------------------------------------------------------------*/
//
void TriangularGridLinear::writeFrame(Job* job, Serializer* serializer){
    return;
}

std::string TriangularGridLinear::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

int TriangularGridLinear::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//
void TriangularGridLinear::writeHeader(Job* job, Body* body, Serializer* serializer, std::ofstream& nfile, int SPEC){
    if (SPEC != Serializer::VTK){
        std::cerr << "ERROR: Unknown file SPEC in writeHeader: " << SPEC  << "! Exiting." << std::endl;
        exit(0);
    }

    int nlen = body->nodes->x.size();

    nfile << "ASCII\n";
    nfile << "DATASET UNSTRUCTURED_GRID\n";

    nfile << "POINTS " << nlen << " double\n";
    for (int i=0;i<nlen;i++){
        //vtk files require x,y,z
        for (int pos = 0; pos < 3; pos++){
            if (pos < GRID_DIM && (body->nodes->active(i) != 0) && std::isfinite(body->nodes->x(i,pos))){
                nfile << body->nodes->x(i,pos) << " ";
            } else {
                nfile << "0 ";
            }
        }
        nfile << "\n";
    }

    nfile << "POINTS " << node_count << " double\n";
    for (int i=0;i<node_count;i++){
        //vtk requires 3D position
        nfile << x_n(i,0) << " " << x_n(i,1) << " " << 0 << "\n";
    }
    nfile << "CELLS " << element_count << " " << 4*element_count << "\n";
    for (int e=0;e<element_count;e++){
        nfile << "3 " << nodeIDs(e,0) << " " << nodeIDs(e,1) << " " << nodeIDs(e,2) << "\n";
    }

    nfile << "CELL_TYPES " << element_count << "\n";
    for (int e=0;e<element_count;e++){
        nfile << "5\n";
    }
    nfile << "POINT_DATA " << nlen << "\n";
} //write cell types


/*----------------------------------------------------------------------------*/
//
int TriangularGridLinear::whichElement(Job* job, KinematicVector& xIN) {
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "whichElement failed");
    //find search cell
    int searchID = whichSearchCell(xIN);
    if (searchID < 0) {
        return -1;
    }

    //search in that cell
    int elementID;
    for (int cellID = search_offsets[searchID]; cellID < search_offsets[searchID + 1]; cellID++) {
        elementID = search_cells[cellID];
        if (inElement(job, xIN, elementID)) {
            return elementID;
        }
    }
    return -1;
}

bool TriangularGridLinear::inDomain(Job* job, KinematicVector& xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "inDomain failed");

    if (whichSearchCell(xIN) == -1){
        return false;
    } else {
        return true;
    }
}

KinematicVector TriangularGridLinear::nodeIDToPosition(Job* job, int idIN){
    return x_n(idIN);
}

/*----------------------------------------------------------------------------*/
//
void TriangularGridLinear::evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL){
    int elementID = whichElement(job,xIN);
    if (elementID < 0){
        return;
    }

    KinematicVector xi = KinematicVector(job->JOB_TYPE);

    //xi = Ainv * (x - x0)
    xi = Ainv(elementID)*(xIN - x_n(nodeIDs(elementID,0)));

    //fill shape function vals
    double tmp;
    tmp = 1 - xi(0) - xi(1);
    nID.push_back(nodeIDs(elementID,0));
    nVAL.push_back(tmp);

    tmp = xi(0);
    nID.push_back(nodeIDs(elementID,1));
    nVAL.push_back(tmp);

    tmp = xi(1);
    nID.push_back(nodeIDs(elementID,2));
    nVAL.push_back(tmp);

    return;
}

void TriangularGridLinear::evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD){
    int elementID = whichElement(job,xIN);
    if (elementID < 0){
        return;
    }

    KinematicVector xi = KinematicVector(job->JOB_TYPE);
    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    KinematicTensor map = A(elementID);
    KinematicTensor inv_map = Ainv(elementID);

    //xi = Ainv * (x - x0)
    xi = inv_map*(xIN - x_n(nodeIDs(elementID,0)));

    //fill shape function vals
    tmpVec(0) = -1; tmpVec(1) = -1;
    tmpVec = map*tmpVec;
    nID.push_back(nodeIDs(elementID,0));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 1; tmpVec(1) = 0;
    tmpVec = map*tmpVec;
    nID.push_back(nodeIDs(elementID,1));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 0; tmpVec(1) = 1;
    tmpVec = map*tmpVec;
    nID.push_back(nodeIDs(elementID,2));
    nGRAD.push_back(tmpVec);

    return;
}

double TriangularGridLinear::nodeVolume(Job* job, int idIN){
    return v_n(idIN);
}

double TriangularGridLinear::elementVolume(Job* job, int idIN){
    return v_e(idIN);
}

int TriangularGridLinear::nodeTag(Job* job, int idIN){
    return -1;
}

double TriangularGridLinear::nodeSurfaceArea(Job *job, int idIN) {
    //need to implement this at some point
    std::cout << "WARNING: TriangularGridLinear::nodeSurfaceArea() isn't implemented!" << std::endl;
    return 0;
}
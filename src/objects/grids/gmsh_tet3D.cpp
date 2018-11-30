//
// Created by aaron on 11/30/18.
// gmsh_tet3D.cpp
//

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
int TetrahedralGridLinear::whichSearchCell(const KinematicVector &xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "whichSearchCell failed; invalid KinematicVector TYPE");
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

bool TetrahedralGridLinear::inElement(Job* job, const KinematicVector& xIN, int idIN){
    //check if point is interior to element
    KinematicVector xi = KinematicVector(job->JOB_TYPE);
    KinematicTensor tmpMat = KinematicTensor(job->JOB_TYPE);

    //xi = Ainv * (x - x0)
    xi = Ainv(idIN)*(xIN - x_n(nodeIDs(idIN,0)));

    return (xi(0) >= 0 && xi(1) >= 0 && xi(2) >= 0 && (xi(0) + xi(1) + xi(2)) <= 1);
}


/*----------------------------------------------------------------------------*/
//
void TetrahedralGridLinear::init(Job* job){
    if (GRID_DIM != 3){
        std::cerr << "TetrahedralGridLinear requires GRID_DIM = 3, got: " << GRID_DIM << std::endl;
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
        std::vector<std::string> headers= {"$Nodes","$Elements"};
        while (std::getline(fin,line)) {
            switch(Parser::findStringID(headers,line)){
                case 0:
                    //nodes
                    std::getline(fin,line);
                    len = std::stoi(line); //number of nodes;

                    //nodeTags.resize(len);
                    x_n = KinematicVectorArray(len,job->JOB_TYPE);
                    v_n.resize(len);
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
                        v_e.resize(len); //element volume

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
        std::cerr << "ERROR: Cannot open " << msh_filename << "! Exiting." << std::endl;
        exit(0);
    }

    //call initializing function
    hiddenInit(job);

    std::cout << "Grid Initialized." << std::endl;

    return;
}
void TetrahedralGridLinear::hiddenInit(Job* job){

    //called from loading or initializing functions

    v_e.resize(nodeIDs.rows());
    v_n.resize(x_n.size());
    v_n.setZero();

    //element volume
    KinematicVector a, b, c;
    for (int e=0; e<nodeIDs.rows(); e++){
        a = x_n(nodeIDs(e,1)) - x_n(nodeIDs(e,0));
        b = x_n(nodeIDs(e,2)) - x_n(nodeIDs(e,0));
        c = x_n(nodeIDs(e,3)) - x_n(nodeIDs(e,0));

        v_e(e) = (a.dot((b.cross(c))))/6.0;

        //nodal volume
        for (int c=0;c<nodeIDs.cols();c++) {
            v_n(nodeIDs(e, c)) += v_e(e) / nodeIDs.cols();
        }
    }

    //setup A matrices
    A = KinematicTensorArray(nodeIDs.rows(), job->JOB_TYPE);
    Ainv = KinematicTensorArray(nodeIDs.rows(), job->JOB_TYPE);

    double s;
    for (int e=0;e<nodeIDs.rows();e++){
        //A
        A(e,0,0) = x_n(nodeIDs(e,1),0) - x_n(nodeIDs(e,0),0);
        A(e,0,1) = x_n(nodeIDs(e,1),1) - x_n(nodeIDs(e,0),0);
        A(e,0,2) = x_n(nodeIDs(e,1),2) - x_n(nodeIDs(e,0),0);

        A(e,1,0) = x_n(nodeIDs(e,2),0) - x_n(nodeIDs(e,0),0);
        A(e,1,1) = x_n(nodeIDs(e,2),1) - x_n(nodeIDs(e,0),0);
        A(e,1,2) = x_n(nodeIDs(e,2),2) - x_n(nodeIDs(e,0),0);

        A(e,2,0) = x_n(nodeIDs(e,3),0) - x_n(nodeIDs(e,0),0);
        A(e,2,1) = x_n(nodeIDs(e,3),1) - x_n(nodeIDs(e,0),0);
        A(e,2,2) = x_n(nodeIDs(e,3),2) - x_n(nodeIDs(e,0),0);

        //Ainv
        Ainv(e) = A(e).inverse();
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


    Eigen::VectorXi ijk;
    int sc;

    //form min/max list
    element_min_max.resize(6*element_count);
    for (int e=0; e<element_count; e++){
        //find search cell
        sc = whichSearchCell(x_n(nodeIDs(e,0)));
        //convert to ijk coords
        ijk = n_to_ijk(job, sc);

        //assign min/max from first node
        for (int pos = 0; pos < ijk.rows(); pos++) {
            element_min_max[6 * e + pos] = ijk[pos];
            element_min_max[6 * e + pos + 3] = ijk[pos];
        }

        //repeat for other 3 nodes
        for (int i=1; i<npe; i++){
            //find search cell
            sc = whichSearchCell(x_n(nodeIDs(e,i)));
            //convert to ijk coords
            ijk = n_to_ijk(job, sc);
            for (int pos=0; pos<ijk.rows(); pos++){
                //if ijk < min, assign to min
                if (ijk(pos) < element_min_max[6*e + pos]){
                    element_min_max[6*e + pos] = ijk[pos];
                } else if (ijk(pos) > element_min_max[6*e + 3 + pos]) {
                    element_min_max[6*e + 3 + pos] = ijk[pos];
                }
            }
        }
    }

    //insert elements into search cell list
    for (int i=0;i<len;i++){
        search_offsets.push_back(search_cells.size());
        //find ijk of search cell
        ijk = n_to_ijk(job, i);
        //fill cell
        for (int e = 0; e<element_count; e++){
            for (int pos = 0; pos < ijk.rows(); pos++){
                //insert element if search cell is within min/max range
                if (ijk(pos) >= element_min_max[6*e + pos] || ijk(pos) <= element_min_max[6*e + 3 + pos]){
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
void TetrahedralGridLinear::writeHeader(Job* job, Body* body, Serializer* serializer, std::ofstream& nfile, int SPEC){
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
        nfile << x_n(i,0) << " " << x_n(i,1) << " " << x_n(i,2) << "\n";
    }
    nfile << "CELLS " << element_count << " " << 4*element_count << "\n";
    for (int e=0;e<element_count;e++){
        nfile << "4 " << nodeIDs(e,0) << " " << nodeIDs(e,1) << " " << nodeIDs(e,2) << " " << nodeIDs(e,3) << "\n";
    }

    nfile << "CELL_TYPES " << element_count << "\n";
    for (int e=0;e<element_count;e++){
        nfile << "10\n";
    }
    nfile << "POINT_DATA " << nlen << "\n";
} //write cell types


/*----------------------------------------------------------------------------*/
//
int TetrahedralGridLinear::whichElement(Job* job, KinematicVector& xIN) {
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

bool TetrahedralGridLinear::inDomain(Job* job, KinematicVector& xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "inDomain failed");

    if (whichSearchCell(xIN) == -1){
        return false;
    } else {
        return true;
    }
}

KinematicVector TetrahedralGridLinear::nodeIDToPosition(Job* job, int idIN){
    return x_n(idIN);
}

/*----------------------------------------------------------------------------*/
//
void TetrahedralGridLinear::evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL){
    int elementID = whichElement(job,xIN);
    if (elementID < 0){
        return;
    }

    KinematicVector xi = KinematicVector(job->JOB_TYPE);

    //xi = Ainv * (x - x0)
    xi = Ainv(elementID)*(xIN - x_n(nodeIDs(elementID,0)));

    //fill shape function vals
    double tmp;
    tmp = 1 - xi(0) - xi(1) - xi(2);
    nID.push_back(nodeIDs(elementID,0));
    nVAL.push_back(tmp);

    tmp = xi(0);
    nID.push_back(nodeIDs(elementID,1));
    nVAL.push_back(tmp);

    tmp = xi(1);
    nID.push_back(nodeIDs(elementID,2));
    nVAL.push_back(tmp);

    tmp = xi(2);
    nID.push_back(nodeIDs(elementID,3));
    nVAL.push_back(tmp);

    return;
}

void TetrahedralGridLinear::evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD){
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
    tmpVec(0) = -1; tmpVec(1) = -1; tmpVec(2) = -1;
    tmpVec = map*tmpVec;
    nID.push_back(nodeIDs(elementID,0));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 1; tmpVec(1) = 0; tmpVec(2) = 0;
    tmpVec = map*tmpVec;
    nID.push_back(nodeIDs(elementID,1));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 0; tmpVec(1) = 1; tmpVec(2) = 0;
    tmpVec = map*tmpVec;
    nID.push_back(nodeIDs(elementID,2));
    nGRAD.push_back(tmpVec);

    tmpVec(0) = 0; tmpVec(1) = 0; tmpVec(2) = 1;
    tmpVec = map*tmpVec;
    nID.push_back(nodeIDs(elementID,3));
    nGRAD.push_back(tmpVec);

    return;
}

double TetrahedralGridLinear::nodeVolume(Job* job, int idIN){
    return v_n(idIN);
}

double TetrahedralGridLinear::elementVolume(Job* job, int idIN){
    return v_e(idIN);
}

int TetrahedralGridLinear::nodeTag(Job* job, int idIN){
    return -1;
}

double TetrahedralGridLinear::nodeSurfaceArea(Job *job, int idIN) {
    //need to implement this at some point
    std::cout << "WARNING: TetrahedralGridLinear::nodeSurfaceArea() isn't implemented!" << std::endl;
    return 0;
}
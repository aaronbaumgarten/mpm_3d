//
// Created by aaron on 7/5/17.
// gmsh_tri2D.cpp
//

//
// Created by aaron on 5/24/17.
// cartesian.cpp
//


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "stringparser.hpp"

#include "job.hpp"
#include "serializer.hpp"

#include "grid.hpp"

#include "body.hpp"
#include "nodes.hpp"

double lc; //characteristic length scale
std::string msh_filename;
Eigen::VectorXd Lx(0,1); //linear dimensions
Eigen::VectorXi Nx(0,1); //linear elements (not nodes)
Eigen::VectorXd hx(0,1);
Eigen::MatrixXi nodeIDs(0,0); //element to node map
//Eigen::MatrixXi nodeTags(0,1); //node tags

size_t npe = 3; //nodes per element
Eigen::MatrixXd x_n(0,0);

Eigen::VectorXd v_n(0,1); //nodal volume
Eigen::VectorXd v_e(0,1); //element volume

Eigen::MatrixXd A(0,0); //map xi to x
Eigen::MatrixXd Ainv(0,0); //map x to xi

std::vector<int> search_cells; //search grid (cell to element map)
std::vector<int> search_offsets; //search offsets


extern "C" void gridInit(Job* job);

extern "C" void gridWriteFrame(Job* job, Serializer* serializer);
extern "C" std::string gridSaveState(Job* job, Serializer* serializer,std::string filepath);
extern "C" int gridLoadState(Job* job, Serializer* serializer, std::string fullpath);

extern "C" int gridWhichElement(Job* job, Eigen::VectorXd xIN);
extern "C" bool gridInDomain(Job* job, Eigen::VectorXd xIN);
extern "C" Eigen::VectorXd gridNodeIDToPosition(Job* job, int idIN);
extern "C" void gridEvaluateShapeFnValue(Job* job, Eigen::VectorXd xIN, std::vector<int>& nID, std::vector<double>& nVAL);
extern "C" void gridEvaluateShapeFnGradient(Job* job, Eigen::VectorXd xIN, std::vector<int>& nID, std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& nGRAD);
extern "C" double gridNodalVolume(Job* job, int idIN);
extern "C" double gridElementVolume(Job* job, int idIN);
extern "C" double gridNodalTag(Job* job, int idIN);

/*----------------------------------------------------------------------------*/

int ijk_to_n(Eigen::VectorXd ijk){
    int idOUT = 0;
    for (size_t pos=0;pos<ijk.rows();pos++){
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

Eigen::VectorXi n_to_ijk(Job* job, int idIN){
    Eigen::VectorXi ijk = job->jobVector<int>();
    int tmp = idIN;
    //find i,j,k representation of node id
    for (size_t i=0;i<ijk.rows();i++){
        ijk(i) = tmp % (Nx(i));
        tmp = tmp/(Nx(i));
    }
    return ijk;
}

/*----------------------------------------------------------------------------*/

bool line_segment_intersect(Eigen::VectorXd s0_p0, Eigen::VectorXd s0_p1, Eigen::VectorXd s1_p0, Eigen::VectorXd s1_p1){
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
    xx00 = (xy_vec - s0_p0).norm();
    xx01 = (xy_vec - s0_p1).norm();
    xx10 = (xy_vec - s1_p0).norm();
    xx11 = (xy_vec - s1_p1).norm();

    return ( (std::max(xx00,xx01) <= (s0_p1 - s0_p0).norm()) && (std::max(xx10,xx11) <= (s1_p1 - s1_p0).norm()) );
}

/*----------------------------------------------------------------------------*/

int whichSearchCell(Eigen::VectorXd xIN) {
    bool inDomain;
    int elementID = floor(xIN[0] / hx[0]); //id in x-dimension
    for (size_t i = 0; i < xIN.size(); i++) {
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

/*----------------------------------------------------------------------------*/

bool inElement(Job* job, Eigen::VectorXd xIN, int idIN){
    //check if point is interior to element
    Eigen::VectorXd xi = job->jobVector<double>();
    Eigen::VectorXd tmpVec = job->jobVector<double>();

    tmpVec = Ainv.row(idIN);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> inv_map = job->jobTensor<double>(tmpVec.data());

    //xi = Ainv * (x - x0)
    xi = inv_map*(xIN - x_n.row(nodeIDs(idIN,0)).transpose());

    return (xi(0) >= 0 && xi(1) >= 0 && (xi(0) + xi(1)) <= 1);
}

/*----------------------------------------------------------------------------*/

void hiddenInit(Job* job){
    //called from loading or initializing functions

    v_e.resize(nodeIDs.rows());
    v_n.resize(x_n.rows());
    v_n.setZero();

    //element volume
    //herons formula
    double a, b, c, s;
    for (size_t e=0; e<nodeIDs.rows(); e++){
        a = (x_n.row(nodeIDs(e,0)) - x_n.row(nodeIDs(e,1))).norm();
        b = (x_n.row(nodeIDs(e,1)) - x_n.row(nodeIDs(e,2))).norm();
        c = (x_n.row(nodeIDs(e,2)) - x_n.row(nodeIDs(e,0))).norm();
        s = (a+b+c)/2.0;

        v_e(e) = std::sqrt(s*(s-a)*(s-b)*(s-c));

        //nodal volume
        for (size_t c=0;c<nodeIDs.cols();c++) {
            v_n(nodeIDs(e, c)) += v_e(e) / nodeIDs.cols();
        }
    }

    //setup A matrices
    A = job->jobTensorArray<double>(nodeIDs.rows());
    Ainv = job->jobTensorArray<double>(nodeIDs.rows());


    for (size_t e=0;e<nodeIDs.rows();e++){
        //A = [-x0+x1, -x0+x2; -y0+y1, -y0+y2]
        A(e,0) = -x_n(nodeIDs(e,0),0) + x_n(nodeIDs(e,1),0);
        A(e,1) = -x_n(nodeIDs(e,0),0) + x_n(nodeIDs(e,2),0);
        A(e,2) = -x_n(nodeIDs(e,0),1) + x_n(nodeIDs(e,1),1);
        A(e,3) = -x_n(nodeIDs(e,0),1) + x_n(nodeIDs(e,2),1);

        //Ainv = 1 / (ad - bc) [d, -b; -c, a]
        s = 1 / (A(e,0)*A(e,3) - A(e,1)*A(e,2));
        Ainv(e,0) = s*A(e,3);
        Ainv(e,1) = -s*A(e,1);
        Ainv(e,2) = -s*A(e,2);
        Ainv(e,3) = s*A(e,0);
    }


    //setup search grid
    Lx = x_n.colwise().maxCoeff();
    Nx = job->jobVector<int>();
    for (size_t pos=0;pos<Nx.rows();pos++){
        Nx(pos) = floor(Lx(pos)/lc);
    }
    hx = Lx.array() / Nx.cast<double>().array();

    size_t len = 1; //number of cells in search grid
    for (size_t pos=0;pos<Nx.rows();pos++){
        len *= Nx(pos);
    }

    Eigen::VectorXd ijk;
    Eigen::Vector2d s1_p0, s1_p1;
    std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>> list_p0;
    std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>> list_p1;

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

    for (size_t i=0;i<len;i++){
        search_offsets.push_back(search_cells.size());
        //fill search_cells
        //check every side for intersect with every edge
        ijk = (n_to_ijk(job,i)).cast<double>();
        for (size_t side = 0; side < 4; side++){
            s1_p0 = (ijk + list_p0[side]).array()*hx.array();
            s1_p1 = (ijk + list_p1[side]).array()*hx.array();
            for (int e = 0; e<job->grid.element_count; e++){
                //check corners within cell first
                if (i == whichSearchCell(x_n.row(nodeIDs(e, 0)))) {
                    search_cells.push_back(e);
                    continue;
                } else if (i == whichSearchCell(x_n.row(nodeIDs(e, 1)))) {
                    search_cells.push_back(e);
                    continue;
                } else if (i == whichSearchCell(x_n.row(nodeIDs(e, 2)))) {
                    search_cells.push_back(e);
                    continue;
                }

                //check 3-edges for intersect
                if (line_segment_intersect(x_n.row(nodeIDs(e,0)), x_n.row(nodeIDs(e,1)), s1_p0, s1_p1)){
                    search_cells.push_back(e);
                    continue;
                } else if (line_segment_intersect(x_n.row(nodeIDs(e,1)), x_n.row(nodeIDs(e,2)), s1_p0, s1_p1)){
                    search_cells.push_back(e);
                    continue;
                } else if (line_segment_intersect(x_n.row(nodeIDs(e,2)), x_n.row(nodeIDs(e,0)), s1_p0, s1_p1)){
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

void gridInit(Job* job){
    if (job->DIM != 2){
        std::cerr << "gmsh_tri2D.so requires DIM = 2, got: " << job->DIM << std::endl;
        exit(0);
    }

    if (job->grid.fp64_props.size() < 1|| job->grid.str_props.size() < 1){
        std::cout << job->grid.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need 1 filename and 1 length scale defined.\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        msh_filename = job->grid.str_props[0];
        lc = job->grid.fp64_props[0];

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
            switch(StringParser::stringFindStringID(headers,line)){
                case 0:
                    //nodes
                    std::getline(fin,line);
                    len = std::stoi(line); //number of nodes;

                    //nodeTags.resize(len);
                    x_n = job->jobVectorArray<double>(len);
                    v_n.resize(len);
                    job->grid.node_count = len;

                    //loop to read in all nodes
                    for (size_t i=0;i<len;i++){
                        std::getline(fin,line);
                        lvec = StringParser::stringSplitString(line,' ');
                        //lvec[0] gives gmsh id (1-indexed)
                        x_n(i,0) = std::stod(lvec[1]); //x-coord
                        x_n(i,1) = std::stod(lvec[2]); //y-coord
                    }
                    break;
                case 1:
                    //elements
                    if (true) {
                        std::getline(fin, line);
                        len = std::stoi(line); //number of potential elements

                        nodeIDs.resize(len, npe); //element to node map
                        v_e.resize(len); //element volume

                        size_t i = -1;
                        while (std::getline(fin, line)) {
                            if (line.compare("$EndElements") == 0) {
                                break;
                            }

                            lvec = StringParser::stringSplitString(line, ' ');
                            //check that element type is triangle
                            if (std::stoi(lvec[1]) == 2) {
                                i++; //increment i
                                len = std::stoi(lvec[2]); //num tags
                                nodeIDs(i, 0) = std::stoi(lvec[3 + len]) - 1;
                                nodeIDs(i, 1) = std::stoi(lvec[4 + len]) - 1;
                                nodeIDs(i, 2) = std::stoi(lvec[5 + len]) - 1;
                            }
                        }

                        job->grid.element_count = i + 1;
                        nodeIDs.conservativeResize(i + 1, 3);
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

    //open point file
    std::ofstream ffile(msh_filename+".vtk", std::ios::trunc);

    if (ffile.is_open()){
        ffile << "# vtk DataFile Version 3.0\n";
        ffile << "Triangles from Grid\n";
        ffile << "ASCII\n";
        ffile << "DATASET UNSTRUCTURED_GRID\n";

        ffile << "POINTS " << job->grid.node_count << " double\n";
        for (size_t i=0;i<job->grid.node_count;i++){
            //vtk requires 3D position
            ffile << x_n(i,0) << " " << x_n(i,1) << " " << 0 << "\n";
        }
        ffile << "CELLS " << job->grid.element_count << " " << 4*job->grid.element_count << "\n";
        for (size_t e=0;e<job->grid.element_count;e++){
            ffile << "3 " << nodeIDs(e,0) << " " << nodeIDs(e,1) << " " << nodeIDs(e,2) << "\n";
        }

        ffile << "CELL_TYPES " << job->grid.element_count << "\n";
        for (size_t e=0;e<job->grid.element_count;e++){
            ffile << "5\n";
        }

        ffile.close();

        //pfile << "POINT_DATA " << plen << "\n";
        //scalars, vectors and tensors here
    } else {
        std::cerr << "Could not open file: " << msh_filename+".vtk" << " !" << std::endl;
    }

    std::cout << "Grid Initialized." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void gridWriteFrame(Job* job, Serializer* serializer){
    //nothing to report
    return;
}

/*----------------------------------------------------------------------------*/

std::string gridSaveState(Job* job, Serializer* serializer,std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename
    std::stringstream s;
    s << "mpm_v2.grid." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    if (ffile.is_open()){
        ffile << "# mpm_v2 grids/cartesian.so\n";
        ffile << msh_filename << "\n";
        ffile << lc << "\n";
        ffile << x_n.rows() << "\n"; //num nodes
        job->jobVectorArrayToFile(x_n,ffile);
        ffile << nodeIDs.rows() << "\n"; //num elements
        ffile << nodeIDs.cols() << "\n"; //num corners
        for (size_t e=0;e<nodeIDs.rows();e++){
            for (size_t c=0;c<nodeIDs.cols();c++){
                if (c>0){
                    ffile << ' ';
                }
                ffile << nodeIDs(e,c);
            }
            ffile << "\n";
        }
        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filepath+filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Grid Saved." << std::endl;

    return filename;
}

/*----------------------------------------------------------------------------*/

int gridLoadState(Job* job, Serializer* serializer, std::string fullpath){
    int len;

    std::string line;
    std::vector<std::string> lvec;
    std::stringstream ss;
    std::ifstream fin(fullpath);

    if(fin.is_open()){
        std::getline(fin,line); //first line

        std::getline(fin,line); //filename
        msh_filename = line;
        std::getline(fin,line); //lc
        lc = std::stod(line);

        std::getline(fin,line); //num nodes
        len = std::stoi(line);

        x_n = job->jobVectorArray<double>(len);
        job->jobVectorArrayFromFile(x_n,fin); //node position from file
        job->grid.node_count = len;

        std::getline(fin,line); // elements
        len = std::stoi(line);
        job->grid.element_count = len;
        std::getline(fin,line); // corners
        npe = std::stoi(line);

        nodeIDs.resize(len,npe);
        for (size_t e=0;e<len;e++){
            std::getline(fin,line);
            lvec = StringParser::stringSplitString(line,' ');
            for (size_t c=0; c<npe; c++){
                nodeIDs(e,c) = std::stoi(lvec[c]);
            }
        }
        std::cout << "Grid properties (filename = " << msh_filename << ", lc = " << lc << ")." << std::endl;

        //call hidden initializer
        hiddenInit(job);

        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Grid Loaded." << std::endl;
    return 1;
}

/*----------------------------------------------------------------------------*/

int gridWhichElement(Job* job, Eigen::VectorXd xIN){
    //find search cell
    int searchID = whichSearchCell(xIN);
    if (searchID < 0){
        return -1;
    }

    //search in that cell
    int elementID;
    for (int cellID = search_offsets[searchID]; cellID<search_offsets[searchID+1]; cellID++){
        elementID = search_cells[cellID];
        if (inElement(job,xIN,elementID)){
            return elementID;
        }
    }
    return -1;
}

/*----------------------------------------------------------------------------*/

bool gridInDomain(Job* job, Eigen::VectorXd xIN){
    if (whichSearchCell(xIN) == -1){
        return false;
    } else {
        return true;
    }
}

/*----------------------------------------------------------------------------*/

Eigen::VectorXd gridNodeIDToPosition(Job* job, int idIN){
    return x_n.row(idIN).transpose();
}

/*----------------------------------------------------------------------------*/

void gridEvaluateShapeFnValue(Job* job, Eigen::VectorXd xIN, std::vector<int>& nID, std::vector<double>& nVAL){

    int elementID = gridWhichElement(job,xIN);
    if (elementID < 0){
        return;
    }

    Eigen::VectorXd tmpVec;
    Eigen::VectorXd xi = job->jobVector<double>();

    //tmpVec = A.row(elementID);
    //Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> map = job->jobTensor<double>(tmpVec.data());

    tmpVec = Ainv.row(elementID);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> inv_map = job->jobTensor<double>(tmpVec.data());

    //xi = Ainv * (x - x0)
    xi = inv_map*(xIN - x_n.row(nodeIDs(elementID,0)).transpose());

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

/*----------------------------------------------------------------------------*/

void gridEvaluateShapeFnGradient(Job* job, Eigen::VectorXd xIN, std::vector<int>& nID, std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& nGRAD){
    int elementID = gridWhichElement(job,xIN);
    if (elementID < 0){
        return;
    }

    Eigen::VectorXd tmpVec;
    Eigen::VectorXd xi = job->jobVector<double>();

    tmpVec = A.row(elementID);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> map = job->jobTensor<double>(tmpVec.data());

    tmpVec = Ainv.row(elementID);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> inv_map = job->jobTensor<double>(tmpVec.data());

    //xi = Ainv * (x - x0)
    tmpVec = job->jobVector<double>();
    xi = inv_map*(xIN - x_n.row(nodeIDs(elementID,0)).transpose());

    //fill shape function vals
    tmpVec << -1, -1;
    tmpVec = map*tmpVec;
    nID.push_back(nodeIDs(elementID,0));
    nGRAD.push_back(tmpVec);

    tmpVec << 1, 0;
    tmpVec = map*tmpVec;
    nID.push_back(nodeIDs(elementID,1));
    nGRAD.push_back(tmpVec);

    tmpVec << 0, 1;
    tmpVec = map*tmpVec;
    nID.push_back(nodeIDs(elementID,2));
    nGRAD.push_back(tmpVec);

    return;
}

/*----------------------------------------------------------------------------*/

double gridNodalVolume(Job* job, int idIN){
    return v_n(idIN);
}

/*----------------------------------------------------------------------------*/

double gridElementVolume(Job* job, int idIN){
    return v_e(idIN);
}

/*----------------------------------------------------------------------------*/

double gridNodalTag(Job* job, int idIN){
    return -1;
}

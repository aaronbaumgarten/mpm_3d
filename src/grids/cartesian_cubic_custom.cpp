//
// Created by aaron on 9/28/17.
// cartesian_cubic_custom.cpp
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "job.hpp"
#include "serializer.hpp"

#include "grid.hpp"

#include "body.hpp"
#include "nodes.hpp"

Eigen::VectorXd Lx(0,1); //linear dimensions
Eigen::VectorXi Nx(0,1); //linear elements (not nodes)
Eigen::VectorXd hx(0,1);
Eigen::MatrixXi nodeIDs(0,0); //element to node map
Eigen::MatrixXi A(0,0); //0,1 for directions in
size_t npe; //nodes per element
Eigen::MatrixXd x_n(0,0);
Eigen::MatrixXi edge_n(0,0); //0-edge, 1-next to boundary, 2-interior

Eigen::VectorXd v_n(0,1); //nodal volume
double v_e = 1; //element volume

Eigen::VectorXi periodic_props(0,1); //boundary props
Eigen::VectorXi nntoni(0,1); //node number to node id

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

double s(double x, double h){
    if (x < -2*h){
        return 0;
    } else if (x < -h){
        return x*x*x/(6.0 * h*h*h) + x*x/(h*h) + 2.0*x/h + 4.0/3.0;
    } else if (x < 0){
        return -1.0*x*x*x/(2.0*h*h*h) - x*x/(h*h) + 2.0/3.0;
    } else if (x < h){
        return x*x*x/(2.0*h*h*h) - x*x/(h*h) + 2.0/3.0;
    } else if (x < 2*h){
        return -1.0*x*x*x/(6.0*h*h*h) + x*x/(h*h) - 2*x/h + 4.0/3.0;
    } else {
        return 0;
    }
}

double g(double x, double h){
    if (x < -2*h){
        return 0;
    } else if (x < -h){
        return x*x/(2.0 * h*h*h) + 2.0*x/(h*h) + 2.0/h;
    } else if (x < 0){
        return -3.0*x*x/(2.0*h*h*h) - 2.0*x/(h*h);
    } else if (x < h){
        return 3.0*x*x/(2.0*h*h*h) - 2.0*x/(h*h);
    } else if (x < 2*h){
        return -1.0*x*x/(2.0*h*h*h) + 2.0*x/(h*h) - 2.0/h;
    } else {
        return 0;
    }
}

/*----------------------------------------------------------------------------*/

void hiddenInit(Job* job){
    //called from loading or initializing functions
    //needs to have Lx,Nx,hx defined

    //count nodes
    job->grid.node_count = 1;
    job->grid.element_count = 1;
    npe = 1;
    for (size_t i=0;i<Nx.rows();i++){
        job->grid.node_count *= (Nx(i)+1);
        job->grid.element_count *= Nx(i);
        npe *= 4; //4 linear node contributions per element
    }

    x_n = job->jobVectorArray<double>(job->grid.node_count);
    for (size_t i=0;i<x_n.rows();i++){
        x_n.row(i) = (job->grid.gridNodeIDToPosition(job,i)).transpose();
    }

    //initialize A matrix
    //for mapping position in cube to id
    //0 -> -1,-1,-1
    //1 -> +0,-1,-1
    //...
    //64 -> +2,+2,+2
    Eigen::VectorXi i_rel = job->jobVector<int>();
    A = job->jobVectorArray<int>(npe);
    i_rel.setOnes();
    i_rel = -1.0*i_rel;
    for (size_t n=0; n<npe;n++){
        for (size_t i=0;i<i_rel.rows();i++){
            A(n,i) = i_rel(i);
        }
        for (size_t i=0;i<i_rel.rows();i++) {
            if (i_rel(i) < 2){
                i_rel(i) += 1;
                break;
            } else {
                i_rel(i) = -1;
            }
        }
    }

    //setup element to node map
    Eigen::VectorXi ijk = job->jobVector<int>();
    int tmp;
    nodeIDs.resize(job->grid.element_count,npe);
    nodeIDs.setZero();
    for (size_t e=0;e<nodeIDs.rows();e++){
        tmp = e;
        //find i,j,k count for element position
        for (size_t i=0;i<ijk.rows();i++){
            ijk(i) = tmp % Nx(i);
            tmp = tmp / Nx(i);
        }

        //find node ids for element
        for (size_t n=0;n<nodeIDs.cols();n++){
            for (size_t i=0;i<ijk.rows();i++) {
                //n = i + j*imax + k*imax*jmax
                //hardcode
                tmp = ijk(i) + A(n,i); // first guess
                if (tmp < 0 || tmp > Nx(i)) {
                    if (periodic_props(i) == 1 && tmp == -1) {
                        tmp = Nx(i)-1;
                    } else if (periodic_props(i) == 1 && tmp == Nx(i)+1){
                        tmp = 1;
                    } else {
                        nodeIDs(e, n) = -1;
                        break;
                    }
                }

                if (i==0) {
                    nodeIDs(e, n) += tmp;
                } else if (i==1) {
                    nodeIDs(e, n) += tmp * (Nx(i-1)+1);
                } else if (i==2) {
                    nodeIDs(e, n) += tmp * (Nx(i-1)+1) * (Nx(i-2)+1);
                }
            }
        }
    }

    //setup node number to node id map
    nntoni.resize(job->grid.node_count);
    for (size_t i=0; i<job->grid.node_count; i++){
        tmp = i;
        //find i,j,k for node position
        for (size_t i=0;i<ijk.rows();i++){
            ijk(i) = tmp % (Nx(i)+1);
            tmp = tmp / (Nx(i)+1);
        }

        if (job->DIM >= 1 && ijk(0) == Nx(0) && periodic_props(0) == 1){
            ijk(0) = 0;
        }
        if (job->DIM >= 2 && ijk(1) == Nx(1) && periodic_props(1) == 1){
            ijk(1) = 0;
        }
        if (job->DIM == 3 && ijk(2) == Nx(2) && periodic_props(2) == 1){
            ijk(2) = 0;
        }

        tmp = 0;
        //create map
        for (size_t pos=0;pos<ijk.rows();pos++){
            if (pos==0) {
                tmp += ijk(pos);
            } else if (pos==1) {
                tmp += ijk(pos) * (Nx(pos-1) + 1);
            } else if (pos==2) {
                tmp += ijk(pos) * (Nx(pos-2) + 1) * (Nx(pos-1) + 1);
            }
        }
        nntoni(i) = tmp;
    }

    //element volume
    v_e = 1;
    for (size_t pos=0;pos<hx.rows();pos++){
        v_e *= hx(pos);
    }

    //nodal volume and edges
    edge_n = job->jobVectorArray<int>(x_n.rows());
    v_n.resize(x_n.rows());
    v_n.setZero();
    for (size_t n=0;n<x_n.rows();n++){
        tmp = n;
        v_n(n) = 1;
        for (size_t i=0;i<ijk.rows();i++){
            ijk(i) = tmp % (Nx(i)+1);
            tmp = tmp / (Nx(i)+1);

            if ((ijk(i) > 1 && ijk(i) < (Nx(i)-1)) || periodic_props(i) == 1){
                v_n(n) *= hx(i);
                if (ijk(i) == 1) {
                    edge_n(n, i) = 3; //next to periodic boundary at start
                } else if (ijk(i) == Nx(i)-1) {
                    edge_n(n, i) = -3; //next to periodic boundary at end
                } else {
                    edge_n(n, i) = 2; //interior
                }
            } else if (ijk(i) == 1 || ijk(i) == (Nx(i)-1)) {
                v_n(n) *= hx(i)*11.0/12.0;
                if (ijk(i) == 1) {
                    edge_n(n, i) = 1; //next to boundary at start
                } else {
                    edge_n(n, i) = -1; //next to boundary at end
                }
            } else {
                v_n(n) *= hx(i)*7.0/12.0;
                edge_n(n,i) = 0; //boundary
            }
        }
    }

    return;
}

/*----------------------------------------------------------------------------*/

void gridInit(Job* job){
    if (job->grid.fp64_props.size() < job->DIM || job->grid.int_props.size() < 2*job->DIM){
        std::cout << job->grid.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i dimensions and %i properties defined.\n",
                __FILE__, __func__, job->DIM, job->DIM);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        Lx = job->jobVector<double>(job->grid.fp64_props.data());
        Nx = job->jobVector<int>(job->grid.int_props.data());
        hx = Lx.array() / Nx.cast<double>().array();

        periodic_props = job->jobVector<int>();
        for (size_t pos=0;pos<job->DIM;pos++){
            periodic_props(pos) = job->grid.int_props[job->DIM + pos];
        }
        //print grid properties
        std::cout << "Grid properties (Lx = { ";
        for (size_t i=0;i<Lx.rows();i++){
            std::cout << Lx(i) << " ";
        }
        std::cout << "}, Nx = { ";
        for (size_t i=0;i<Nx.rows();i++){
            std::cout << Nx(i) << " ";
        }
        std::cout << "})." << std::endl;
    }

    //call initializing function
    hiddenInit(job);

    std::cout << "Grid Initialized." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void gridWriteFrame(Job* job, Serializer* serializer){
    Eigen::MatrixXd tmp = edge_n.cast<double>();
    serializer->serializerWriteVectorArray(tmp,"edge_id");
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
        ffile << "# mpm_v2 grids/cartesian_cubic_custom.so\n";
        ffile << Lx.transpose() << "\n" << Nx.transpose() << "\n"; //vectors
        ffile << periodic_props.transpose() << "\n";
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
    //job object should be loaded first
    Lx = job->jobVector<double>();
    Nx = job->jobVector<int>();
    hx = job->jobVector<double>();

    std::string line;
    std::stringstream ss;
    std::ifstream fin(fullpath);

    if(fin.is_open()){
        std::getline(fin,line); //first line

        std::getline(fin,line); //Lx
        ss = std::stringstream(line);
        for (size_t i=0;i<Lx.rows();i++){
            ss >> Lx(i);
        }

        std::getline(fin,line); //Nx
        ss = std::stringstream(line);
        for (size_t i=0;i<Nx.rows();i++){
            ss >> Nx(i);
        }

        hx = Lx.array() / Nx.cast<double>().array();

        std::getline(fin,line); //periodic_props
        ss = std::stringstream(line);
        for (size_t i=0;i<periodic_props.rows();i++){
            ss >> periodic_props(i);
        }

        //print grid properties
        std::cout << "Grid properties (Lx = { ";
        for (size_t i=0;i<Lx.rows();i++){
            std::cout << Lx(i) << " ";
        }
        std::cout << "}, Nx = { ";
        for (size_t i=0;i<Nx.rows();i++){
            std::cout << Nx(i) << " ";
        }
        std::cout << "})." << std::endl;

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

void fixPosition(Job* job, Eigen::VectorXd& xIN){
    //wrap position vector around relevant axes

    //int tmp = floor(xIN[1]/hx[1]); //element i,j,k
    //double rem = xIN[1] - tmp*hx[1]; //position relative to element
    //int true_elem = tmp%Nx[1]; //adjusted element
    //double true_pos = true_elem*hx[1] + rem; //adjusted position


    int tmp;
    for (size_t pos=0; pos<periodic_props.size(); pos++) {
        if (periodic_props(pos) == 1) {
            tmp = floor(xIN[pos] / hx[pos]);
            xIN[pos] += (tmp % Nx[pos] * hx[pos]) - tmp * hx[pos];
            if (xIN[pos] < 0) {
                xIN[pos] += Lx[pos];
            }
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

int gridWhichElement(Job* job, Eigen::VectorXd xIN){
    fixPosition(job,xIN);

    bool inDomain;
    int elementID = floor(xIN[0]/hx[0]); //id in x-dimension
    for (size_t i=0;i<xIN.size();i++){
        inDomain = (xIN[i] < Lx[i] && xIN[i] >= 0);
        if (!inDomain) { //if xIn is outside domain, return -1
            return -1;
        }
        if (i == 1){
            //add number of elements in next layer of id in higher dimensions
            elementID += floor(xIN[i]/hx[i])*(Nx[i-1]);
        } else if (i == 2){
            elementID += floor(xIN[i]/hx[i])*(Nx[i-1])*(Nx[i-2]);
        }
    }
    return elementID;
}

/*----------------------------------------------------------------------------*/

bool gridInDomain(Job* job, Eigen::VectorXd xIN){
    fixPosition(job,xIN);

    for (size_t i=0;i<xIN.size();i++){
        if (!(xIN[i] <= Lx[i] && xIN[i] >= 0)) { //if xIn is outside domain, return -1
            return false;
        }
    }
    return true;
}

/*----------------------------------------------------------------------------*/

Eigen::VectorXd gridNodeIDToPosition(Job* job, int idIN){
    Eigen::VectorXi ijk = job->jobVector<int>();
    Eigen::VectorXd tmpVec = job->jobVector<double>();
    int tmp = idIN;
    //find i,j,k representation of node id
    for (size_t i=0;i<ijk.rows();i++){
        ijk(i) = tmp % (Nx(i)+1);
        tmp = tmp/(Nx(i)+1);
    }
    tmpVec = hx.array() * ijk.cast<double>().array();
    return tmpVec;
}

/*----------------------------------------------------------------------------*/

void gridEvaluateShapeFnValue(Job* job, Eigen::VectorXd xIN, std::vector<int>& nID, std::vector<double>& nVAL){
    fixPosition(job,xIN);

    Eigen::VectorXd rst = job->jobVector<double>();
    double tmp = 1.0;
    int elementID = gridWhichElement(job,xIN);
    if (elementID < 0){
        return;
    }

    //double test_sum = 0;
    for (size_t n=0;n<nodeIDs.cols();n++){
        if (nodeIDs(elementID,n) == -1){
            continue; //break out if node id is -1;
        }
        tmp = 1.0;

        //find local coordinates relative to nodal position
        //r = (x_p - x_n)/hx
        rst = xIN - x_n.row(nodeIDs(elementID,n)).transpose();
        for (size_t i=0;i<xIN.rows();i++){
            //check proximity to edge
            if (edge_n(nodeIDs(elementID,n),i) == 2) {
                tmp *= s(rst(i), hx(i));
            } else if (edge_n(nodeIDs(elementID,n),i) == 1) {
                tmp *= (s(rst(i), hx(i)) - s(rst(i)+2*hx(i), hx(i)));
            } else if (edge_n(nodeIDs(elementID,n),i) == -1) {
                tmp *= (s(rst(i), hx(i)) - s(rst(i)-2*hx(i), hx(i)));
            } else if (edge_n(nodeIDs(elementID,n),i) == 0) {
                tmp *= (s(rst(i), hx(i)) + 2.0*s(std::abs(rst(i))+hx(i), hx(i)));
            } else if (edge_n(nodeIDs(elementID,n),i) == 3) {
                tmp *= (s(rst(i), hx(i)) + s(rst(i) - Lx(i), hx(i))); //sum of two shape functions
            } else if (edge_n(nodeIDs(elementID,n),i) == -3) {
                tmp *= (s(rst(i), hx(i)) + s(rst(i) + Lx(i), hx(i))); //sum of two
            }
        }
        //if (tmp < 1e-10){
        //std::cout << "VERY SMALL!!! " << tmp << std::endl;
        //avoid small numbers for numerical stability
        //continue;
        //}

        //test_sum += tmp;
        nID.push_back(nntoni(nodeIDs(elementID,n)));
        nVAL.push_back(tmp);
    }
    //std::cout << elementID << " : " << test_sum << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void gridEvaluateShapeFnGradient(Job* job, Eigen::VectorXd xIN, std::vector<int>& nID, std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& nGRAD){
    fixPosition(job,xIN);

    Eigen::VectorXd rst = job->jobVector<double>();
    Eigen::VectorXd tmpVec = job->jobVector<double>(Job::ONES);
    double tmp = 1.0;
    int elementID = gridWhichElement(job,xIN);
    if (elementID < 0){
        return;
    }

    //Eigen::VectorXd test_grad = job->jobVector<double>(Job::ZERO);
    for (size_t n=0;n<nodeIDs.cols();n++){
        if (nodeIDs(elementID,n) == -1){
            continue; //break out if node id is -1;
        }

        //find local coordinates relative to nodal position
        //r = (x_p - x_n)
        rst = xIN - x_n.row(nodeIDs(elementID,n)).transpose();
        for (size_t i=0;i<xIN.rows();i++){
            //check proximity to edge
            if (edge_n(nodeIDs(elementID,n),i) == 2) {
                tmpVec(i) = g(rst(i), hx(i));
            } else if (edge_n(nodeIDs(elementID,n),i) == 1) {
                tmpVec(i) = (g(rst(i), hx(i)) - g(rst(i)+2*hx(i), hx(i)));
            } else if (edge_n(nodeIDs(elementID,n),i) == -1) {
                tmpVec(i) = (g(rst(i), hx(i)) - g(rst(i)-2*hx(i), hx(i)));
            } else if (edge_n(nodeIDs(elementID,n),i) == 0) {
                if (rst(i) >= 0) {
                    tmpVec(i) = (g(rst(i), hx(i)) + 2.0 * g(rst(i) + hx(i), hx(i)));
                } else if (rst(i) < 0){
                    tmpVec(i) = (g(rst(i), hx(i)) + 2.0 * g(rst(i) - hx(i), hx(i)));
                }
            } else if (edge_n(nodeIDs(elementID,n),i) == 3) {
                tmpVec(i) = (g(rst(i), hx(i)) + g(rst(i) - Lx(i), hx(i))); //sum of two shape functions
            } else if (edge_n(nodeIDs(elementID,n),i) == -3) {
                tmpVec(i) = (g(rst(i), hx(i)) + g(rst(i) + Lx(i), hx(i))); //sum of two
            }
        }

        for (size_t i=0;i<xIN.rows();i++){
            //check proximity to edge
            if (edge_n(nodeIDs(elementID,n),i) == 2) {
                tmpVec *= s(rst(i), hx(i));
                tmpVec(i) *= 1.0/s(rst(i), hx(i));
            } else if (edge_n(nodeIDs(elementID,n),i) == 1) {
                tmpVec *= (s(rst(i), hx(i)) - s(rst(i)+2*hx(i), hx(i)));
                tmpVec(i) *= 1.0/(s(rst(i), hx(i)) - s(rst(i)+2*hx(i), hx(i)));
            } else if (edge_n(nodeIDs(elementID,n),i) == -1) {
                tmpVec *= (s(rst(i), hx(i)) - s(rst(i)-2*hx(i), hx(i)));
                tmpVec(i) *= 1.0/(s(rst(i), hx(i)) - s(rst(i)-2*hx(i), hx(i)));
            } else if (edge_n(nodeIDs(elementID,n),i) == 0) {
                tmpVec *= (s(rst(i), hx(i)) + 2.0*s(std::abs(rst(i))+hx(i), hx(i)));
                tmpVec(i) *= 1.0/(s(rst(i), hx(i)) + 2.0*s(std::abs(rst(i))+hx(i), hx(i)));
            } else if (edge_n(nodeIDs(elementID,n),i) == 3) {
                tmpVec *= (s(rst(i), hx(i)) + s(rst(i) - Lx(i), hx(i))); //sum of two shape functions
                tmpVec(i) /= (s(rst(i), hx(i)) + s(rst(i) - Lx(i), hx(i))); //sum of two shape functions
            } else if (edge_n(nodeIDs(elementID,n),i) == -3) {
                tmpVec *= (s(rst(i), hx(i)) + s(rst(i) + Lx(i), hx(i))); //sum of two
                tmpVec(i) /= (s(rst(i), hx(i)) + s(rst(i) + Lx(i), hx(i))); //sum of two
            }
        }

        //test_grad = test_grad + tmpVec;
        nID.push_back(nntoni(nodeIDs(elementID,n)));
        nGRAD.push_back(tmpVec);
        tmpVec.setZero();
    }

    //std::cout << elementID << " : " << test_grad.transpose() << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

double gridNodalVolume(Job* job, int idIN){
    return v_n(idIN);
}

/*----------------------------------------------------------------------------*/

double gridElementVolume(Job* job, int idIN){
    return v_e;
}

/*----------------------------------------------------------------------------*/

double gridNodalTag(Job* job, int idIN){
    return -1;
}

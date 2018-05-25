//
// Created by aaron on 5/23/18.
// cartesian_cubic.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "grids.hpp"

/*----------------------------------------------------------------------------*/

double CartesianCubic::s(double x, double h){
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

double CartesianCubic::g(double x, double h){
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

/*----------------------------------------------------------------------------*/
//initialize grid assuming that job is set up
void CartesianCubic::init(Job* job){
    if (fp64_props.size() < job->DIM || int_props.size() < job->DIM){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i dimensions defined.\n",
                __FILE__, __func__, job->DIM);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        Lx = KinematicVector(job->JOB_TYPE);
        Nx = Eigen::VectorXi(job->DIM);
        hx = KinematicVector(job->JOB_TYPE);

        for (int pos=0;pos<hx.size();pos++){
            Lx(pos) = fp64_props[pos];
            Nx(pos) = int_props[pos];
            hx(pos) = Lx(pos) / Nx(pos);
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
//hidden initializer
void CartesianCubic::hiddenInit(Job* job){//called from loading or initializing functions
    //needs to have Lx,Nx,hx defined

    //count nodes
    node_count = 1;
    element_count = 1;
    npe = 1;
    for (size_t i=0;i<Nx.rows();i++){
        node_count *= (Nx(i)+1);
        element_count *= Nx(i);
        npe *= 4; //4 linear node contributions per element
    }

    x_n = KinematicVectorArray(node_count,job->JOB_TYPE);
    for (size_t i=0;i<x_n.size();i++){
        x_n(i) = nodeIDToPosition(job,i);
    }

    //initialize A matrix
    //for mapping position in cube to id
    //0 -> -1,-1,-1
    //1 -> +0,-1,-1
    //...
    //64 -> +2,+2,+2
    Eigen::VectorXi i_rel = Eigen::VectorXi(job->DIM);
    A = Eigen::MatrixXi(npe,job->DIM);
    i_rel.setOnes();
    i_rel = -1*i_rel;
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
    Eigen::VectorXi ijk = Eigen::VectorXi(job->DIM);
    int tmp;
    nodeIDs.resize(element_count,npe);
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
                if ((ijk(i) + A(n,i) < 0) || (ijk(i) + A(n,i) > Nx(i))) {
                    nodeIDs(e, n) = -1;
                    break;
                }

                if (i==0) {
                    nodeIDs(e, n) += (ijk(i) + A(n,i));
                } else if (i==1) {
                    nodeIDs(e, n) += (ijk(i) + A(n,i)) * (Nx(i-1)+1);
                } else if (i==2) {
                    nodeIDs(e, n) += (ijk(i) + A(n,i)) * (Nx(i-1)+1) * (Nx(i-2)+1);
                }
            }
        }
    }

    //element volume
    v_e = 1;
    for (size_t pos=0;pos<hx.rows();pos++){
        v_e *= hx(pos);
    }

    //nodal volume and edges
    edge_n = KinematicVectorArray(x_n.size(), job->JOB_TYPE);
    v_n.resize(x_n.size());
    v_n.setZero();
    for (size_t n=0;n<x_n.size();n++){
        tmp = n;
        v_n(n) = 1;
        for (size_t i=0;i<ijk.rows();i++){
            ijk(i) = tmp % (Nx(i)+1);
            tmp = tmp / (Nx(i)+1);

            if (ijk(i) > 1 && ijk(i) < (Nx(i)-1)){
                v_n(n) *= hx(i);
                edge_n(n,i) = 2; //interior
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
//
void CartesianCubic::writeFrame(Job* job, Serializer* serializer){
    serializer->writeVectorArray(edge_n,"edge_id");
    return;
}


/*----------------------------------------------------------------------------*/
//
std::string CartesianCubic::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}


/*----------------------------------------------------------------------------*/
//
int CartesianCubic::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}



/*----------------------------------------------------------------------------*/
//
int CartesianCubic::whichElement(Job* job, KinematicVector& xIN){
    return CartesianLinear::cartesianWhichElement(job, xIN, Lx, hx, Nx);
}

/*----------------------------------------------------------------------------*/
//
bool CartesianCubic::inDomain(Job* job, KinematicVector& xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "inDomain failed");
    for (int i=0;i<xIN.size();i++){
        if (!(xIN[i] <= Lx[i] && xIN[i] >= 0)) { //if xIn is outside domain, return -1
            return false;
        }
    }
    return true;
}


/*----------------------------------------------------------------------------*/
//
KinematicVector CartesianCubic::nodeIDToPosition(Job* job, int idIN){
    Eigen::VectorXi ijk(job->DIM);
    KinematicVector tmpVec(job->JOB_TYPE);
    int tmp = idIN;
    //find i,j,k representation of node id
    for (int i=0;i<ijk.rows();i++){
        ijk(i) = tmp % (Nx(i)+1);
        tmp = tmp/(Nx(i)+1);
    }
    for (int i=0;i<tmpVec.DIM;i++){
        tmpVec[i] = hx(i)*ijk(i);
    }
    return tmpVec;
}



/*----------------------------------------------------------------------------*/
//
void CartesianCubic::evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "evaluateShapeFnValue failed");
    KinematicVector rst(job->JOB_TYPE);
    double tmp = 1.0;
    int elementID = whichElement(job,xIN);
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
        rst = xIN - x_n(nodeIDs(elementID,n));
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
            }
        }
        //if (tmp < 1e-10){
        //std::cout << "VERY SMALL!!! " << tmp << std::endl;
        //avoid small numbers for numerical stability
        //continue;
        //}

        //test_sum += tmp;
        nID.push_back(nodeIDs(elementID,n));
        nVAL.push_back(tmp);
    }
    //std::cout << elementID << " : " << test_sum << std::endl;
    return;
}


/*----------------------------------------------------------------------------*/
//
void CartesianCubic::evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && nGRAD.VECTOR_TYPE == Lx.VECTOR_TYPE && "evaluateShapeFnGradient failed");
    KinematicVector rst(job->JOB_TYPE);
    KinematicVector tmpVec(job->JOB_TYPE);
    double tmp = 1.0;
    int elementID = whichElement(job,xIN);
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
        rst = xIN - x_n(nodeIDs(elementID,n));
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
            }
        }

        //test_grad = test_grad + tmpVec;
        nID.push_back(nodeIDs(elementID,n));
        nGRAD.push_back(tmpVec);
        tmpVec.setZero();
    }

    //std::cout << elementID << " : " << test_grad.transpose() << std::endl;
    return;
}


/*----------------------------------------------------------------------------*/
//
double CartesianCubic::nodeVolume(Job* job, int idIN){
    return v_n(idIN);
}


/*----------------------------------------------------------------------------*/
//
double CartesianCubic::elementVolume(Job* job, int idIN){
    return v_e;
}


/*----------------------------------------------------------------------------*/
//
int CartesianCubic::nodeTag(Job* job, int idIN){
    return -1;
}

/*----------------------------------------------------------------------------*/
//
void CartesianCubic::writeHeader(Job *job, Body *body, Serializer *serializer, std::ofstream &nfile, int SPEC) {
    int nlen = body->nodes->x.size();

    nfile << "ASCII\n";
    nfile << "DATASET UNSTRUCTURED_GRID\n";

    nfile << "POINTS " << nlen << " double\n";
    for (int i=0;i<nlen;i++){
        //vtk files require x,y,z
        for (int pos = 0; pos < 3; pos++){
            if (pos < body->nodes->x.DIM && (body->nodes->active(i) != 0) && std::isfinite(body->nodes->x(i,pos))){
                nfile << body->nodes->x(i,pos) << " ";
            } else {
                nfile << "0 ";
            }
        }
        nfile << "\n";
    }

    if (job->DIM == 1) {
        //use lines
        nfile << "CELLS " << element_count << " " << 3 * element_count << "\n";
        for (int e = 0; e < element_count; e++) {
            nfile << "2 " << nodeIDs(e,0) << " " << nodeIDs(e,1) << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e = 0; e < element_count; e++) {
            nfile << "3\n";
        }
    } else if (job->DIM == 2){
        nfile << "CELLS " << element_count << " " << 5 * element_count << "\n";
        for (int e = 0; e < element_count; e++){
            nfile << "4 " << nodeIDs(e,0) << " " << nodeIDs(e,1) << " " << nodeIDs(e,2) << " " << nodeIDs(e,3) << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e=0; e<element_count; e++){
            nfile << "8\n";
        }
    } else if (job->DIM == 3){
        nfile << "CELLS " << element_count << " " << 9 * element_count << "\n";
        for (int e = 0; e < element_count; e++){
            nfile << "8 ";
            nfile << nodeIDs(e, 0) << " " << nodeIDs(e, 1) << " ";
            nfile << nodeIDs(e, 2) << " " << nodeIDs(e, 3) << " ";
            nfile << nodeIDs(e, 4) << " " << nodeIDs(e, 5) << " ";
            nfile << nodeIDs(e, 6) << " " << nodeIDs(e, 7) << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e=0; e<element_count; e++){
            nfile << "11\n";
        }
    }

    nfile << "POINT_DATA " << nlen << "\n";
}
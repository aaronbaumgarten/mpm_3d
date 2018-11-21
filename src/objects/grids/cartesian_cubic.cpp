//
// Created by aaron on 5/23/18.
// cartesian_cubic.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <iomanip>
#include <fstream>

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
    if (fp64_props.size() < GRID_DIM || int_props.size() < GRID_DIM){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i dimensions defined.\n",
                __FILE__, __func__, GRID_DIM);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        Lx = KinematicVector(job->JOB_TYPE);
        Nx = Eigen::VectorXi(job->DIM);
        hx = KinematicVector(job->JOB_TYPE);

        for (int pos=0;pos<GRID_DIM;pos++){
            Lx(pos) = fp64_props[pos];
            Nx(pos) = int_props[pos];
            hx(pos) = Lx(pos) / Nx(pos);
        }

        for (int pos=GRID_DIM;pos<hx.size();pos++){
            Lx(pos) = 0;
            Nx(pos) = 0;
            hx(pos) = 0;
        }

        //print grid properties
        std::cout << "Grid properties (Lx = { ";
        for (int i=0;i<Lx.size();i++){
            std::cout << Lx(i) << " ";
        }
        std::cout << "}, Nx = { ";
        for (int i=0;i<Nx.size();i++){
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
    for (int i=0;i<GRID_DIM;i++){
        node_count *= (Nx(i)+1);
        element_count *= Nx(i);
        npe *= 4; //4 linear node contributions per element
    }

    x_n = KinematicVectorArray(node_count,job->JOB_TYPE);
    for (int i=0;i<x_n.size();i++){
        x_n(i) = nodeIDToPosition(job,i);
    }

    //initialize A matrix
    //for mapping position in cube to id
    //0 -> -1,-1,-1
    //1 -> +0,-1,-1
    //...
    //64 -> +2,+2,+2
    Eigen::VectorXi i_rel = Eigen::VectorXi(GRID_DIM);
    A = Eigen::MatrixXi(npe,GRID_DIM);
    i_rel.setOnes();
    i_rel = -1*i_rel;
    for (int n=0; n<npe;n++){
        for (int i=0;i<i_rel.rows();i++){
            A(n,i) = i_rel(i);
        }
        for (int i=0;i<i_rel.rows();i++) {
            if (i_rel(i) < 2){
                i_rel(i) += 1;
                break;
            } else {
                i_rel(i) = -1;
            }
        }
    }

    //setup element to node map
    Eigen::VectorXi ijk = Eigen::VectorXi(GRID_DIM);
    int tmp;
    nodeIDs.resize(element_count,npe);
    nodeIDs.setZero();
    for (int e=0;e<nodeIDs.rows();e++){
        tmp = e;
        //find i,j,k count for element position
        for (int i=0;i<ijk.rows();i++){
            ijk(i) = tmp % Nx(i);
            tmp = tmp / Nx(i);
        }

        //find node ids for element
        for (int n=0;n<nodeIDs.cols();n++){
            for (int i=0;i<ijk.rows();i++) {
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
    for (int pos=0;pos<GRID_DIM;pos++){
        v_e *= hx(pos);
    }

    //nodal volume and edges
    edge_n = KinematicVectorArray(x_n.size(), job->JOB_TYPE);
    v_n.resize(x_n.size());
    v_n.setZero();
    for (int n=0;n<x_n.size();n++){
        tmp = n;
        v_n(n) = 1;
        for (int i=0;i<ijk.rows();i++){
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

        for (int i=ijk.rows(); i<edge_n.DIM; i++){
            edge_n(n,i) = 2; //interior (but not really important)
        }
    }

    //node mapping matrix
    S_grid = MPMScalarSparseMatrix(node_count, node_count);
    std::vector<int> nvec(0);
    std::vector<double> valvec(0);
    for (int i=0; i<node_count; i++) {

        KinematicVector x_i = x_n(i);
        for (int pos = 0; pos < GRID_DIM; pos++){
            if (x_i(pos) == Lx(pos)){
                //adjust point just a little bit for post-processing
                x_i(pos) = Lx(pos) - hx(pos)/100.0;
            }
        }

        nvec.resize(0);
        valvec.resize(0);
        evaluateBasisFnValue(job, x_i, nvec, valvec);
        for (int j = 0; j < nvec.size(); j++) {
            S_grid.push_back(nvec[j], i, valvec[j]); //node, point, value
        }
    }

    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    //nodal surface integral
    s_n.resize(x_n.size());
    s_n.setZero();
    for (int n=0;n<x_n.size();n++){
        tmp = n;
        for (int i=0;i<ijk.rows();i++){
            ijk(i) = tmp % (Nx(i)+1);
            tmp = tmp/(Nx(i)+1);

            if (ijk(i) == 0 || ijk(i) == Nx(i)){
                s_n(n) = 1;
                for (int pos=0;pos<(GRID_DIM-1);pos++){
                    s_n(n) *= hx(pos);
                }

                //approximate adjustment for axisymetric case
                if (job->JOB_TYPE == job->JOB_AXISYM){
                    tmpVec = nodeIDToPosition(job, n);
                    s_n(n) *= tmpVec(0);
                }
                break;
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
    return CartesianLinear::cartesianWhichElement(job, xIN, Lx, hx, Nx, GRID_DIM);
}

/*----------------------------------------------------------------------------*/
//
bool CartesianCubic::inDomain(Job* job, KinematicVector& xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "inDomain failed");
    for (int i=0;i<GRID_DIM;i++){
        if (!(xIN[i] <= Lx[i] && xIN[i] >= 0)) { //if xIn is outside domain, return -1
            return false;
        }
    }
    return true;
}


/*----------------------------------------------------------------------------*/
//
KinematicVector CartesianCubic::nodeIDToPosition(Job* job, int idIN){
    Eigen::VectorXi ijk(GRID_DIM);
    KinematicVector tmpVec(job->JOB_TYPE);
    int tmp = idIN;
    //find i,j,k representation of node id
    for (int i=0;i<ijk.rows();i++){
        ijk(i) = tmp % (Nx(i)+1);
        tmp = tmp/(Nx(i)+1);
    }
    for (int i=0;i<GRID_DIM;i++){
        tmpVec[i] = hx(i)*ijk(i);
    }
    for (int i=GRID_DIM;i<hx.size();i++){
        tmpVec[i] = 0;
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
    for (int n=0;n<nodeIDs.cols();n++){
        if (nodeIDs(elementID,n) == -1){
            continue; //break out if node id is -1;
        }
        tmp = 1.0;

        //find local coordinates relative to nodal position
        //r = (x_p - x_n)/hx
        rst = xIN - x_n(nodeIDs(elementID,n));
        for (int i=0;i<GRID_DIM;i++){
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
    for (int n=0;n<nodeIDs.cols();n++){
        if (nodeIDs(elementID,n) == -1){
            continue; //break out if node id is -1;
        }

        //find local coordinates relative to nodal position
        //r = (x_p - x_n)
        rst = xIN - x_n(nodeIDs(elementID,n));
        for (int i=0;i<GRID_DIM;i++){
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

        for (int i=0;i<GRID_DIM;i++){
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

        for (int i=GRID_DIM;i<tmpVec.size();i++){
            tmpVec(i) = 0;
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
double CartesianCubic::nodeSurfaceArea(Job *job, int idIN) {
    return s_n(idIN);
}

/*----------------------------------------------------------------------------*/
//
void CartesianCubic::writeHeader(Job *job, Body *body, Serializer *serializer, std::ofstream &nfile, int SPEC) {
    //filter the velocity field only for now
    KinematicVectorArray grid_velocity = S_grid.operate(body->nodes->x_t, MPMSparseMatrixBase::TRANSPOSED);

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
            if (pos < body->nodes->x.DIM && (body->nodes->active(i) != 0) && std::isfinite(body->nodes->x(i,pos))){
                nfile << body->nodes->x(i,pos) << " ";
            } else {
                nfile << "0 ";
            }
        }
        nfile << "\n";
    }

    if (GRID_DIM == 1) {
        //use lines
        nfile << "CELLS " << element_count << " " << 3 * element_count << "\n";
        for (int e = 0; e < element_count; e++) {
            nfile << "2 " << nodeIDs(e,1) << " " << nodeIDs(e,2) << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e = 0; e < element_count; e++) {
            nfile << "3\n";
        }
    } else if (GRID_DIM == 2){
        nfile << "CELLS " << element_count << " " << 5 * element_count << "\n";
        for (int e = 0; e < element_count; e++){
            nfile << "4 " << nodeIDs(e,5) << " " << nodeIDs(e,6) << " " << nodeIDs(e,9) << " " << nodeIDs(e,10) << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e=0; e<element_count; e++){
            nfile << "8\n";
        }
    } else if (GRID_DIM == 3){
        nfile << "CELLS " << element_count << " " << 9 * element_count << "\n";
        for (int e = 0; e < element_count; e++){
            nfile << "8 ";
            nfile << nodeIDs(e, 21) << " " << nodeIDs(e, 22) << " ";
            nfile << nodeIDs(e, 25) << " " << nodeIDs(e, 26) << " ";
            nfile << nodeIDs(e, 37) << " " << nodeIDs(e, 38) << " ";
            nfile << nodeIDs(e, 41) << " " << nodeIDs(e, 42) << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e=0; e<element_count; e++){
            nfile << "11\n";
        }
    }

    nfile << "POINT_DATA " << nlen << "\n";

    //write out grid_velocity data
    nfile << "VECTORS grid_velocity double\n";
    for (int i = 0; i < nlen; i++){
        //vtk format requires x,y,z
        for (int pos = 0; pos < 3; pos++){
            if ((body->nodes->active(i) == 1) && std::isfinite(grid_velocity(i,pos))){
                nfile << grid_velocity(i,pos) << " ";
            } else {
                nfile << "0 ";
            }
        }
        nfile << "\n";
    }
}
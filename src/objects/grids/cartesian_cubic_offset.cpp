//
// Created by aaron on 10/11/18.
// cartesian_cubic_offset.cpp
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
//initialize grid assuming that job is set up
void CartesianCubic_Offset::init(Job* job){
    if (job->JOB_TYPE != job->JOB_AXISYM){
        std::cerr << "ERROR: CartesianCubic_Offset requires JOB_TYPE = JOB_AXISYM. Exiting." << std::endl;
    }

    if (fp64_props.size() < GRID_DIM+1 || int_props.size() < GRID_DIM){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i dimensions and 1 offset defined.\n",
                __FILE__, __func__, GRID_DIM);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        Lx = KinematicVector(job->JOB_TYPE);
        Nx = Eigen::VectorXi(job->DIM);
        hx = KinematicVector(job->JOB_TYPE);
        x_offset = fp64_props[GRID_DIM];

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
        std::cout << "}, x_offset = " << x_offset << ", Nx = { ";
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
void CartesianCubic_Offset::hiddenInit(Job* job){//called from loading or initializing functions
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
int CartesianCubic_Offset::whichElement(Job* job, KinematicVector& xIN){
    //adjust input position for regular cartesian mapping
    //xIN[0] -= x_offset;
    KinematicVector tmpVec = xIN;
    tmpVec[0] -= x_offset;
    
    return CartesianLinear::cartesianWhichElement(job, tmpVec, Lx, hx, Nx, GRID_DIM);
}

/*----------------------------------------------------------------------------*/
//
bool CartesianCubic_Offset::inDomain(Job* job, KinematicVector& xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "inDomain failed");

    if (!(xIN[0] <= Lx[0]+x_offset && xIN[0] >= x_offset)) { //if xIn is outside domain, return -1
        return false;
    }

    for (int i=1;i<GRID_DIM;i++){
        if (!(xIN[i] <= Lx[i] && xIN[i] >= 0)) { //if xIn is outside domain, return -1
            return false;
        }
    }
    return true;
}


/*----------------------------------------------------------------------------*/
//
KinematicVector CartesianCubic_Offset::nodeIDToPosition(Job* job, int idIN){
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

    //adjust x-position
    tmpVec[0] += x_offset;

    return tmpVec;
}


/*----------------------------------------------------------------------------*/
//
std::string CartesianCubic_Offset::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}


/*----------------------------------------------------------------------------*/
//
int CartesianCubic_Offset::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

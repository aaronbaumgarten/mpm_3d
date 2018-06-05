//
// Created by aaron on 5/25/18.
// cartesian_cubic_custom.cpp
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
//
void CartesianCubicCustom::init(Job* job){
    if (fp64_props.size() < job->DIM || int_props.size() < 2*job->DIM){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i dimensions and %i properties defined.\n",
                __FILE__, __func__, job->DIM, job->DIM);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        //store length, number of linear nodes, and deltas
        Lx = KinematicVector(job->JOB_TYPE);
        Nx = Eigen::VectorXi(job->DIM);
        hx = KinematicVector(job->JOB_TYPE);

        for (int pos=0;pos<hx.size();pos++){
            Lx(pos) = fp64_props[pos];
            Nx(pos) = int_props[pos];
            hx(pos) = Lx(pos) / Nx(pos);
        }

        periodic_props = Eigen::VectorXi(job->DIM);
        for (int pos=0;pos<job->DIM;pos++){
            periodic_props(pos) = int_props[job->DIM + pos];
        }
        //print grid properties
        std::cout << "Grid properties (Lx = { ";
        for (int i=0;i<Lx.rows();i++){
            std::cout << Lx(i) << " ";
        }
        std::cout << "}, Nx = { ";
        for (int i=0;i<Nx.rows();i++){
            std::cout << Nx(i) << " ";
        }
        std::cout << "})." << std::endl;
    }

    //call initializing function
    hiddenInit(job);

    std::cout << "Grid Initialized." << std::endl;

    return;
}

void CartesianCubicCustom::hiddenInit(Job* job){
    //called from loading or initializing functions
    //needs to have Lx,Nx,hx defined

    //count nodes
    node_count = 1;
    element_count = 1;
    npe = 1;
    for (int i=0;i<Nx.rows();i++){
        node_count *= (Nx(i)+1);
        element_count *= Nx(i);
        npe *= 4; //4 linear node contributions per element
    }

    x_n = KinematicVectorArray(node_count,job->JOB_TYPE);
    for (int i=0;i<x_n.size();i++){
        x_n(i) = (nodeIDToPosition(job,i));
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
    Eigen::VectorXi ijk = Eigen::VectorXi(job->DIM);
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
    nntoni.resize(node_count);
    for (int i=0; i<node_count; i++){
        tmp = i;
        //find i,j,k for node position
        for (int i=0;i<ijk.rows();i++){
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
        for (int pos=0;pos<ijk.rows();pos++){
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
    for (int pos=0;pos<hx.rows();pos++){
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
//
std::string CartesianCubicCustom::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

int CartesianCubicCustom::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//
void CartesianCubicCustom::fixPosition(Job* job, KinematicVector& xIN){
    //wrap position vector around relevant axes
    int tmp;
    for (int pos=0; pos<periodic_props.size(); pos++) {
        if (periodic_props(pos) == PERIODIC) {
            tmp = floor(xIN[pos] / hx[pos]);
            xIN[pos] += (tmp % Nx[pos] * hx[pos]) - tmp * hx[pos];
            if (xIN[pos] < 0) {
                xIN[pos] += Lx[pos];
            }
        }
    }
    return;
}

int CartesianCubicCustom::whichElement(Job* job, KinematicVector& xIN){
    KinematicVector tmpX = xIN;
    fixPosition(job, tmpX);
    //return standard cartesian which element
    return CartesianLinear::cartesianWhichElement(job, tmpX, Lx, hx, Nx);
}

bool CartesianCubicCustom::inDomain(Job* job, KinematicVector& xIN){
    //adjust xIN
    KinematicVector tmpX = xIN;
    fixPosition(job,tmpX);

    for (int i=0;i<tmpX.DIM;i++){
        if (!(tmpX[i] <= Lx[i] && tmpX[i] >= 0)) { //if xIn is outside domain, return -1
            return false;
        }
    }
    return true;
}


/*----------------------------------------------------------------------------*/
//
void CartesianCubicCustom::evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL){
    KinematicVector tmpX = xIN;
    fixPosition(job,tmpX);

    KinematicVector rst = KinematicVector(job->JOB_TYPE);
    double tmp = 1.0;
    int elementID = whichElement(job,tmpX);
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
        rst = tmpX - x_n(nodeIDs(elementID,n));
        for (int i=0;i<tmpX.rows();i++){
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

void CartesianCubicCustom::evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD){
    KinematicVector tmpX = xIN;
    fixPosition(job,tmpX);

    KinematicVector rst = KinematicVector(job->JOB_TYPE);
    KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
    tmpVec.setOnes();

    double tmp = 1.0;
    int elementID = whichElement(job,tmpX);
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
        rst = tmpX - x_n(nodeIDs(elementID,n));
        for (int i=0;i<tmpX.rows();i++){
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

        for (int i=0;i<tmpX.rows();i++){
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
//
// Created by aaron on 5/31/18.
// regular_2D_taylor_couetta.cpp
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
double Regular2DTaylorCouetteCell::s_linear(double x, double h){
    if (x < -h){
        return 0;
    } else if (x < 0){
        return (1.0 + x/h);
    } else if (x < h){
        return (1.0 - x/h);
    } else {
        return 0;
    }
}

double Regular2DTaylorCouetteCell::g_linear(double x, double h){
    if (x < -h){
        return 0;
    } else if (x < 0){
        return (1.0/h);
    } else if (x < h){
        return (-1.0/h);
    } else {
        return 0;
    }
}

double Regular2DTaylorCouetteCell::s_cubic(double x, double h){
    return CartesianCubic::s(x,h);
}

double Regular2DTaylorCouetteCell::g_cubic(double x, double h){
    return CartesianCubic::g(x,h);
}

/*----------------------------------------------------------------------------*/
//map a point in cartesian space to radial coordinates
KinematicVector Regular2DTaylorCouetteCell::cPoint_to_rPoint(const KinematicVector &xyz){
    KinematicVector tmp = xyz; //copy vector (z value remains unchanged)
    tmp(0) = std::sqrt(xyz(0)*xyz(0) + xyz(1)*xyz(1));
    if (xyz(0) == 0 && xyz(1) >= 0){
        tmp(1) = pi/2;
    } else if (xyz(0) == 0 && xyz(1) < 0) {
        tmp(1) = -pi / 2;
    } else if (xyz(0) > 0){
        tmp(1) = std::atan(xyz(1)/xyz(0));
    } else if (xyz(0) < 0){
        tmp(1) = pi + std::atan(xyz(1)/xyz(0));
    }

    //theta between 0 and 2*pi
    if (tmp(1) < 0){
        tmp(1) += 2*pi;
    }

    return tmp;
}

//map radial coordinates to a point in cartesian space
KinematicVector Regular2DTaylorCouetteCell::rPoint_to_cPoint(const KinematicVector &rtz){
    KinematicVector tmp = rtz; //copy vector (z value remains unchanged)
    tmp(0) = rtz(0)*std::cos(rtz(1));
    tmp(1) = rtz(0)*std::sin(rtz(1));
    return tmp;
}

//map vector at given position in cartesian space to a vector in r, theta space
KinematicVector Regular2DTaylorCouetteCell::cVector_to_rVector(const KinematicVector &vec, const KinematicVector &xyz){
    KinematicVector rtz = cPoint_to_rPoint(xyz); //change to radial coordinates
    KinematicVector tmp = vec; //copy vector (z value remains unchanged)
    KinematicVector er(xyz.VECTOR_TYPE), et(xyz.VECTOR_TYPE); //unit directions in cartesian space

    er(0) = std::cos(rtz(1)); er(1) = std::sin(rtz(1));
    et(0) = -std::sin(rtz(1)); et(1) = std::cos(rtz(1));

    tmp(0) = vec.dot(er);
    tmp(1) = vec.dot(et);

    return tmp;
}

//map vector at given radial coordinate and return vector in cartesian space
KinematicVector Regular2DTaylorCouetteCell::rVector_to_cVector(const KinematicVector &vec, const KinematicVector &rtz){
    KinematicVector tmp = vec; //copy vector(z value remains unchanged)
    tmp(0) = vec(0)*std::cos(rtz(1)) - vec(1)*std::sin(rtz(1));
    tmp(1) = vec(0)*std::sin(rtz(1)) + vec(1)*std::cos(rtz(1));
    return tmp;
}

/*----------------------------------------------------------------------------*/
//
void Regular2DTaylorCouetteCell::init(Job* job){
    if (job->JOB_TYPE != Job::JOB_2D){
        std::cerr << "ERROR: Regual2DTaylorCouetteCell requires JOB_TYPE = JOB_2D, got: " << job->JOB_TYPE << ". Exiting." << std::endl;
        exit(0);
    }
    if (fp64_props.size() < 2 || int_props.size() < 2){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i dimensions defined.\n",
                __FILE__, __func__, 2);
        exit(0);
    } else {
        //store length, number of linear nodes, and deltas
        Lx = KinematicVector(job->JOB_TYPE);
        Nx = Eigen::VectorXi(job->DIM);
        hx = KinematicVector(job->JOB_TYPE);

        Ri = fp64_props[0];
        Ro = fp64_props[1];
        Lx(0) = Ro - Ri;    //radial dimension
        Lx(1) = 2*pi;       //theta dimension

        for (int pos=0;pos<hx.size();pos++){
            Nx(pos) = int_props[pos];
            hx(pos) = Lx(pos) / Nx(pos);
        }

        if (int_props.size() >= 3){
            order = int_props[2];
        }

        //print grid properties
        std::cout << "Grid properties (Ri = " << Ri << ", Ro = " << Ro << ", Nx = { ";
        for (int i=0;i<Nx.size();i++){
            std::cout << Nx(i) << " ";
        }
        std::cout << "}, order = " << order << ")." << std::endl;
    }

    //call initializing function
    if (order == LINEAR){
        linearInit(job);
    } else if (order == CUBIC){
        cubicInit(job);
    } else {
        std::cerr << "ERROR: Regular2DTaylorCouetteCell does not have order " << order << " basis fn's implemented. Exiting." << std::endl;
        exit(0);
    }

    std::cout << "Grid Initialized." << std::endl;

    return;
}

void Regular2DTaylorCouetteCell::linearInit(Job *job){
    //called from loading or initializing functions
    //needs to have Lx,Nx,hx defined

    //count nodes (we know this is only 2D)
    //theta-direction nodes are wrapped around
    node_count = (Nx(0)+1)*Nx(1);
    element_count = Nx(0)*Nx(1);
    npe = 4;

    x_n = KinematicVectorArray(node_count,job->JOB_TYPE);
    for (int i=0;i<x_n.size();i++){
        x_n(i) = nodeIDToPosition(job,i);
    }

    //initialize A matrix
    //for mapping position in cube to id
    //0 -> +0,+0,+0
    //1 -> +1,+0,+0
    //...
    //8 -> +1,+1,+1
    Eigen::VectorXi onoff(job->DIM);
    A = Eigen::MatrixXi(npe,job->DIM);
    onoff.setZero();
    for (int n=0; n<npe;n++){
        for (int i=0;i<onoff.rows();i++){
            A(n,i) = onoff(i);
        }
        for (int i=0;i<onoff.rows();i++) {
            if (onoff(i) == 0){
                onoff(i) = 1;
                break;
            } else {
                onoff(i) = 0;
            }
        }
    }

    //setup element to node map
    /*
     * r - direction
     * 0 --- 1 --- 2 --- ... N-1 --- N --- N+1
     *
     * theta - direction
     * /- 0 --- 1 --- 2 --- ... N-1 --- N -/
     */
    Eigen::VectorXi ijk(job->DIM);
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
                if (i==0) {
                    //r-direction
                    nodeIDs(e, n) += (ijk(i) + A(n,i));
                } else if (i==1) {
                    //theta-direction
                    if (ijk(i) + A(n,i) < Nx(i)) {
                        nodeIDs(e, n) += (ijk(i) + A(n, i)) * (Nx(i - 1) + 1);
                    } else {
                        //do nothing, same as j = 0
                    }
                }
                /*
                else if (i==2) {
                    nodeIDs(e, n) += (ijk(i) + A(n,i)) * (Nx(i-1)+1) * (Nx(i-2)+1);
                }
                */
            }
        }
    }

    //element volume
    double ro, ri;
    KinematicVector rtz;
    v_e.resize(element_count);
    v_e.setOnes();
    for (int e=0;e<element_count;e++){
        //volume = pi*(ro^2 - ri^2) / (number of theta slices)
        rtz = cPoint_to_rPoint(x_n(nodeIDs(e,0)));
        ri = rtz(0);
        rtz = cPoint_to_rPoint(x_n(nodeIDs(e,1)));
        ro = rtz(0);
        v_e(e) = pi*(ro*ro - ri*ri) / Nx(1);
    }

    //nodal volume
    v_n.resize(x_n.size());
    v_n.setZero();

    //use 3-pt Gaussian Quadrature to determine nodal value (up to 4th order bases)
    double x0; //element minimum radius
    double xtmp, rtmp; //gauss point in radial coords
    double x1 = 0;
    double x2 = std::sqrt(3.0/5.0);
    double x3 = -x2;
    double w1 = 8.0/9.0;
    double w2 = 5.0/9.0;
    double w3 = w2;
    KinematicVector diff_x_rtz;

    for (int e=0;e<nodeIDs.rows();e++){
        //x0 = x_n(nodeIDs(e,0))[0];
        for (int c=0;c<nodeIDs.cols();c++) {
            rtz = cPoint_to_rPoint(x_n(nodeIDs(e,0)));
            diff_x_rtz = cPoint_to_rPoint(x_n(nodeIDs(e,0))) - cPoint_to_rPoint(x_n(nodeIDs(e,c)));

            rtmp = rtz[0] + hx(0)/2 * (x1 + 1);
            xtmp = diff_x_rtz[0] + hx(0)/2 * (x1 + 1);
            v_n(nodeIDs(e, c)) += hx(0)/2 * w1 * rtmp * s_linear(xtmp, hx(0));

            rtmp = rtz[0] + hx(0)/2 * (x2 + 1);
            xtmp = diff_x_rtz[0] + hx(0)/2 * (x2 + 1);
            v_n(nodeIDs(e, c)) += hx(0)/2 * w2 * rtmp * s_linear(xtmp, hx(0));

            rtmp = rtz[0] + hx(0)/2 * (x3 + 1);
            xtmp = diff_x_rtz[0] + hx(0)/2 * (x3 + 1);
            v_n(nodeIDs(e, c)) += hx(0)/2 * w3 * rtmp * s_linear(xtmp, hx(0));
        }
    }

    //scale v_n by theta contribution
    //each integral will be repeated twice
    v_n *= hx(1)/2.0;

    //edges
    node_tag = Eigen::VectorXi(x_n.size());
    node_tag.setZero();
    for (int n=0;n<x_n.size();n++){
        tmp = n;
        for (int i=0;i<1;i++){
            ijk(i) = tmp % (Nx(i)+1);
            tmp = tmp / (Nx(i)+1);

            if (ijk(i) > 0 && ijk(i) < Nx(i)){
                //node_tag(n) = 0; //interior
            } else {
                if (ijk(i) == 0){
                    node_tag(n) = -1;
                } else {
                    node_tag(n) = 1; //boundary
                }
            }
        }
    }

    return;
}

void Regular2DTaylorCouetteCell::cubicInit(Job *job) {//called from loading or initializing functions
    //needs to have Lx,Nx,hx defined

    //count nodes (we know this is only 2D)
    //theta-direction nodes are wrapped around
    node_count = (Nx(0)+1)*Nx(1);
    element_count = Nx(0)*Nx(1);
    npe = 16;

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
                if (i==0) {
                    if (ijk(i) + A(n,i) < 0 || ijk(i) + A(n,i) >= (Nx(i)+1)){
                        nodeIDs(e,n) = -1;
                        break;
                    }
                    //r-direction
                    nodeIDs(e, n) += (ijk(i) + A(n,i));
                } else if (i==1) {
                    //theta-direction
                    if (ijk(i) + A(n,i) == -1) {
                        nodeIDs(e, n) += (Nx(1) - 1) * (Nx(i-1) + 1);
                    } else if (ijk(i) + A(n,i) < Nx(i)) {
                        nodeIDs(e, n) += (ijk(i) + A(n, i)) * (Nx(i - 1) + 1);
                    } else if (ijk(i) + A(n,i) == Nx(i)){
                        //do nothing, same as j = 0
                    } else if (ijk(i) + A(n,i) == Nx(i) + 1) {
                        nodeIDs(e, n) += Nx(i - 1) + 1;
                    }
                }
            }
        }
    }

    //element volume
    //element volume
    double ro, ri;
    KinematicVector rtz;
    v_e.resize(x_n.size());
    v_e.setOnes();
    for (int e=0;e<element_count;e++){
        //volume = pi*(ro^2 - ri^2) / (number of theta slices)
        rtz = cPoint_to_rPoint(x_n(nodeIDs(e,0)));
        ri = rtz(0);
        rtz = cPoint_to_rPoint(x_n(nodeIDs(e,1)));
        ro = rtz(0);
        v_e(e) = pi*(ro*ro - ri*ri) / Nx(1);
    }

    //edges
    edge_n = KinematicVectorArray(x_n.size(), job->JOB_TYPE);
    node_tag = Eigen::VectorXi(x_n.size());
    node_tag.setZero();
    for (int n=0;n<x_n.size();n++){
        tmp = n;
        for (int i=0;i<2;i++){
            if (i == 1) {
                edge_n(n,i) = 2;
                continue;
            }
            ijk(i) = tmp % (Nx(i)+1);
            tmp = tmp / (Nx(i)+1);

            if (ijk(i) > 1 && ijk(i) < (Nx(i)-1)){
                edge_n(n,i) = 2; //interior
            } else if (ijk(i) == 1 || ijk(i) == (Nx(i)-1)) {
                if (ijk(i) == 1) {
                    edge_n(n, i) = 1; //next to boundary at start
                } else {
                    edge_n(n, i) = -1; //next to boundary at end
                }
            } else {
                edge_n(n,i) = 0; //boundary
                if (ijk(i) == 0){
                    node_tag(n) = -1;
                } else {
                    node_tag(n) = 1;
                }
            }
        }
    }

    //nodal volume
    v_n.resize(x_n.size());
    v_n.setZero();

    //use 3-pt Gaussian Quadrature to determine nodal value (up to 4th order bases)
    double x0; //element minimum radius
    double xtmp, rtmp; //gauss point in radial coords
    double x1 = 0;
    double x2 = std::sqrt(3.0/5.0);
    double x3 = -x2;
    double w1 = 8.0/9.0;
    double w2 = 5.0/9.0;
    double w3 = w2;
    KinematicVector diff_x_rtz;

    for (int e=0;e<nodeIDs.rows();e++){
        //x0 = x_n(nodeIDs(e,0))[0];
        for (int c=0;c<nodeIDs.cols();c++) {
            if (nodeIDs(e,c) < 0){
                //ignore these
                continue;
            }

            rtz = cPoint_to_rPoint(x_n(nodeIDs(e,5)));
            diff_x_rtz = cPoint_to_rPoint(x_n(nodeIDs(e,5))) - cPoint_to_rPoint(x_n(nodeIDs(e,c)));

            xtmp = diff_x_rtz[0] + hx(0)/2 * (x1 + 1);
            rtmp = rtz[0] + hx(0)/2 * (x1 + 1);
            if (edge_n(nodeIDs(e,c),0) == 2) {
                //in bulk
                v_n(nodeIDs(e, c)) += hx(0)/2 * w1 * rtmp * s_cubic(xtmp, hx(0));
            } else if (edge_n(nodeIDs(e,c),0) == 1) {
                //next to outer edge
                v_n(nodeIDs(e, c)) += hx(0) / 2 * w1 * rtmp * (s_cubic(xtmp, hx(0)) - s_cubic(xtmp + 2 * hx(0), hx(0)));
            } else if (edge_n(nodeIDs(e,c),0) == -1) {
                //next to inner edge
                v_n(nodeIDs(e, c)) += hx(0) / 2 * w1 * rtmp * (s_cubic(xtmp, hx(0)) - s_cubic(xtmp - 2 * hx(0), hx(0)));
            } else if (edge_n(nodeIDs(e,c),0) == 0) {
                //edge
                v_n(nodeIDs(e, c)) += hx(0) / 2 * w1 * rtmp * (s_cubic(xtmp, hx(0)) + 2.0*s_cubic(std::abs(xtmp) + hx(0), hx(0)));
            }

            xtmp = diff_x_rtz[0] + hx(0)/2 * (x2 + 1);
            rtmp = rtz[0] + hx(0)/2 * (x2 + 1);
            if (edge_n(nodeIDs(e,c),0) == 2) {
                //in bulk
                v_n(nodeIDs(e, c)) += hx(0)/2 * w2 * rtmp * s_cubic(xtmp, hx(0));
            } else if (edge_n(nodeIDs(e,c),0) == 1) {
                //next to outer edge
                v_n(nodeIDs(e, c)) += hx(0) / 2 * w2 * rtmp * (s_cubic(xtmp, hx(0)) - s_cubic(xtmp + 2 * hx(0), hx(0)));
            } else if (edge_n(nodeIDs(e,c),0) == -1) {
                //next to inner edge
                v_n(nodeIDs(e, c)) += hx(0) / 2 * w2 * rtmp * (s_cubic(xtmp, hx(0)) - s_cubic(xtmp - 2 * hx(0), hx(0)));
            } else if (edge_n(nodeIDs(e,c),0) == 0) {
                //edge
                v_n(nodeIDs(e, c)) += hx(0) / 2 * w2 * rtmp * (s_cubic(xtmp, hx(0)) + 2.0*s_cubic(std::abs(xtmp) + hx(0), hx(0)));
            }

            xtmp = diff_x_rtz[0] + hx(0)/2 * (x3 + 1);
            rtmp = rtz[0] + hx(0)/2 * (x3 + 1);
            if (edge_n(nodeIDs(e,c),0) == 2) {
                //in bulk
                v_n(nodeIDs(e, c)) += hx(0)/2 * w3 * rtmp * s_cubic(xtmp, hx(0));
            } else if (edge_n(nodeIDs(e,c),0) == 1) {
                //next to outer edge
                v_n(nodeIDs(e, c)) += hx(0) / 2 * w3 * rtmp * (s_cubic(xtmp, hx(0)) - s_cubic(xtmp + 2 * hx(0), hx(0)));
            } else if (edge_n(nodeIDs(e,c),0) == -1) {
                //next to inner edge
                v_n(nodeIDs(e, c)) += hx(0) / 2 * w3 * rtmp * (s_cubic(xtmp, hx(0)) - s_cubic(xtmp - 2 * hx(0), hx(0)));
            } else if (edge_n(nodeIDs(e,c),0) == 0) {
                //edge
                v_n(nodeIDs(e, c)) += hx(0) / 2 * w3 * rtmp * (s_cubic(xtmp, hx(0)) + 2.0*s_cubic(std::abs(xtmp) + hx(0), hx(0)));
            }

            //std::cout << xtmp/hx(0) << std::endl;
        }
    }

    /*
    for (int i=0; i<21; i++){
        std::cout << nodeIDs.row(i) << std::endl;
        std::cout << v_n(nodeIDs(i,5)) << " ";
        std::cout << v_n(nodeIDs(i,6)) << " ";
        std::cout << v_n(nodeIDs(i,9)) << " ";
        std::cout << v_n(nodeIDs(i,10)) << std::endl;
    }
    */

    //scale v_n by theta contribution
    //each node will have integral repeated 4x
    v_n *= hx(1)/4.0;

    return;
}

/*----------------------------------------------------------------------------*/
//
void Regular2DTaylorCouetteCell::writeFrame(Job* job, Serializer* serializer){
    serializer->writeVectorArray(edge_n,"edge_id");
    return;
}

std::string Regular2DTaylorCouetteCell::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

int Regular2DTaylorCouetteCell::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

/*----------------------------------------------------------------------------*/
//
void Regular2DTaylorCouetteCell::writeHeader(Job* job, Body* body, Serializer* serializer, std::ofstream& nfile, int SPEC){
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

    if (order == LINEAR) {
        nfile << "CELLS " << element_count << " " << 5 * element_count << "\n";
        for (int e = 0; e < element_count; e++) {
            nfile << "4 " << nodeIDs(e, 0) << " " << nodeIDs(e, 1) << " " << nodeIDs(e, 3) << " " << nodeIDs(e, 2)
                  << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e = 0; e < element_count; e++) {
            nfile << "9\n";
        }
    } else if (order == CUBIC){
        nfile << "CELLS " << element_count << " " << 5 * element_count << "\n";
        for (int e = 0; e < element_count; e++){
            nfile << "4 " << nodeIDs(e,5) << " " << nodeIDs(e,6) << " " << nodeIDs(e,10) << " " << nodeIDs(e,9) << "\n";
        }

        nfile << "CELL_TYPES " << element_count << "\n";
        for (int e=0; e<element_count; e++){
            nfile << "9\n";
        }
    }

    nfile << "POINT_DATA " << nlen << "\n";
}

/*----------------------------------------------------------------------------*/
//
int Regular2DTaylorCouetteCell::whichElement(Job* job, KinematicVector& xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && xIN.VECTOR_TYPE == hx.VECTOR_TYPE && "whichElement failed");
    KinematicVector rtz = cPoint_to_rPoint(xIN); //convert to radial coords

    bool inDomain = (rtz[0] <= Ro && rtz[0] >= Ri);
    if (!inDomain) { //if xIn is outside domain, return -1
        return -1;
    }

    return floor((rtz[0] - Ri) / hx[0]) + floor(rtz[1] / hx[1]) * (Nx[0]); //id in x-dimension
}

bool Regular2DTaylorCouetteCell::inDomain(Job* job, KinematicVector& xIN){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && xIN.VECTOR_TYPE == hx.VECTOR_TYPE && "whichElement failed");
    KinematicVector rtz = cPoint_to_rPoint(xIN); //convert to radial coords

    return (rtz[0] <= Ro && rtz[0] >= Ri);
}

KinematicVector Regular2DTaylorCouetteCell::nodeIDToPosition(Job* job, int idIN){
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

    tmpVec[0] += Ri;

    return rPoint_to_cPoint(tmpVec);
}

/*----------------------------------------------------------------------------*/
//
void Regular2DTaylorCouetteCell::evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && "evaluateShapeFnValue failed");
    KinematicVector rtz = cPoint_to_rPoint(xIN);
    KinematicVector rtz_node, diff_rtz;
    KinematicVector rst(job->JOB_TYPE);
    double tmp = 1.0;
    int elementID = whichElement(job,xIN);
    if (elementID < 0){
        return;
    }
    for (int n=0;n<nodeIDs.cols();n++){
        if (nodeIDs(elementID,n) == -1){
            continue; //break out if node id is -1;
        }
        //find local coordinates relative to nodal position
        rtz_node = cPoint_to_rPoint(x_n(nodeIDs(elementID,n)));
        diff_rtz = rtz - rtz_node;

        //check for wrapping
        if (std::abs(rtz[1] - rtz_node[1] + 2*pi) < std::abs(diff_rtz[1])){
            diff_rtz[1] = rtz[1] - rtz_node[1] + 2*pi;
        } else if (std::abs(rtz[1] - rtz_node[1] - 2*pi) < std::abs(diff_rtz[1])){
            diff_rtz[1] = rtz[1] - rtz_node[1] - 2*pi;
        }


        if (order == LINEAR) {
            for (int i = 0; i < rtz.DIM; i++) {
                //r = (x_p - x_n)/hx
                rst[i] = (diff_rtz[i]) / hx(i);

                //standard linear hat function
                tmp *= (1 - std::abs(rst[i]));
            }
        } else if (order == CUBIC){
            for (int i = 0; i < rtz.DIM; i++) {
                //r = (x_p - x_n)/hx
                rst[i] = (diff_rtz[i]); //is scaled when input to s_cubic

                if (edge_n(nodeIDs(elementID, n), i) == 2) {
                    tmp *= s_cubic(rst(i), hx(i));
                } else if (edge_n(nodeIDs(elementID, n), i) == 1) {
                    tmp *= (s_cubic(rst(i), hx(i)) - s_cubic(rst(i) + 2 * hx(i), hx(i)));
                } else if (edge_n(nodeIDs(elementID, n), i) == -1) {
                    tmp *= (s_cubic(rst(i), hx(i)) - s_cubic(rst(i) - 2 * hx(i), hx(i)));
                } else if (edge_n(nodeIDs(elementID, n), i) == 0) {
                    tmp *= (s_cubic(rst(i), hx(i)) + 2.0 * s_cubic(std::abs(rst(i)) + hx(i), hx(i)));
                }
            }
        }

        nID.push_back(nodeIDs(elementID,n));
        nVAL.push_back(tmp);
        tmp = 1.0;
    }
    return;
}
void Regular2DTaylorCouetteCell::evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD){
    assert(xIN.VECTOR_TYPE == Lx.VECTOR_TYPE && nGRAD.VECTOR_TYPE == Lx.VECTOR_TYPE && "evaluateShapeFnGradient failed");
    KinematicVector rtz = cPoint_to_rPoint(xIN);
    KinematicVector rtz_node, diff_rtz;
    KinematicVector rst(job->JOB_TYPE);
    KinematicVector tmpVec(job->JOB_TYPE);
    double tmp = 1.0;
    int elementID = whichElement(job,xIN);
    if (elementID < 0){
        return;
    }
    for (int n=0;n<nodeIDs.cols();n++){
        //find local coordinates relative to nodal position
        if (nodeIDs(elementID,n) == -1){
            continue; //break out if node id is -1;
        }

        rtz_node = cPoint_to_rPoint(x_n(nodeIDs(elementID,n)));
        diff_rtz = rtz - rtz_node;

        //check for wrapping
        if (std::abs(rtz[1] - rtz_node[1] + 2*pi) < std::abs(diff_rtz[1])){
            diff_rtz[1] = rtz[1] - rtz_node[1] + 2*pi;
        } else if (std::abs(rtz[1] - rtz_node[1] - 2*pi) < std::abs(diff_rtz[1])){
            diff_rtz[1] = rtz[1] - rtz_node[1] - 2*pi;
        }

        if (order == LINEAR) {
            for (int i = 0; i < rtz.DIM; i++) {
                //r = (x_p - x_n)/hx
                rst[i] = (diff_rtz[i]) / hx(i);

                //standard linear hat function
                //evaluate at point
                tmp *= (1 - std::abs(rst(i)));
            }
            for (int i = 0; i < rtz.rows(); i++) {
                //replace i-direction contribution with sign function
                tmpVec(i) = -tmp / (1 - std::abs(rst(i))) * rst(i) / std::abs(rst(i)) / hx(i);
            }
        } else if (order = CUBIC){
            //find local coordinates relative to nodal position
            //r = (x_p - x_n)
            rst = diff_rtz;
            for (int i=0;i<xIN.rows();i++){
                //check proximity to edge
                if (edge_n(nodeIDs(elementID,n),i) == 2) {
                    tmpVec(i) = g_cubic(rst(i), hx(i));
                } else if (edge_n(nodeIDs(elementID,n),i) == 1) {
                    tmpVec(i) = (g_cubic(rst(i), hx(i)) - g_cubic(rst(i)+2*hx(i), hx(i)));
                } else if (edge_n(nodeIDs(elementID,n),i) == -1) {
                    tmpVec(i) = (g_cubic(rst(i), hx(i)) - g_cubic(rst(i)-2*hx(i), hx(i)));
                } else if (edge_n(nodeIDs(elementID,n),i) == 0) {
                    if (rst(i) >= 0) {
                        tmpVec(i) = (g_cubic(rst(i), hx(i)) + 2.0 * g_cubic(rst(i) + hx(i), hx(i)));
                    } else if (rst(i) < 0){
                        tmpVec(i) = (g_cubic(rst(i), hx(i)) + 2.0 * g_cubic(rst(i) - hx(i), hx(i)));
                    }
                }
            }

            for (int i=0;i<xIN.rows();i++){
                //check proximity to edge
                if (edge_n(nodeIDs(elementID,n),i) == 2) {
                    tmpVec *= s_cubic(rst(i), hx(i));
                    tmpVec(i) *= 1.0/s_cubic(rst(i), hx(i));
                } else if (edge_n(nodeIDs(elementID,n),i) == 1) {
                    tmpVec *= (s_cubic(rst(i), hx(i)) - s_cubic(rst(i)+2*hx(i), hx(i)));
                    tmpVec(i) *= 1.0/(s_cubic(rst(i), hx(i)) - s_cubic(rst(i)+2*hx(i), hx(i)));
                } else if (edge_n(nodeIDs(elementID,n),i) == -1) {
                    tmpVec *= (s_cubic(rst(i), hx(i)) - s_cubic(rst(i)-2*hx(i), hx(i)));
                    tmpVec(i) *= 1.0/(s_cubic(rst(i), hx(i)) - s_cubic(rst(i)-2*hx(i), hx(i)));
                } else if (edge_n(nodeIDs(elementID,n),i) == 0) {
                    tmpVec *= (s_cubic(rst(i), hx(i)) + 2.0*s_cubic(std::abs(rst(i))+hx(i), hx(i)));
                    tmpVec(i) *= 1.0/(s_cubic(rst(i), hx(i)) + 2.0*s_cubic(std::abs(rst(i))+hx(i), hx(i)));
                }
            }
        }

        //convert from dtheta to dx_theta
        tmpVec(1) /= rtz(0);

        //convert from local radial coords to cartesian frame
        tmpVec = rVector_to_cVector(tmpVec, rtz);

        nID.push_back(nodeIDs(elementID,n));
        nGRAD.push_back(tmpVec);
        tmp = 1.0;
        tmpVec.setZero();
    }
    return;
}

double Regular2DTaylorCouetteCell::nodeVolume(Job* job, int idIN){
    return v_n(idIN);
}

double Regular2DTaylorCouetteCell::elementVolume(Job* job, int idIN){
    return v_e(idIN);
}

int Regular2DTaylorCouetteCell::nodeTag(Job* job, int idIN){
    return node_tag(idIN);
}
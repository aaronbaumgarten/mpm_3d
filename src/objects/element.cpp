//
// Created by aaron on 12/20/16.
// element.cpp
//

#include <stdlib.h>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <cmath>

#include "node.hpp"
#include "body.hpp"
#include "particle.hpp"
#include "element.hpp"
#include "process.hpp"

Elements::Elements(size_t e, size_t npe){
    numElements = e;
    numNodesPerElement = npe;
    nodeID.resize(e,npe);
}

void Elements::resizeElements(size_t e, size_t npe){
    numElements = e;
    numNodesPerElement = npe;
    nodeID.resize(e,npe);
}

void Elements::addElement(Eigen::VectorXi nodeIDsIn, size_t idIn) {
    //assert that node id vector has same number as number of nodes per element
    assert(nodeIDsIn.size()==this->nodeID.cols());
    assert(idIn<nodeID.rows());

    this->nodeID.block(idIn,0,1,nodeIDsIn.size()) << nodeIDsIn.transpose();

    return;
}

void Elements::calculatePhic(Body *body, size_t ide, size_t idp, size_t idc, int use_cpdi){
    //Phic = 1/8*A_i(r^p_c)
    //gPhic = hx^3/8*A_i(r^p_c)*gradA_i(r*^p)
    Eigen::VectorXi nID(8);
    nID << this->nodeID.block(ide,0,1,8).transpose();
    double r,s,t;
    double dx = body->nodes.x[nID[7]] - body->nodes.x[nID[0]];
    double dy = body->nodes.y[nID[7]] - body->nodes.y[nID[0]];
    double dz = body->nodes.z[nID[7]] - body->nodes.z[nID[0]];
    if (dx <= 0 || dy<=0 || dz<=0){
        std::cout << "Error! h <= 0 in element: " << ide << "\n";
    } else {
        r = (body->particles.corner_x(idp,idc) - body->nodes.x[nID[0]])/dx;
        s = (body->particles.corner_y(idp,idc) - body->nodes.y[nID[0]])/dy;
        t = (body->particles.corner_z(idp,idc) - body->nodes.z[nID[0]])/dz;

        if (r<0 || s<0 || t<0) {
            std::cout << "Error! <r,s,t> has negative component!\n";
        }

        double a = body->particles.a[idp];
        double v = body->particles.v_averaging[idp];
        double Ap[3][8] = {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};

        body->PhiTriplets.emplace_back(nID[0],idp,1.0/8.0*S_L(r)*S_L(s)*S_L(t));//(1-r)*(1-s)*(1-t));
        body->PhiTriplets.emplace_back(nID[1],idp,1.0/8.0*S_R(r)*S_L(s)*S_L(t));//(1+r)*(1-s)*(1-t));
        body->PhiTriplets.emplace_back(nID[2],idp,1.0/8.0*S_L(r)*S_R(s)*S_L(t));//(1-r)*(1+s)*(1-t));
        body->PhiTriplets.emplace_back(nID[3],idp,1.0/8.0*S_R(r)*S_R(s)*S_L(t));//(1+r)*(1+s)*(1-t));
        body->PhiTriplets.emplace_back(nID[4],idp,1.0/8.0*S_L(r)*S_L(s)*S_R(t));//(1-r)*(1-s)*(1+t));
        body->PhiTriplets.emplace_back(nID[5],idp,1.0/8.0*S_R(r)*S_L(s)*S_R(t));//(1+r)*(1-s)*(1+t));
        body->PhiTriplets.emplace_back(nID[6],idp,1.0/8.0*S_L(r)*S_R(s)*S_R(t));//(1-r)*(1+s)*(1+t));
        body->PhiTriplets.emplace_back(nID[7],idp,1.0/8.0*S_R(r)*S_R(s)*S_R(t));//(1+r)*(1+s)*(1+t));

        if (use_cpdi==0){
            //calculate gradient explicitly
            body->gradPhiXTriplets.emplace_back(nID[0],idp , -0.125/dx*S_L(s)*S_L(t));//(1-r)*(1-s)*(1-t)*(1-sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nID[1],idp , 0.125/dx*S_L(s)*S_L(t));//(1+r)*(1-s)*(1-t)*(1-sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nID[2],idp , -0.125/dx*S_R(s)*S_L(t));//(1-r)*(1+s)*(1-t)*(1+sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nID[3],idp , 0.125/dx*S_R(s)*S_L(t));//(1+r)*(1+s)*(1-t)*(1+sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nID[4],idp , -0.125/dx*S_L(s)*S_R(t));//(1-r)*(1-s)*(1+t)*(1-sp)*(1+tp));
            body->gradPhiXTriplets.emplace_back(nID[5],idp , 0.125/dx*S_L(s)*S_R(t));//(1+r)*(1-s)*(1+t)*(1-sp)*(1+tp));
            body->gradPhiXTriplets.emplace_back(nID[6],idp , -0.125/dx*S_R(s)*S_R(t));//(1-r)*(1+s)*(1+t)*(1+sp)*(1+tp));
            body->gradPhiXTriplets.emplace_back(nID[7],idp , 0.125/dx*S_R(s)*S_R(t));//(1+r)*(1+s)*(1+t)*(1+sp)*(1+tp));

            body->gradPhiYTriplets.emplace_back(nID[0],idp , S_L(r)*-0.125/dy*S_L(t));//(1-r)*(1-s)*(1-t)*(1-rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nID[1],idp , S_R(r)*-0.125/dy*S_L(t));//(1+r)*(1-s)*(1-t)*(1+rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nID[2],idp , S_L(r)*0.125/dy*S_L(t));//(1-r)*(1+s)*(1-t)*(1-rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nID[3],idp , S_R(r)*0.125/dy*S_L(t));//(1+r)*(1+s)*(1-t)*(1+rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nID[4],idp , S_L(r)*-0.125/dy*S_R(t));//(1-r)*(1-s)*(1+t)*(1-rp)*(1+tp));
            body->gradPhiYTriplets.emplace_back(nID[5],idp , S_R(r)*-0.125/dy*S_R(t));//(1+r)*(1-s)*(1+t)*(1+rp)*(1+tp));
            body->gradPhiYTriplets.emplace_back(nID[6],idp , S_L(r)*0.125/dy*S_R(t));//(1-r)*(1+s)*(1+t)*(1-rp)*(1+tp));
            body->gradPhiYTriplets.emplace_back(nID[7],idp , S_R(r)*0.125/dy*S_R(t));//(1+r)*(1+s)*(1+t)*(1+rp)*(1+tp));

            body->gradPhiZTriplets.emplace_back(nID[0],idp , S_L(r)*S_L(s)*-0.125/dz);//(1-r)*(1-s)*(1-t)*(1-sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nID[1],idp , S_R(r)*S_L(s)*-0.125/dz);//(1+r)*(1-s)*(1-t)*(1-sp)*(1+rp));
            body->gradPhiZTriplets.emplace_back(nID[2],idp , S_L(r)*S_R(s)*-0.125/dz);//(1-r)*(1+s)*(1-t)*(1+sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nID[3],idp , S_R(r)*S_R(s)*-0.125/dz);//(1+r)*(1+s)*(1-t)*(1+sp)*(1+rp));
            body->gradPhiZTriplets.emplace_back(nID[4],idp , S_L(r)*S_L(s)*0.125/dz);//(1-r)*(1-s)*(1+t)*(1-sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nID[5],idp , S_R(r)*S_L(s)*0.125/dz);//(1+r)*(1-s)*(1+t)*(1-sp)*(1+rp));
            body->gradPhiZTriplets.emplace_back(nID[6],idp , S_L(r)*S_R(s)*0.125/dz);//(1-r)*(1+s)*(1+t)*(1+sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nID[7],idp , S_R(r)*S_R(s)*0.125/dz);//(1+r)*(1+s)*(1+t)*(1+sp)*(1+rp));
        } else {
            //use cpdi method of element crossover
            body->gradPhiXTriplets.emplace_back(nID[0],idp , 1/(8*a)*S_L(r)*S_L(s)*S_L(t)*Ap[0][idc]);//(1-r)*(1-s)*(1-t)*(1-sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nID[1],idp , 1/(8*a)*S_R(r)*S_L(s)*S_L(t)*Ap[0][idc]);//(1+r)*(1-s)*(1-t)*(1-sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nID[2],idp , 1/(8*a)*S_L(r)*S_R(s)*S_L(t)*Ap[0][idc]);//(1-r)*(1+s)*(1-t)*(1+sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nID[3],idp , 1/(8*a)*S_R(r)*S_R(s)*S_L(t)*Ap[0][idc]);//(1+r)*(1+s)*(1-t)*(1+sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nID[4],idp , 1/(8*a)*S_L(r)*S_L(s)*S_R(t)*Ap[0][idc]);//(1-r)*(1-s)*(1+t)*(1-sp)*(1+tp));
            body->gradPhiXTriplets.emplace_back(nID[5],idp , 1/(8*a)*S_R(r)*S_L(s)*S_R(t)*Ap[0][idc]);//(1+r)*(1-s)*(1+t)*(1-sp)*(1+tp));
            body->gradPhiXTriplets.emplace_back(nID[6],idp , 1/(8*a)*S_L(r)*S_R(s)*S_R(t)*Ap[0][idc]);//(1-r)*(1+s)*(1+t)*(1+sp)*(1+tp));
            body->gradPhiXTriplets.emplace_back(nID[7],idp , 1/(8*a)*S_R(r)*S_R(s)*S_R(t)*Ap[0][idc]);//(1+r)*(1+s)*(1+t)*(1+sp)*(1+tp));

            body->gradPhiYTriplets.emplace_back(nID[0],idp , 1/(8*a)*S_L(r)*S_L(s)*S_L(t)*Ap[1][idc]);//(1-r)*(1-s)*(1-t)*(1-rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nID[1],idp , 1/(8*a)*S_R(r)*S_L(s)*S_L(t)*Ap[1][idc]);//(1+r)*(1-s)*(1-t)*(1+rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nID[2],idp , 1/(8*a)*S_L(r)*S_R(s)*S_L(t)*Ap[1][idc]);//(1-r)*(1+s)*(1-t)*(1-rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nID[3],idp , 1/(8*a)*S_R(r)*S_R(s)*S_L(t)*Ap[1][idc]);//(1+r)*(1+s)*(1-t)*(1+rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nID[4],idp , 1/(8*a)*S_L(r)*S_L(s)*S_R(t)*Ap[1][idc]);//(1-r)*(1-s)*(1+t)*(1-rp)*(1+tp));
            body->gradPhiYTriplets.emplace_back(nID[5],idp , 1/(8*a)*S_R(r)*S_L(s)*S_R(t)*Ap[1][idc]);//(1+r)*(1-s)*(1+t)*(1+rp)*(1+tp));
            body->gradPhiYTriplets.emplace_back(nID[6],idp , 1/(8*a)*S_L(r)*S_R(s)*S_R(t)*Ap[1][idc]);//(1-r)*(1+s)*(1+t)*(1-rp)*(1+tp));
            body->gradPhiYTriplets.emplace_back(nID[7],idp , 1/(8*a)*S_R(r)*S_R(s)*S_R(t)*Ap[1][idc]);//(1+r)*(1+s)*(1+t)*(1+rp)*(1+tp));

            body->gradPhiZTriplets.emplace_back(nID[0],idp , 1/(8*a)*S_L(r)*S_L(s)*S_L(t)*Ap[2][idc]);//(1-r)*(1-s)*(1-t)*(1-sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nID[1],idp , 1/(8*a)*S_R(r)*S_L(s)*S_L(t)*Ap[2][idc]);//(1+r)*(1-s)*(1-t)*(1-sp)*(1+rp));
            body->gradPhiZTriplets.emplace_back(nID[2],idp , 1/(8*a)*S_L(r)*S_R(s)*S_L(t)*Ap[2][idc]);//(1-r)*(1+s)*(1-t)*(1+sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nID[3],idp , 1/(8*a)*S_R(r)*S_R(s)*S_L(t)*Ap[2][idc]);//(1+r)*(1+s)*(1-t)*(1+sp)*(1+rp));
            body->gradPhiZTriplets.emplace_back(nID[4],idp , 1/(8*a)*S_L(r)*S_L(s)*S_R(t)*Ap[2][idc]);//(1-r)*(1-s)*(1+t)*(1-sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nID[5],idp , 1/(8*a)*S_R(r)*S_L(s)*S_R(t)*Ap[2][idc]);//(1+r)*(1-s)*(1+t)*(1-sp)*(1+rp));
            body->gradPhiZTriplets.emplace_back(nID[6],idp , 1/(8*a)*S_L(r)*S_R(s)*S_R(t)*Ap[2][idc]);//(1-r)*(1+s)*(1+t)*(1+sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nID[7],idp , 1/(8*a)*S_R(r)*S_R(s)*S_R(t)*Ap[2][idc]);//(1+r)*(1+s)*(1+t)*(1+sp)*(1+rp));
        }
    }
    return;
}
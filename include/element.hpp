//
// Created by aaron on 8/26/16.
// element.hpp
//

#ifndef MPM_3D_ELEMENT_HPP
#define MPM_3D_ELEMENT_HPP

#include <stdlib.h>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <cmath>

#define B_r(x) ((std::abs(x)<=1.0)?(1.0-std::abs(x)):0.0) //right hand shape function
#define B_l(x) ((std::abs(x-1.0)<=1.0)?(1.0-std::abs(x-1.0)):0.0) //left hand shape function

class Element{
public:
    //unique id
    size_t id;

    //node info
    size_t numNodes;
    std::vector<size_t> nodeID;

    //particle info
    size_t numParticles;
    std::vector<size_t> particleID;

    //element color (for threading)
    //int color;

    Element() {}

    Element(size_t nn, size_t* nodeIDs, size_t idIn):
            id(idIn),
            numNodes(nn),
            nodeID(nn)
    {
        for (int i=0;i<nn;i++){
            nodeID[i] = nodeIDs[i];
        }
    }

    template<class bodyT, class particleT>
    void calculateSipc(bodyT *body, particleT *particle, size_t idc){
        //Sipc = 1/8*A_i(r^p_c)
        //gSipc = hx^3/8*A_i(r^p_c)*gradA_i(r*^p)
        double r,s,t;
        double rp,sp,tp;
        double r1,r2,s1,s2,t1,t2;
        double dx = body->nodes[nodeID[1]].x[0] - body->nodes[nodeID[0]].x[0];
        if (dx <= 0){
            std::cout << "Error! dx <= 0 in element.calculateSipc()\n";
        } else {
            r = (particle->corner[idc][0] - body->nodes[nodeID[0]].x[0])/dx;
            s = (particle->corner[idc][1] - body->nodes[nodeID[0]].y[0])/dx;
            t = (particle->corner[idc][2] - body->nodes[nodeID[0]].z[0])/dx;

            if (r<0 || s<0 || t<0) {
                std::cout << "Error! <r,s,t> has negative component!\n";
            }

            rp = (particle->x[0] - body->nodes[nodeID[0]].x[0])/dx;
            sp = (particle->y[0] - body->nodes[nodeID[0]].y[0])/dx;
            tp = (particle->z[0] - body->nodes[nodeID[0]].z[0])/dx;

            r1 = std::fmin(r,2.0*rp-r);
            r2 = std::fmax(r,2.0*rp-r);
            s1 = std::fmin(s,2.0*sp-s);
            s2 = std::fmax(s,2.0*sp-s);
            t1 = std::fmin(t,2.0*tp-t);
            t2 = std::fmax(t,2.0*tp-t);

            //if (rp<0 || sp<0 || tp<0) {
            //    std::cout << "Error! <rp,sp,tp> has negative component!\n";
            //}

            body->SipTriplets.emplace_back(nodeID[0],particle->id,1.0/8.0*B_l(r)*B_l(s)*B_l(t));//(1-r)*(1-s)*(1-t));
            body->SipTriplets.emplace_back(nodeID[1],particle->id,1.0/8.0*B_r(r)*B_l(s)*B_l(t));//(1+r)*(1-s)*(1-t));
            body->SipTriplets.emplace_back(nodeID[2],particle->id,1.0/8.0*B_l(r)*B_r(s)*B_l(t));//(1-r)*(1+s)*(1-t));
            body->SipTriplets.emplace_back(nodeID[3],particle->id,1.0/8.0*B_r(r)*B_r(s)*B_l(t));//(1+r)*(1+s)*(1-t));
            body->SipTriplets.emplace_back(nodeID[4],particle->id,1.0/8.0*B_l(r)*B_l(s)*B_r(t));//(1-r)*(1-s)*(1+t));
            body->SipTriplets.emplace_back(nodeID[5],particle->id,1.0/8.0*B_r(r)*B_l(s)*B_r(t));//(1+r)*(1-s)*(1+t));
            body->SipTriplets.emplace_back(nodeID[6],particle->id,1.0/8.0*B_l(r)*B_r(s)*B_r(t));//(1-r)*(1+s)*(1+t));
            body->SipTriplets.emplace_back(nodeID[7],particle->id,1.0/8.0*B_r(r)*B_r(s)*B_r(t));//(1+r)*(1+s)*(1+t));

            body->gradSipXTriplets.emplace_back(nodeID[0],particle->id , -dx*dx*dx/8.0*B_l(r)*B_l(s)*B_l(t)*(1.0/4.0)*(B_l(r1)-B_l(r2))*(B_l(s1)+B_l(s2))*(B_l(t1)+B_l(t2)));//(1-r)*(1-s)*(1-t)*(1-sp)*(1-tp));
            body->gradSipXTriplets.emplace_back(nodeID[1],particle->id , -dx*dx*dx/8.0*B_r(r)*B_l(s)*B_l(t)*(1.0/4.0)*(B_r(r1)-B_r(r2))*(B_l(s1)+B_l(s2))*(B_l(t1)+B_l(t2)));//(1+r)*(1-s)*(1-t)*(1-sp)*(1-tp));
            body->gradSipXTriplets.emplace_back(nodeID[2],particle->id , -dx*dx*dx/8.0*B_l(r)*B_r(s)*B_l(t)*(1.0/4.0)*(B_l(r1)-B_l(r2))*(B_r(s1)+B_r(s2))*(B_l(t1)+B_l(t2)));//(1-r)*(1+s)*(1-t)*(1+sp)*(1-tp));
            body->gradSipXTriplets.emplace_back(nodeID[3],particle->id , -dx*dx*dx/8.0*B_r(r)*B_r(s)*B_l(t)*(1.0/4.0)*(B_r(r1)-B_r(r2))*(B_r(s1)+B_r(s2))*(B_l(t1)+B_l(t2)));//(1+r)*(1+s)*(1-t)*(1+sp)*(1-tp));
            body->gradSipXTriplets.emplace_back(nodeID[4],particle->id , -dx*dx*dx/8.0*B_l(r)*B_l(s)*B_r(t)*(1.0/4.0)*(B_l(r1)-B_l(r2))*(B_l(s1)+B_l(s2))*(B_r(t1)+B_r(t2)));//(1-r)*(1-s)*(1+t)*(1-sp)*(1+tp));
            body->gradSipXTriplets.emplace_back(nodeID[5],particle->id , -dx*dx*dx/8.0*B_r(r)*B_l(s)*B_r(t)*(1.0/4.0)*(B_r(r1)-B_r(r2))*(B_l(s1)+B_l(s2))*(B_r(t1)+B_r(t2)));//(1+r)*(1-s)*(1+t)*(1-sp)*(1+tp));
            body->gradSipXTriplets.emplace_back(nodeID[6],particle->id , -dx*dx*dx/8.0*B_l(r)*B_r(s)*B_r(t)*(1.0/4.0)*(B_l(r1)-B_l(r2))*(B_r(s1)+B_r(s2))*(B_r(t1)+B_r(t2)));//(1-r)*(1+s)*(1+t)*(1+sp)*(1+tp));
            body->gradSipXTriplets.emplace_back(nodeID[7],particle->id , -dx*dx*dx/8.0*B_r(r)*B_r(s)*B_r(t)*(1.0/4.0)*(B_r(r1)-B_r(r2))*(B_r(s1)+B_r(s2))*(B_r(t1)+B_r(t2)));//(1+r)*(1+s)*(1+t)*(1+sp)*(1+tp));

            body->gradSipYTriplets.emplace_back(nodeID[0],particle->id , -dx*dx*dx/8.0*B_l(r)*B_l(s)*B_l(t)*(1.0/4.0)*(B_l(r1)+B_l(r2))*(B_l(s1)-B_l(s2))*(B_l(t1)+B_l(t2)));//(1-r)*(1-s)*(1-t)*(1-rp)*(1-tp));
            body->gradSipYTriplets.emplace_back(nodeID[1],particle->id , -dx*dx*dx/8.0*B_r(r)*B_l(s)*B_l(t)*(1.0/4.0)*(B_r(r1)+B_r(r2))*(B_l(s1)-B_l(s2))*(B_l(t1)+B_l(t2)));//(1+r)*(1-s)*(1-t)*(1+rp)*(1-tp));
            body->gradSipYTriplets.emplace_back(nodeID[2],particle->id , -dx*dx*dx/8.0*B_l(r)*B_r(s)*B_l(t)*(1.0/4.0)*(B_l(r1)+B_l(r2))*(B_r(s1)-B_r(s2))*(B_l(t1)+B_l(t2)));//(1-r)*(1+s)*(1-t)*(1-rp)*(1-tp));
            body->gradSipYTriplets.emplace_back(nodeID[3],particle->id , -dx*dx*dx/8.0*B_r(r)*B_r(s)*B_l(t)*(1.0/4.0)*(B_r(r1)+B_r(r2))*(B_r(s1)-B_r(s2))*(B_l(t1)+B_l(t2)));//(1+r)*(1+s)*(1-t)*(1+rp)*(1-tp));
            body->gradSipYTriplets.emplace_back(nodeID[4],particle->id , -dx*dx*dx/8.0*B_l(r)*B_l(s)*B_r(t)*(1.0/4.0)*(B_l(r1)+B_l(r2))*(B_l(s1)-B_l(s2))*(B_r(t1)+B_r(t2)));//(1-r)*(1-s)*(1+t)*(1-rp)*(1+tp));
            body->gradSipYTriplets.emplace_back(nodeID[5],particle->id , -dx*dx*dx/8.0*B_r(r)*B_l(s)*B_r(t)*(1.0/4.0)*(B_r(r1)+B_r(r2))*(B_l(s1)-B_l(s2))*(B_r(t1)+B_r(t2)));//(1+r)*(1-s)*(1+t)*(1+rp)*(1+tp));
            body->gradSipYTriplets.emplace_back(nodeID[6],particle->id , -dx*dx*dx/8.0*B_l(r)*B_r(s)*B_r(t)*(1.0/4.0)*(B_l(r1)+B_l(r2))*(B_r(s1)-B_r(s2))*(B_r(t1)+B_r(t2)));//(1-r)*(1+s)*(1+t)*(1-rp)*(1+tp));
            body->gradSipYTriplets.emplace_back(nodeID[7],particle->id , -dx*dx*dx/8.0*B_r(r)*B_r(s)*B_r(t)*(1.0/4.0)*(B_r(r1)+B_r(r2))*(B_r(s1)-B_r(s2))*(B_r(t1)+B_r(t2)));//(1+r)*(1+s)*(1+t)*(1+rp)*(1+tp));

            body->gradSipZTriplets.emplace_back(nodeID[0],particle->id , -dx*dx*dx/8.0*B_l(r)*B_l(s)*B_l(t)*(1.0/4.0)*(B_l(r1)+B_l(r2))*(B_l(s1)+B_l(s2))*(B_l(t1)-B_l(t2)));//(1-r)*(1-s)*(1-t)*(1-sp)*(1-rp));
            body->gradSipZTriplets.emplace_back(nodeID[1],particle->id , -dx*dx*dx/8.0*B_r(r)*B_l(s)*B_l(t)*(1.0/4.0)*(B_r(r1)+B_r(r2))*(B_l(s1)+B_l(s2))*(B_l(t1)-B_l(t2)));//(1+r)*(1-s)*(1-t)*(1-sp)*(1+rp));
            body->gradSipZTriplets.emplace_back(nodeID[2],particle->id , -dx*dx*dx/8.0*B_l(r)*B_r(s)*B_l(t)*(1.0/4.0)*(B_l(r1)+B_l(r2))*(B_r(s1)+B_r(s2))*(B_l(t1)-B_l(t2)));//(1-r)*(1+s)*(1-t)*(1+sp)*(1-rp));
            body->gradSipZTriplets.emplace_back(nodeID[3],particle->id , -dx*dx*dx/8.0*B_r(r)*B_r(s)*B_l(t)*(1.0/4.0)*(B_r(r1)+B_r(r2))*(B_r(s1)+B_r(s2))*(B_l(t1)-B_l(t2)));//(1+r)*(1+s)*(1-t)*(1+sp)*(1+rp));
            body->gradSipZTriplets.emplace_back(nodeID[4],particle->id , -dx*dx*dx/8.0*B_l(r)*B_l(s)*B_r(t)*(1.0/4.0)*(B_l(r1)+B_l(r2))*(B_l(s1)+B_l(s2))*(B_r(t1)-B_r(t2)));//(1-r)*(1-s)*(1+t)*(1-sp)*(1-rp));
            body->gradSipZTriplets.emplace_back(nodeID[5],particle->id , -dx*dx*dx/8.0*B_r(r)*B_l(s)*B_r(t)*(1.0/4.0)*(B_r(r1)+B_r(r2))*(B_l(s1)+B_l(s2))*(B_r(t1)-B_r(t2)));//(1+r)*(1-s)*(1+t)*(1-sp)*(1+rp));
            body->gradSipZTriplets.emplace_back(nodeID[6],particle->id , -dx*dx*dx/8.0*B_l(r)*B_r(s)*B_r(t)*(1.0/4.0)*(B_l(r1)+B_l(r2))*(B_r(s1)+B_r(s2))*(B_r(t1)-B_r(t2)));//(1-r)*(1+s)*(1+t)*(1+sp)*(1-rp));
            body->gradSipZTriplets.emplace_back(nodeID[7],particle->id , -dx*dx*dx/8.0*B_r(r)*B_r(s)*B_r(t)*(1.0/4.0)*(B_r(r1)+B_r(r2))*(B_r(s1)+B_r(s2))*(B_r(t1)-B_r(t2)));//(1+r)*(1+s)*(1+t)*(1+sp)*(1+rp));

        }
        return;
    };

};

#endif //MPM_3D_ELEMENT_HPP

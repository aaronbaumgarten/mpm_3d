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

#define S_l(x) ((std::abs(x)<=1.0)?(1.0-std::abs(x)):0.0) //right hand shape function
#define S_r(x) ((std::abs(x-1.0)<=1.0)?(1.0-std::abs(x-1.0)):0.0) //left hand shape function

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
    void calculatePhic(bodyT *body, particleT *particle, size_t idc){
        //Phic = 1/8*A_i(r^p_c)
        //gPhic = hx^3/8*A_i(r^p_c)*gradA_i(r*^p)
        double r,s,t;
        double rp,sp,tp;
        double r1,r2,s1,s2,t1,t2;
        double dx = body->nodes[nodeID[1]].x[0] - body->nodes[nodeID[0]].x[0];
        if (dx <= 0){
            std::cout << "Error! dx <= 0 in element.calculatePhic()\n";
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

            double a = particle->a[0];
            double v = particle->v_averaging[0];
            double Ap[3][8] = {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};

            //if (rp<0 || sp<0 || tp<0) {
            //    std::cout << "Error! <rp,sp,tp> has negative component!\n";
            //}

            body->PhiTriplets.emplace_back(nodeID[0],particle->id,1.0/8.0*S_l(r)*S_l(s)*S_l(t));//(1-r)*(1-s)*(1-t));
            body->PhiTriplets.emplace_back(nodeID[1],particle->id,1.0/8.0*S_r(r)*S_l(s)*S_l(t));//(1+r)*(1-s)*(1-t));
            body->PhiTriplets.emplace_back(nodeID[2],particle->id,1.0/8.0*S_l(r)*S_r(s)*S_l(t));//(1-r)*(1+s)*(1-t));
            body->PhiTriplets.emplace_back(nodeID[3],particle->id,1.0/8.0*S_r(r)*S_r(s)*S_l(t));//(1+r)*(1+s)*(1-t));
            body->PhiTriplets.emplace_back(nodeID[4],particle->id,1.0/8.0*S_l(r)*S_l(s)*S_r(t));//(1-r)*(1-s)*(1+t));
            body->PhiTriplets.emplace_back(nodeID[5],particle->id,1.0/8.0*S_r(r)*S_l(s)*S_r(t));//(1+r)*(1-s)*(1+t));
            body->PhiTriplets.emplace_back(nodeID[6],particle->id,1.0/8.0*S_l(r)*S_r(s)*S_r(t));//(1-r)*(1+s)*(1+t));
            body->PhiTriplets.emplace_back(nodeID[7],particle->id,1.0/8.0*S_r(r)*S_r(s)*S_r(t));//(1+r)*(1+s)*(1+t));

            body->gradPhiXTriplets.emplace_back(nodeID[0],particle->id , a*a/v*S_l(r)*S_l(s)*S_l(t)*Ap[0][idc]);//(1-r)*(1-s)*(1-t)*(1-sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nodeID[1],particle->id , a*a/v*S_r(r)*S_l(s)*S_l(t)*Ap[0][idc]);//(1+r)*(1-s)*(1-t)*(1-sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nodeID[2],particle->id , a*a/v*S_l(r)*S_r(s)*S_l(t)*Ap[0][idc]);//(1-r)*(1+s)*(1-t)*(1+sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nodeID[3],particle->id , a*a/v*S_r(r)*S_r(s)*S_l(t)*Ap[0][idc]);//(1+r)*(1+s)*(1-t)*(1+sp)*(1-tp));
            body->gradPhiXTriplets.emplace_back(nodeID[4],particle->id , a*a/v*S_l(r)*S_l(s)*S_r(t)*Ap[0][idc]);//(1-r)*(1-s)*(1+t)*(1-sp)*(1+tp));
            body->gradPhiXTriplets.emplace_back(nodeID[5],particle->id , a*a/v*S_r(r)*S_l(s)*S_r(t)*Ap[0][idc]);//(1+r)*(1-s)*(1+t)*(1-sp)*(1+tp));
            body->gradPhiXTriplets.emplace_back(nodeID[6],particle->id , a*a/v*S_l(r)*S_r(s)*S_r(t)*Ap[0][idc]);//(1-r)*(1+s)*(1+t)*(1+sp)*(1+tp));
            body->gradPhiXTriplets.emplace_back(nodeID[7],particle->id , a*a/v*S_r(r)*S_r(s)*S_r(t)*Ap[0][idc]);//(1+r)*(1+s)*(1+t)*(1+sp)*(1+tp));

            body->gradPhiYTriplets.emplace_back(nodeID[0],particle->id , a*a/v*S_l(r)*S_l(s)*S_l(t)*Ap[1][idc]);//(1-r)*(1-s)*(1-t)*(1-rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nodeID[1],particle->id , a*a/v*S_r(r)*S_l(s)*S_l(t)*Ap[1][idc]);//(1+r)*(1-s)*(1-t)*(1+rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nodeID[2],particle->id , a*a/v*S_l(r)*S_r(s)*S_l(t)*Ap[1][idc]);//(1-r)*(1+s)*(1-t)*(1-rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nodeID[3],particle->id , a*a/v*S_r(r)*S_r(s)*S_l(t)*Ap[1][idc]);//(1+r)*(1+s)*(1-t)*(1+rp)*(1-tp));
            body->gradPhiYTriplets.emplace_back(nodeID[4],particle->id , a*a/v*S_l(r)*S_l(s)*S_r(t)*Ap[1][idc]);//(1-r)*(1-s)*(1+t)*(1-rp)*(1+tp));
            body->gradPhiYTriplets.emplace_back(nodeID[5],particle->id , a*a/v*S_r(r)*S_l(s)*S_r(t)*Ap[1][idc]);//(1+r)*(1-s)*(1+t)*(1+rp)*(1+tp));
            body->gradPhiYTriplets.emplace_back(nodeID[6],particle->id , a*a/v*S_l(r)*S_r(s)*S_r(t)*Ap[1][idc]);//(1-r)*(1+s)*(1+t)*(1-rp)*(1+tp));
            body->gradPhiYTriplets.emplace_back(nodeID[7],particle->id , a*a/v*S_r(r)*S_r(s)*S_r(t)*Ap[1][idc]);//(1+r)*(1+s)*(1+t)*(1+rp)*(1+tp));

            body->gradPhiZTriplets.emplace_back(nodeID[0],particle->id , a*a/v*S_l(r)*S_l(s)*S_l(t)*Ap[2][idc]);//(1-r)*(1-s)*(1-t)*(1-sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nodeID[1],particle->id , a*a/v*S_r(r)*S_l(s)*S_l(t)*Ap[2][idc]);//(1+r)*(1-s)*(1-t)*(1-sp)*(1+rp));
            body->gradPhiZTriplets.emplace_back(nodeID[2],particle->id , a*a/v*S_l(r)*S_r(s)*S_l(t)*Ap[2][idc]);//(1-r)*(1+s)*(1-t)*(1+sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nodeID[3],particle->id , a*a/v*S_r(r)*S_r(s)*S_l(t)*Ap[2][idc]);//(1+r)*(1+s)*(1-t)*(1+sp)*(1+rp));
            body->gradPhiZTriplets.emplace_back(nodeID[4],particle->id , a*a/v*S_l(r)*S_l(s)*S_r(t)*Ap[2][idc]);//(1-r)*(1-s)*(1+t)*(1-sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nodeID[5],particle->id , a*a/v*S_r(r)*S_l(s)*S_r(t)*Ap[2][idc]);//(1+r)*(1-s)*(1+t)*(1-sp)*(1+rp));
            body->gradPhiZTriplets.emplace_back(nodeID[6],particle->id , a*a/v*S_l(r)*S_r(s)*S_r(t)*Ap[2][idc]);//(1-r)*(1+s)*(1+t)*(1+sp)*(1-rp));
            body->gradPhiZTriplets.emplace_back(nodeID[7],particle->id , a*a/v*S_r(r)*S_r(s)*S_r(t)*Ap[2][idc]);//(1+r)*(1+s)*(1+t)*(1+sp)*(1+rp));

        }
        return;
    };

};

#endif //MPM_3D_ELEMENT_HPP

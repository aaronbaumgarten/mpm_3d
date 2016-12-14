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

#define S_L(x) (1.0-(x))//right hand shape function
#define S_R(x) (x)//left hand shape function

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
    void calculatePhic(bodyT *body, particleT *particle, size_t idc, int use_cpdi){
        //Phic = 1/8*A_i(r^p_c)
        //gPhic = hx^3/8*A_i(r^p_c)*gradA_i(r*^p)
        double r,s,t;
        double rp,sp,tp;
        double dx = body->nodes[nodeID[7]].x[0] - body->nodes[nodeID[0]].x[0];
        double dy = body->nodes[nodeID[7]].y[0] - body->nodes[nodeID[0]].y[0];
        double dz = body->nodes[nodeID[7]].z[0] - body->nodes[nodeID[0]].z[0];
        if (dx <= 0 || dy<=0 || dz<=0){
            std::cout << "Error! h <= 0 in element: " << id << "\n";
            //std::cout << "Element: " << id << "\n";
            //std::cout << "nodeIDs: " << nodeID[0] << ", " << nodeID[1] << "\n";
        } else {
            r = (particle->corner[idc][0] - body->nodes[nodeID[0]].x[0])/dx;
            s = (particle->corner[idc][1] - body->nodes[nodeID[0]].y[0])/dy;
            t = (particle->corner[idc][2] - body->nodes[nodeID[0]].z[0])/dz;

            if (r<0 || s<0 || t<0) {
                std::cout << "Error! <r,s,t> has negative component!\n";
            }

            rp = (particle->x[0] - body->nodes[nodeID[0]].x[0])/dx;
            sp = (particle->y[0] - body->nodes[nodeID[0]].y[0])/dy;
            tp = (particle->z[0] - body->nodes[nodeID[0]].z[0])/dz;

            double a = particle->a[0];
            double v = particle->v_averaging[0];
            double Ap[3][8] = {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};

            //if (rp<0 || sp<0 || tp<0) {
            //    std::cout << "Error! <rp,sp,tp> has negative component!\n";
            //}

            body->PhiTriplets.emplace_back(nodeID[0],particle->id,1.0/8.0*S_L(r)*S_L(s)*S_L(t));//(1-r)*(1-s)*(1-t));
            body->PhiTriplets.emplace_back(nodeID[1],particle->id,1.0/8.0*S_R(r)*S_L(s)*S_L(t));//(1+r)*(1-s)*(1-t));
            body->PhiTriplets.emplace_back(nodeID[2],particle->id,1.0/8.0*S_L(r)*S_R(s)*S_L(t));//(1-r)*(1+s)*(1-t));
            body->PhiTriplets.emplace_back(nodeID[3],particle->id,1.0/8.0*S_R(r)*S_R(s)*S_L(t));//(1+r)*(1+s)*(1-t));
            body->PhiTriplets.emplace_back(nodeID[4],particle->id,1.0/8.0*S_L(r)*S_L(s)*S_R(t));//(1-r)*(1-s)*(1+t));
            body->PhiTriplets.emplace_back(nodeID[5],particle->id,1.0/8.0*S_R(r)*S_L(s)*S_R(t));//(1+r)*(1-s)*(1+t));
            body->PhiTriplets.emplace_back(nodeID[6],particle->id,1.0/8.0*S_L(r)*S_R(s)*S_R(t));//(1-r)*(1+s)*(1+t));
            body->PhiTriplets.emplace_back(nodeID[7],particle->id,1.0/8.0*S_R(r)*S_R(s)*S_R(t));//(1+r)*(1+s)*(1+t));

            if (use_cpdi==0 || std::abs(particle->x[0]-particle->corner[idc][0])<1e-7){
                body->gradPhiXTriplets.emplace_back(nodeID[0],particle->id , -0.125/dx*S_L(s)*S_L(t));//(1-r)*(1-s)*(1-t)*(1-sp)*(1-tp));
                body->gradPhiXTriplets.emplace_back(nodeID[1],particle->id , 0.125/dx*S_L(s)*S_L(t));//(1+r)*(1-s)*(1-t)*(1-sp)*(1-tp));
                body->gradPhiXTriplets.emplace_back(nodeID[2],particle->id , -0.125/dx*S_R(s)*S_L(t));//(1-r)*(1+s)*(1-t)*(1+sp)*(1-tp));
                body->gradPhiXTriplets.emplace_back(nodeID[3],particle->id , 0.125/dx*S_R(s)*S_L(t));//(1+r)*(1+s)*(1-t)*(1+sp)*(1-tp));
                body->gradPhiXTriplets.emplace_back(nodeID[4],particle->id , -0.125/dx*S_L(s)*S_R(t));//(1-r)*(1-s)*(1+t)*(1-sp)*(1+tp));
                body->gradPhiXTriplets.emplace_back(nodeID[5],particle->id , 0.125/dx*S_L(s)*S_R(t));//(1+r)*(1-s)*(1+t)*(1-sp)*(1+tp));
                body->gradPhiXTriplets.emplace_back(nodeID[6],particle->id , -0.125/dx*S_R(s)*S_R(t));//(1-r)*(1+s)*(1+t)*(1+sp)*(1+tp));
                body->gradPhiXTriplets.emplace_back(nodeID[7],particle->id , 0.125/dx*S_R(s)*S_R(t));//(1+r)*(1+s)*(1+t)*(1+sp)*(1+tp));

                body->gradPhiYTriplets.emplace_back(nodeID[0],particle->id , S_L(r)*-0.125/dy*S_L(t));//(1-r)*(1-s)*(1-t)*(1-rp)*(1-tp));
                body->gradPhiYTriplets.emplace_back(nodeID[1],particle->id , S_R(r)*-0.125/dy*S_L(t));//(1+r)*(1-s)*(1-t)*(1+rp)*(1-tp));
                body->gradPhiYTriplets.emplace_back(nodeID[2],particle->id , S_L(r)*0.125/dy*S_L(t));//(1-r)*(1+s)*(1-t)*(1-rp)*(1-tp));
                body->gradPhiYTriplets.emplace_back(nodeID[3],particle->id , S_R(r)*0.125/dy*S_L(t));//(1+r)*(1+s)*(1-t)*(1+rp)*(1-tp));
                body->gradPhiYTriplets.emplace_back(nodeID[4],particle->id , S_L(r)*-0.125/dy*S_R(t));//(1-r)*(1-s)*(1+t)*(1-rp)*(1+tp));
                body->gradPhiYTriplets.emplace_back(nodeID[5],particle->id , S_R(r)*-0.125/dy*S_R(t));//(1+r)*(1-s)*(1+t)*(1+rp)*(1+tp));
                body->gradPhiYTriplets.emplace_back(nodeID[6],particle->id , S_L(r)*0.125/dy*S_R(t));//(1-r)*(1+s)*(1+t)*(1-rp)*(1+tp));
                body->gradPhiYTriplets.emplace_back(nodeID[7],particle->id , S_R(r)*0.125/dy*S_R(t));//(1+r)*(1+s)*(1+t)*(1+rp)*(1+tp));

                body->gradPhiZTriplets.emplace_back(nodeID[0],particle->id , S_L(r)*S_L(s)*-0.125/dz);//(1-r)*(1-s)*(1-t)*(1-sp)*(1-rp));
                body->gradPhiZTriplets.emplace_back(nodeID[1],particle->id , S_R(r)*S_L(s)*-0.125/dz);//(1+r)*(1-s)*(1-t)*(1-sp)*(1+rp));
                body->gradPhiZTriplets.emplace_back(nodeID[2],particle->id , S_L(r)*S_R(s)*-0.125/dz);//(1-r)*(1+s)*(1-t)*(1+sp)*(1-rp));
                body->gradPhiZTriplets.emplace_back(nodeID[3],particle->id , S_R(r)*S_R(s)*-0.125/dz);//(1+r)*(1+s)*(1-t)*(1+sp)*(1+rp));
                body->gradPhiZTriplets.emplace_back(nodeID[4],particle->id , S_L(r)*S_L(s)*0.125/dz);//(1-r)*(1-s)*(1+t)*(1-sp)*(1-rp));
                body->gradPhiZTriplets.emplace_back(nodeID[5],particle->id , S_R(r)*S_L(s)*0.125/dz);//(1+r)*(1-s)*(1+t)*(1-sp)*(1+rp));
                body->gradPhiZTriplets.emplace_back(nodeID[6],particle->id , S_L(r)*S_R(s)*0.125/dz);//(1-r)*(1+s)*(1+t)*(1+sp)*(1-rp));
                body->gradPhiZTriplets.emplace_back(nodeID[7],particle->id , S_R(r)*S_R(s)*0.125/dz);//(1+r)*(1+s)*(1+t)*(1+sp)*(1+rp));
            } else {
                body->gradPhiXTriplets.emplace_back(nodeID[0],particle->id , a*a/v*S_L(r)*S_L(s)*S_L(t)*Ap[0][idc]);//(1-r)*(1-s)*(1-t)*(1-sp)*(1-tp));
                body->gradPhiXTriplets.emplace_back(nodeID[1],particle->id , a*a/v*S_R(r)*S_L(s)*S_L(t)*Ap[0][idc]);//(1+r)*(1-s)*(1-t)*(1-sp)*(1-tp));
                body->gradPhiXTriplets.emplace_back(nodeID[2],particle->id , a*a/v*S_L(r)*S_R(s)*S_L(t)*Ap[0][idc]);//(1-r)*(1+s)*(1-t)*(1+sp)*(1-tp));
                body->gradPhiXTriplets.emplace_back(nodeID[3],particle->id , a*a/v*S_R(r)*S_R(s)*S_L(t)*Ap[0][idc]);//(1+r)*(1+s)*(1-t)*(1+sp)*(1-tp));
                body->gradPhiXTriplets.emplace_back(nodeID[4],particle->id , a*a/v*S_L(r)*S_L(s)*S_R(t)*Ap[0][idc]);//(1-r)*(1-s)*(1+t)*(1-sp)*(1+tp));
                body->gradPhiXTriplets.emplace_back(nodeID[5],particle->id , a*a/v*S_R(r)*S_L(s)*S_R(t)*Ap[0][idc]);//(1+r)*(1-s)*(1+t)*(1-sp)*(1+tp));
                body->gradPhiXTriplets.emplace_back(nodeID[6],particle->id , a*a/v*S_L(r)*S_R(s)*S_R(t)*Ap[0][idc]);//(1-r)*(1+s)*(1+t)*(1+sp)*(1+tp));
                body->gradPhiXTriplets.emplace_back(nodeID[7],particle->id , a*a/v*S_R(r)*S_R(s)*S_R(t)*Ap[0][idc]);//(1+r)*(1+s)*(1+t)*(1+sp)*(1+tp));

                body->gradPhiYTriplets.emplace_back(nodeID[0],particle->id , a*a/v*S_L(r)*S_L(s)*S_L(t)*Ap[1][idc]);//(1-r)*(1-s)*(1-t)*(1-rp)*(1-tp));
                body->gradPhiYTriplets.emplace_back(nodeID[1],particle->id , a*a/v*S_R(r)*S_L(s)*S_L(t)*Ap[1][idc]);//(1+r)*(1-s)*(1-t)*(1+rp)*(1-tp));
                body->gradPhiYTriplets.emplace_back(nodeID[2],particle->id , a*a/v*S_L(r)*S_R(s)*S_L(t)*Ap[1][idc]);//(1-r)*(1+s)*(1-t)*(1-rp)*(1-tp));
                body->gradPhiYTriplets.emplace_back(nodeID[3],particle->id , a*a/v*S_R(r)*S_R(s)*S_L(t)*Ap[1][idc]);//(1+r)*(1+s)*(1-t)*(1+rp)*(1-tp));
                body->gradPhiYTriplets.emplace_back(nodeID[4],particle->id , a*a/v*S_L(r)*S_L(s)*S_R(t)*Ap[1][idc]);//(1-r)*(1-s)*(1+t)*(1-rp)*(1+tp));
                body->gradPhiYTriplets.emplace_back(nodeID[5],particle->id , a*a/v*S_R(r)*S_L(s)*S_R(t)*Ap[1][idc]);//(1+r)*(1-s)*(1+t)*(1+rp)*(1+tp));
                body->gradPhiYTriplets.emplace_back(nodeID[6],particle->id , a*a/v*S_L(r)*S_R(s)*S_R(t)*Ap[1][idc]);//(1-r)*(1+s)*(1+t)*(1-rp)*(1+tp));
                body->gradPhiYTriplets.emplace_back(nodeID[7],particle->id , a*a/v*S_R(r)*S_R(s)*S_R(t)*Ap[1][idc]);//(1+r)*(1+s)*(1+t)*(1+rp)*(1+tp));

                body->gradPhiZTriplets.emplace_back(nodeID[0],particle->id , a*a/v*S_L(r)*S_L(s)*S_L(t)*Ap[2][idc]);//(1-r)*(1-s)*(1-t)*(1-sp)*(1-rp));
                body->gradPhiZTriplets.emplace_back(nodeID[1],particle->id , a*a/v*S_R(r)*S_L(s)*S_L(t)*Ap[2][idc]);//(1+r)*(1-s)*(1-t)*(1-sp)*(1+rp));
                body->gradPhiZTriplets.emplace_back(nodeID[2],particle->id , a*a/v*S_L(r)*S_R(s)*S_L(t)*Ap[2][idc]);//(1-r)*(1+s)*(1-t)*(1+sp)*(1-rp));
                body->gradPhiZTriplets.emplace_back(nodeID[3],particle->id , a*a/v*S_R(r)*S_R(s)*S_L(t)*Ap[2][idc]);//(1+r)*(1+s)*(1-t)*(1+sp)*(1+rp));
                body->gradPhiZTriplets.emplace_back(nodeID[4],particle->id , a*a/v*S_L(r)*S_L(s)*S_R(t)*Ap[2][idc]);//(1-r)*(1-s)*(1+t)*(1-sp)*(1-rp));
                body->gradPhiZTriplets.emplace_back(nodeID[5],particle->id , a*a/v*S_R(r)*S_L(s)*S_R(t)*Ap[2][idc]);//(1+r)*(1-s)*(1+t)*(1-sp)*(1+rp));
                body->gradPhiZTriplets.emplace_back(nodeID[6],particle->id , a*a/v*S_L(r)*S_R(s)*S_R(t)*Ap[2][idc]);//(1-r)*(1+s)*(1+t)*(1+sp)*(1-rp));
                body->gradPhiZTriplets.emplace_back(nodeID[7],particle->id , a*a/v*S_R(r)*S_R(s)*S_R(t)*Ap[2][idc]);//(1+r)*(1+s)*(1+t)*(1+sp)*(1+rp));
            }
        }
        return;
    };

};

#endif //MPM_3D_ELEMENT_HPP

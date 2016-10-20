//
// Created by aaron on 8/26/16.
// element.hpp
//

#ifndef MPM_3D_ELEMENT_HPP
#define MPM_3D_ELEMENT_HPP

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

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

            if (rp<0 || sp<0 || tp<0) {
                std::cout << "Error! <rp,sp,tp> has negative component!\n";
            }

            body->SipTriplets.emplace_back(Eigen::Triplet<double>(nodeID[0],particle->id,1.0/8.0*(1-r)*(1-s)*(1-t)));
            body->SipTriplets.emplace_back(Eigen::Triplet<double>(nodeID[1],particle->id,1.0/8.0*(1+r)*(1-s)*(1-t)));
            body->SipTriplets.emplace_back(Eigen::Triplet<double>(nodeID[2],particle->id,1.0/8.0*(1-r)*(1+s)*(1-t)));
            body->SipTriplets.emplace_back(Eigen::Triplet<double>(nodeID[3],particle->id,1.0/8.0*(1+r)*(1+s)*(1-t)));
            body->SipTriplets.emplace_back(Eigen::Triplet<double>(nodeID[4],particle->id,1.0/8.0*(1-r)*(1-s)*(1+t)));
            body->SipTriplets.emplace_back(Eigen::Triplet<double>(nodeID[5],particle->id,1.0/8.0*(1+r)*(1-s)*(1+t)));
            body->SipTriplets.emplace_back(Eigen::Triplet<double>(nodeID[6],particle->id,1.0/8.0*(1-r)*(1+s)*(1+t)));
            body->SipTriplets.emplace_back(Eigen::Triplet<double>(nodeID[7],particle->id,1.0/8.0*(1+r)*(1+s)*(1+t)));

            body->gradSipXTriplets.emplace_back(Eigen::Triplet<double>(nodeID[0],particle->id , -dx*dx*dx/8.0*(1-r)*(1-s)*(1-t)*(1-sp)*(1-tp)));
            body->gradSipXTriplets.emplace_back(Eigen::Triplet<double>(nodeID[1],particle->id , dx*dx*dx/8.0*(1+r)*(1-s)*(1-t)*(1-sp)*(1-tp)));
            body->gradSipXTriplets.emplace_back(Eigen::Triplet<double>(nodeID[2],particle->id , -dx*dx*dx/8.0*(1-r)*(1+s)*(1-t)*(1+sp)*(1-tp)));
            body->gradSipXTriplets.emplace_back(Eigen::Triplet<double>(nodeID[3],particle->id , dx*dx*dx/8.0*(1+r)*(1+s)*(1-t)*(1+sp)*(1-tp)));
            body->gradSipXTriplets.emplace_back(Eigen::Triplet<double>(nodeID[4],particle->id , -dx*dx*dx/8.0*(1-r)*(1-s)*(1+t)*(1-sp)*(1+tp)));
            body->gradSipXTriplets.emplace_back(Eigen::Triplet<double>(nodeID[5],particle->id , dx*dx*dx/8.0*(1+r)*(1-s)*(1+t)*(1-sp)*(1+tp)));
            body->gradSipXTriplets.emplace_back(Eigen::Triplet<double>(nodeID[6],particle->id , -dx*dx*dx/8.0*(1-r)*(1+s)*(1+t)*(1+sp)*(1+tp)));
            body->gradSipXTriplets.emplace_back(Eigen::Triplet<double>(nodeID[7],particle->id , dx*dx*dx/8.0*(1+r)*(1+s)*(1+t)*(1+sp)*(1+tp)));

            body->gradSipYTriplets.emplace_back(Eigen::Triplet<double>(nodeID[0],particle->id , dx*dx*dx/8.0*(1-r)*(1-s)*(1-t)*(1-rp)*(1-tp)));
            body->gradSipYTriplets.emplace_back(Eigen::Triplet<double>(nodeID[1],particle->id , dx*dx*dx/8.0*(1+r)*(1-s)*(1-t)*(1+rp)*(1-tp)));
            body->gradSipYTriplets.emplace_back(Eigen::Triplet<double>(nodeID[2],particle->id , -dx*dx*dx/8.0*(1-r)*(1+s)*(1-t)*(1-rp)*(1-tp)));
            body->gradSipYTriplets.emplace_back(Eigen::Triplet<double>(nodeID[3],particle->id , -dx*dx*dx/8.0*(1+r)*(1+s)*(1-t)*(1+rp)*(1-tp)));
            body->gradSipYTriplets.emplace_back(Eigen::Triplet<double>(nodeID[4],particle->id , dx*dx*dx/8.0*(1-r)*(1-s)*(1+t)*(1-rp)*(1+tp)));
            body->gradSipYTriplets.emplace_back(Eigen::Triplet<double>(nodeID[5],particle->id , dx*dx*dx/8.0*(1+r)*(1-s)*(1+t)*(1+rp)*(1+tp)));
            body->gradSipYTriplets.emplace_back(Eigen::Triplet<double>(nodeID[6],particle->id , -dx*dx*dx/8.0*(1-r)*(1+s)*(1+t)*(1-rp)*(1+tp)));
            body->gradSipYTriplets.emplace_back(Eigen::Triplet<double>(nodeID[7],particle->id , -dx*dx*dx/8.0*(1+r)*(1+s)*(1+t)*(1+rp)*(1+tp)));

            body->gradSipZTriplets.emplace_back(Eigen::Triplet<double>(nodeID[0],particle->id , dx*dx*dx/8.0*(1-r)*(1-s)*(1-t)*(1-sp)*(1-rp)));
            body->gradSipZTriplets.emplace_back(Eigen::Triplet<double>(nodeID[1],particle->id , dx*dx*dx/8.0*(1+r)*(1-s)*(1-t)*(1-sp)*(1+rp)));
            body->gradSipZTriplets.emplace_back(Eigen::Triplet<double>(nodeID[2],particle->id , dx*dx*dx/8.0*(1-r)*(1+s)*(1-t)*(1+sp)*(1-rp)));
            body->gradSipZTriplets.emplace_back(Eigen::Triplet<double>(nodeID[3],particle->id , dx*dx*dx/8.0*(1+r)*(1+s)*(1-t)*(1+sp)*(1+rp)));
            body->gradSipZTriplets.emplace_back(Eigen::Triplet<double>(nodeID[4],particle->id , -dx*dx*dx/8.0*(1-r)*(1-s)*(1+t)*(1-sp)*(1-rp)));
            body->gradSipZTriplets.emplace_back(Eigen::Triplet<double>(nodeID[5],particle->id , -dx*dx*dx/8.0*(1+r)*(1-s)*(1+t)*(1-sp)*(1+rp)));
            body->gradSipZTriplets.emplace_back(Eigen::Triplet<double>(nodeID[6],particle->id , -dx*dx*dx/8.0*(1-r)*(1+s)*(1+t)*(1+sp)*(1-rp)));
            body->gradSipZTriplets.emplace_back(Eigen::Triplet<double>(nodeID[7],particle->id , -dx*dx*dx/8.0*(1+r)*(1+s)*(1+t)*(1+sp)*(1+rp)));

            /*body->Sip.coeffRef(nodeID[0],particle->id) += 1.0/8.0*(1-r)*(1-s)*(1-t);
            body->Sip.coeffRef(nodeID[1],particle->id) += 1.0/8.0*(1+r)*(1-s)*(1-t);
            body->Sip.coeffRef(nodeID[2],particle->id) += 1.0/8.0*(1-r)*(1+s)*(1-t);
            body->Sip.coeffRef(nodeID[3],particle->id) += 1.0/8.0*(1+r)*(1+s)*(1-t);
            body->Sip.coeffRef(nodeID[4],particle->id) += 1.0/8.0*(1-r)*(1-s)*(1+t);
            body->Sip.coeffRef(nodeID[5],particle->id) += 1.0/8.0*(1+r)*(1-s)*(1+t);
            body->Sip.coeffRef(nodeID[6],particle->id) += 1.0/8.0*(1-r)*(1+s)*(1+t);
            body->Sip.coeffRef(nodeID[7],particle->id) += 1.0/8.0*(1+r)*(1+s)*(1+t);

            body->gradSipX.coeffRef(nodeID[0],particle->id) += -dx*dx*dx/8.0*(1-r)*(1-s)*(1-t)*(1-sp)*(1-tp);
            body->gradSipX.coeffRef(nodeID[1],particle->id) += dx*dx*dx/8.0*(1+r)*(1-s)*(1-t)*(1-sp)*(1-tp);
            body->gradSipX.coeffRef(nodeID[2],particle->id) += -dx*dx*dx/8.0*(1-r)*(1+s)*(1-t)*(1+sp)*(1-tp);
            body->gradSipX.coeffRef(nodeID[3],particle->id) += dx*dx*dx/8.0*(1+r)*(1+s)*(1-t)*(1+sp)*(1-tp);
            body->gradSipX.coeffRef(nodeID[4],particle->id) += -dx*dx*dx/8.0*(1-r)*(1-s)*(1+t)*(1-sp)*(1+tp);
            body->gradSipX.coeffRef(nodeID[5],particle->id) += dx*dx*dx/8.0*(1+r)*(1-s)*(1+t)*(1-sp)*(1+tp);
            body->gradSipX.coeffRef(nodeID[6],particle->id) += -dx*dx*dx/8.0*(1-r)*(1+s)*(1+t)*(1+sp)*(1+tp);
            body->gradSipX.coeffRef(nodeID[7],particle->id) += dx*dx*dx/8.0*(1+r)*(1+s)*(1+t)*(1+sp)*(1+tp);

            body->gradSipY.coeffRef(nodeID[0],particle->id) += dx*dx*dx/8.0*(1-r)*(1-s)*(1-t)*(1-rp)*(1-tp);
            body->gradSipY.coeffRef(nodeID[1],particle->id) += dx*dx*dx/8.0*(1+r)*(1-s)*(1-t)*(1+rp)*(1-tp);
            body->gradSipY.coeffRef(nodeID[2],particle->id) += -dx*dx*dx/8.0*(1-r)*(1+s)*(1-t)*(1-rp)*(1-tp);
            body->gradSipY.coeffRef(nodeID[3],particle->id) += -dx*dx*dx/8.0*(1+r)*(1+s)*(1-t)*(1+rp)*(1-tp);
            body->gradSipY.coeffRef(nodeID[4],particle->id) += dx*dx*dx/8.0*(1-r)*(1-s)*(1+t)*(1-rp)*(1+tp);
            body->gradSipY.coeffRef(nodeID[5],particle->id) += dx*dx*dx/8.0*(1+r)*(1-s)*(1+t)*(1+rp)*(1+tp);
            body->gradSipY.coeffRef(nodeID[6],particle->id) += -dx*dx*dx/8.0*(1-r)*(1+s)*(1+t)*(1-rp)*(1+tp);
            body->gradSipY.coeffRef(nodeID[7],particle->id) += -dx*dx*dx/8.0*(1+r)*(1+s)*(1+t)*(1+rp)*(1+tp);

            body->gradSipZ.coeffRef(nodeID[0],particle->id) += dx*dx*dx/8.0*(1-r)*(1-s)*(1-t)*(1-sp)*(1-rp);
            body->gradSipZ.coeffRef(nodeID[1],particle->id) += dx*dx*dx/8.0*(1+r)*(1-s)*(1-t)*(1-sp)*(1+rp);
            body->gradSipZ.coeffRef(nodeID[2],particle->id) += dx*dx*dx/8.0*(1-r)*(1+s)*(1-t)*(1+sp)*(1-rp);
            body->gradSipZ.coeffRef(nodeID[3],particle->id) += dx*dx*dx/8.0*(1+r)*(1+s)*(1-t)*(1+sp)*(1+rp);
            body->gradSipZ.coeffRef(nodeID[4],particle->id) += -dx*dx*dx/8.0*(1-r)*(1-s)*(1+t)*(1-sp)*(1-rp);
            body->gradSipZ.coeffRef(nodeID[5],particle->id) += -dx*dx*dx/8.0*(1+r)*(1-s)*(1+t)*(1-sp)*(1+rp);
            body->gradSipZ.coeffRef(nodeID[6],particle->id) += -dx*dx*dx/8.0*(1-r)*(1+s)*(1+t)*(1+sp)*(1-rp);
            body->gradSipZ.coeffRef(nodeID[7],particle->id) += -dx*dx*dx/8.0*(1+r)*(1+s)*(1+t)*(1+sp)*(1+rp);*/
        }
    };

};

#endif //MPM_3D_ELEMENT_HPP

//
// Created by aaron on 8/26/16.
// element.hpp
//

#ifndef MPM_3D_ELEMENT_HPP
#define MPM_3D_ELEMENT_HPP

#include <vector>

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

};

#endif //MPM_3D_ELEMENT_HPP

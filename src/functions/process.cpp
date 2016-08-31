//
// Created by aaron on 8/27/16.
// process.cpp
//

#include <iostream>
#include <stdlib.h>

#include "process.hpp"
#include "body.hpp"

//hard coded for now. need to change
job_t::job_t()//:
        //bodies(2,Body(0,0,0))
{ std::cout << "Job created.\n";}


void job_t::createBody(Body *bd, size_t nn, size_t np, size_t ne, size_t id) {
    bd = new Body(nn,np,ne,id);
}

inline void job_t::node_number_to_coords(
        double * x, double * y, double * z,
        size_t node_number,
        size_t Nx, double hx
) {
    size_t i = node_number % Nx;
    size_t j = (node_number/Nx) % Nx;
    size_t k = node_number % (Nx*Nx);

    *x = i*hx;
    *y = j*hx;
    *z = k*hx;
}

inline int job_t::ijkton_safe(int i, int j, int k,
                       int imax, int jmax, int kmax){
    if (i>imax || i<0){
        return -1;
    }
    if (j>jmax || j<0){
        return -1;
    }
    if (k>kmax || k<0){
        return -1;
    }
    return i*jmax*kmax + j*kmax + k;
}

int job_t::importNodesandParticles(const char *nfilename, const char *pfilename){
    //only reads particle file for 2 bodies

    FILE *fp;
    int r;

    size_t numNodes;
    size_t numElements;
    size_t numLinearNodes;
    double Lx;
    double hx;
    size_t numParticles;
    size_t numParticles1;
    size_t numParticles2;
    size_t numBodies;

    char s[16384]; //#from mpm-2d-legacy

    //open grid file and parse
    fp = fopen(nfilename, "r");
    if (fp == NULL){
        std::cout << "Cannot parse grid file: " << nfilename << "\n";
        return -1;
    }

    if (NULL==fgets(s, sizeof(s)/sizeof(char), fp)){
        std::cout << "Cannot parse grid file: " << nfilename << "\n";
        return -1;
    }

    r = sscanf(s,"%i",&numLinearNodes);
    numElements = (numLinearNodes-1)*(numLinearNodes-1)*(numLinearNodes-1); //(number of linear nodes - 1) ^3
    numNodes = numLinearNodes*numLinearNodes*numLinearNodes; //3d domain
    if (r == 1) {
        fgets(s,sizeof(s)/sizeof(char), fp);
        r = sscanf(s,"%g",&Lx);
        hx = Lx/(numLinearNodes-1);
    }
    if (r!=1){
        std::cout << "Cannot parse grid file: " << nfilename << "\n";
        return -1;
    }
    //std::cout << "Number of nodes: " << numNodes << "\n";
    //std::cout << "Number of elements: " << numElements << "\n";

    fclose(fp);

    //open particle file and parse header
    fp = fopen(pfilename, "r");
    if (fp == NULL){
        std::cout << "Cannot parse particle file: " << pfilename << "\n";
        return -1;
    }

    if (NULL==fgets(s, sizeof(s)/sizeof(char), fp)){
        std::cout << "Cannot parse particle file: " << pfilename << "\n";
        return -1;
    }

    r = sscanf(s,"%i %i %i",&numParticles,&numParticles1,&numParticles2);
    if (r!=3){
        std::cout << "Cannot parse particle file: " << pfilename << "\n";
        std::cout << "Expected 3 particle counts at file header." << "\n";
        std::cout << "Got: " << r << "\n";
        return -1;
    }
    //std::cout << "Number of particles: " << numParticles << "\n";

    //create body objects
    numBodies = 0;
    if (numParticles != 0){
        if (numParticles1 != 0){
            this->bodies.push_back(Body(numNodes,numParticles1,numElements,++numBodies));
            //job_t::createBody(&(this->bodies[0]),numNodes,numParticles1,++numBodies);
        }
        if (numParticles2 != 0){
            this->bodies.push_back(Body(numNodes,numParticles1,numElements,++numBodies));
            //job_t::createBody(&(this->bodies[1]),numNodes,numParticles2,++numBodies);
        }
    } else {
        return -1;
    }
    std::cout << "Bodies created (" << numBodies << ").\n";

    //assign particle to bodies for each particle in particle file
    size_t pb1 = 0;
    size_t pb2 = 0;
    for (size_t i=0; i<numParticles; i++){
        double b, m, v, x, y, z, x_t, y_t, z_t;
        if (NULL==fgets(s, sizeof(s)/sizeof(char), fp)){
            std::cout << "Mismatch! Numbers of particles do not match!\n";
            continue;
        }
        sscanf(s, "%lg %lg %lg %lg %lg %lg %lg %lg %lg", &b, &m, &v, &x, &y, &z, &x_t, &y_t, &z_t);

        //currently assumes that bodies are 1-indexed
        this->bodies[b-1].addParticle(m,v,x,y,z,x_t,y_t,z_t,pb1); //0-index particles
        pb1+=1;
    }
    //close file
    fclose(fp);

    std::cout << "Particles created (" << numParticles << ").\n";

    //assign nodes to bodies
    for (size_t i=0; i<numBodies; i++){
        for (size_t nn=0; nn<numNodes; nn++){
            double x, y, z;
            job_t::node_number_to_coords(&x,&y,&z,nn,numLinearNodes,hx);
            this->bodies[i].addNode(x,y,z,nn);
        }

    }

    std::cout << "Nodes created (" << numNodes << ").\n";

    //assign elements to bodies and nodes to elements
    for (size_t i=0; i<numBodies;i++){
        for(size_t ne=0; ne<numElements; ne++){
            size_t nodeIDs[8];
            size_t c = ne % (numLinearNodes-1);
            size_t r = (ne/(numLinearNodes-1)) % (numLinearNodes-1);
            size_t l = ne / ((numLinearNodes-1)*(numLinearNodes-1));
            size_t n = ijkton_safe(c,r,l,numLinearNodes-1,numLinearNodes-1,numLinearNodes-1);

            nodeIDs[0] = n + ijkton_safe(0,0,0,numLinearNodes-1,numLinearNodes-1,numLinearNodes-1);
            nodeIDs[1] = n + ijkton_safe(0,0,1,numLinearNodes-1,numLinearNodes-1,numLinearNodes-1);
            nodeIDs[2] = n + ijkton_safe(0,1,0,numLinearNodes-1,numLinearNodes-1,numLinearNodes-1);
            nodeIDs[3] = n + ijkton_safe(0,1,1,numLinearNodes-1,numLinearNodes-1,numLinearNodes-1);
            nodeIDs[4] = n + ijkton_safe(1,0,0,numLinearNodes-1,numLinearNodes-1,numLinearNodes-1);
            nodeIDs[5] = n + ijkton_safe(1,0,1,numLinearNodes-1,numLinearNodes-1,numLinearNodes-1);
            nodeIDs[6] = n + ijkton_safe(1,1,0,numLinearNodes-1,numLinearNodes-1,numLinearNodes-1);
            nodeIDs[7] = n + ijkton_safe(1,1,1,numLinearNodes-1,numLinearNodes-1,numLinearNodes-1);

            this->bodies[i].addElement(nodeIDs,ne);
        }
    }

    std::cout << "Elements created (" << numElements << ").\n";
}
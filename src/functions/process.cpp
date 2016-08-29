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


void job_t::createBody(Body *bd, size_t nn, size_t np, size_t id) {
    bd = new Body(nn,np,id);
}

int job_t::importNodesandParticles(const char *nfilename, const char *pfilename){
    //only reads particle file for 2 bodies

    FILE *fp;
    int r;

    size_t numNodes;
    double Lx;
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

    r = sscanf(s,"%i",&numNodes);
    numNodes = numNodes*numNodes*numNodes; //3d domain
    if (r == 1) {
        fgets(s,sizeof(s)/sizeof(char), fp);
        r = sscanf(s,"%g",&Lx);
    }
    if (r!=1){
        std::cout << "Cannot parse grid file: " << nfilename << "\n";
        return -1;
    }
    std::cout << "Number of nodes: " << numNodes << "\n";

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
    std::cout << "Number of particles: " << numParticles << "\n";

    //create body objects
    numBodies = 0;
    if (numParticles != 0){
        if (numParticles1 != 0){
            this->bodies.push_back(Body(numNodes,numParticles1,++numBodies));
            //job_t::createBody(&(this->bodies[0]),numNodes,numParticles1,++numBodies);
        }
        if (numParticles2 != 0){
            this->bodies.push_back(Body(numNodes,numParticles1,++numBodies));
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
    std::cout << "Particles created (" << numParticles << ").\n";

    fclose(fp);
}
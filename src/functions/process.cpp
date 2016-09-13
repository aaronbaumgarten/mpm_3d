//
// Created by aaron on 8/27/16.
// process.cpp
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "process.hpp"
#include "body.hpp"


//hard coded for now. need to change
job_t::job_t():
        use_cpdi(1),
        dt(1e-6)
{ std::cout << "Job created.\n";}


/*void job_t::createBody(Body *bd, size_t nn, size_t np, size_t ne, size_t id) {
    bd = new Body(nn,np,ne,id);
}*/

inline void job_t::node_number_to_coords(
        double * x, double * y, double * z,
        size_t node_number,
        size_t Nx, double hx
) {
    size_t i = node_number % Nx;
    size_t j = (node_number/Nx) % Nx;
    size_t k = node_number / (Nx*Nx);

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
    //return i*jmax*kmax + j*kmax + k;
    return i + j*imax + k*jmax*imax;
}

int job_t::importNodesandParticles(const char *nfilename, const char *pfilename){
    //only reads particle file for 2 bodies

    FILE *fp;
    int r = 0;

    size_t numNodes;
    size_t numElements;
    size_t numLinearNodes;
    //double Lx;
    //double hx;
    size_t numParticles;
    size_t numParticles1;
    size_t numParticles2;
    size_t numBodies;

    char s[16384]; //#from mpm-2d-legacy

    //open grid file and parse
    std::ifstream fin(nfilename);
    if (fin.is_open()){
        std::string line;
        if (std::getline(fin,line)){
            std::stringstream ss(line);
            if (!(ss >> numLinearNodes)){
                std::cout << "Cannot parse grid file: " << nfilename << "\n";
                return -1;
            }
            this->Nx = numLinearNodes;
            this->Ny = numLinearNodes;
            this->Nz = numLinearNodes;
            numElements = (numLinearNodes-1)*(numLinearNodes-1)*(numLinearNodes-1); //(number of linear nodes - 1) ^3
            numNodes = numLinearNodes*numLinearNodes*numLinearNodes; //3d domain
            this->num_nodes = numNodes;
            this->num_elements = numElements;
        }
        if (std::getline(fin,line)) {
            std::stringstream ss(line);
            if (!(ss >> this->Lx)) {
                std::cout << "Cannot parse grid file: " << nfilename << "\n";
                return -1;
            }
            this->Ly = this->Lx;
            this->Lz = this->Lx;
            this->hx = this->Lx / (numLinearNodes - 1);
            this->hy = this->hx;
            this->hz = this->hx;
        }
    } else {
        std::cout << "Cannot parse grid file: " << nfilename << "\n";
        return -1;
    }

    fin.close();

    /*//open grid file and parse
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
    num_nodes = numNodes;
    num_elements = numElements;
    if (r == 1) {
        if (NULL==fgets(s,sizeof(s)/sizeof(char), fp)){
            std::cout << "Cannot parse grid file: " << nfilename << "\n";
            return -1;
        }
        r = sscanf(s,"%g",&Lx);
        hx = Lx/(numLinearNodes-1);
    }
    if (r!=1){
        std::cout << "Cannot parse grid file: " << nfilename << "\n";
        return -1;
    }
    //std::cout << "Number of nodes: " << numNodes << "\n";
    //std::cout << "Number of elements: " << numElements << "\n";

    fclose(fp);*/

    //open particle file
    fin.open(pfilename);
    if (fin.is_open()){
        std::string line;

        //parse header
        if (getline(fin,line)){
            std::stringstream ss(line);
            if (!(ss >> numParticles >> numParticles1 >> numParticles2)){
                std::cout << "Cannot parse particle file: " << pfilename << "\n";
                std::cout << "Expected 3 particle counts at file header." << "\n";
                return -1;
            }
        }
        //create body objects
        numBodies = 0;
        if (numParticles != 0){
            if (numParticles1 != 0){
                this->bodies.push_back(Body(numNodes,numParticles1,numElements,++numBodies));
            }
            if (numParticles2 != 0){
                this->bodies.push_back(Body(numNodes,numParticles1,numElements,++numBodies));
            }
        } else {
            return -1;
        }
        this->num_particles = numParticles;
        this->num_bodies = numBodies;
        std::cout << "Bodies created (" << numBodies << ").\n";

        //assign particle to bodies for each particle in particle file
        size_t pb1 = 0;
        size_t pb2 = 0;
        while(getline(fin,line)){
            std::stringstream ss(line);
            double b, m, v, x, y, z, x_t, y_t, z_t;
            size_t idOut;

            if (!(ss >> b >> m >> v >> x >> y >> z >> x_t >> y_t >> z_t)){
                std::cout << "Cannot parse particle file: " << pfilename << "\n";
                return -1;
            }
            //currently assumes that bodies are 1-indexed
            if (b==1){
                idOut = pb1;
                pb1+=1;
            } else {
                idOut = pb2;
                pb2+=1;
            }
            this->bodies[b-1].addParticle(m,v,x,y,z,x_t,y_t,z_t,idOut); //0-index particles
        }
    } else {
        std::cout << "Cannot parse particle file: " << pfilename << "\n";
        return -1;
    }
    fin.close();

    /*//open particle file and parse header
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
    num_particles = numParticles;
    if (r!=3){
        std::cout << "Cannot parse particle file: " << pfilename << "\n";
        std::cout << "Expected 3 particle counts at file header." << "\n";
        std::cout << "Got: z" << r << "\n";
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
    num_bodies = numBodies;
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
    fclose(fp);*/

    std::cout << "Particles created (" << numParticles << ").\n";

    //assign nodes to bodies
    for (size_t i=0; i<numBodies; i++){
        for (size_t nn=0; nn<numNodes; nn++){
            double x, y, z;
            job_t::node_number_to_coords(&x,&y,&z,nn,numLinearNodes,this->hx);
            this->bodies[i].addNode(x,y,z,nn);
        }

    }

    std::cout << "Nodes created (" << numNodes << ").\n";

    //assign elements to bodies and nodes to elements
    for (size_t i=0; i<numBodies;i++){
        for(size_t ne=0; ne<numElements; ne++){
            size_t nodeIDs[8];
            size_t c = ne % (numLinearNodes);
            size_t r = (ne/(numLinearNodes)) % (numLinearNodes);
            size_t l = ne / ((numLinearNodes)*(numLinearNodes));
            size_t n = ijkton_safe(c,r,l,numLinearNodes,numLinearNodes,numLinearNodes);

            nodeIDs[0] = n + ijkton_safe(0,0,0,numLinearNodes,numLinearNodes,numLinearNodes);
            nodeIDs[1] = n + ijkton_safe(1,0,0,numLinearNodes,numLinearNodes,numLinearNodes);
            nodeIDs[2] = n + ijkton_safe(0,1,0,numLinearNodes,numLinearNodes,numLinearNodes);
            nodeIDs[3] = n + ijkton_safe(1,1,0,numLinearNodes,numLinearNodes,numLinearNodes);
            nodeIDs[4] = n + ijkton_safe(0,0,1,numLinearNodes,numLinearNodes,numLinearNodes);
            nodeIDs[5] = n + ijkton_safe(1,0,1,numLinearNodes,numLinearNodes,numLinearNodes);
            nodeIDs[6] = n + ijkton_safe(0,1,1,numLinearNodes,numLinearNodes,numLinearNodes);
            nodeIDs[7] = n + ijkton_safe(1,1,1,numLinearNodes,numLinearNodes,numLinearNodes);

            this->bodies[i].addElement(nodeIDs,ne);
        }
    }

    std::cout << "Elements created (" << numElements << ").\n";
    return 1;
}

int job_t::assignMaterials() {
    for (size_t i=0;i<num_bodies;i++){
        if (i==0) {
            this->bodies[i].material.calculate_stress = material1::calculate_stress;
            this->bodies[i].material.calculate_stress_threaded = material1::calculate_stress_threaded;
            this->bodies[i].material.material_init = material1::material_init;
        } else {
            this->bodies[i].material.calculate_stress = material2::calculate_stress;
            this->bodies[i].material.calculate_stress_threaded = material2::calculate_stress_threaded;
            this->bodies[i].material.material_init = material2::material_init;
        }
    }
    std::cout << "Materials assigned (" << num_bodies << ").\n";
    return 1;
}

int job_t::assignMaterials(const char* matFile1, const char* matFile2){
    // for defining material based on .so files (not implemented as of 9/8/16
    return -1;
}
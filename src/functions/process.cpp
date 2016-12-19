//
// Created by aaron on 8/27/16.
// process.cpp
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <algorithm>

#include "tensor.hpp"
#include "process.hpp"
#include "body.hpp"
#include "loading.hpp"

//mass tolerance
#define TOL 5e-11
//contact friction
#define MU_F 0.4
//squared norm error tolerance
#define R_TOL 1e-5
#define R_MAX 1e10

//hard coded for now. need to change
job_t::job_t():
        use_cpdi(1), //default cpdi
        use_3d(1), //default 3d
        use_implicit(0), //default explicit
        dt(1e-3),
        dt_base(dt),
        dt_minimum(1e-6),
        t(0.0),
        step_start_time(0.0),
        stepcount(0),
        newtonTOL(1e-5),
        linearStepSize(1e-6)
{
    std::cout << "Job created.\n";
    boundary = Boundary();
}


/*void job_t::createBody(Body *bd, size_t nn, size_t np, size_t ne, size_t id) {
    bd = new Body(nn,np,ne,id);
}*/

inline void job_t::node_number_to_coords(
        double * x, double * y, double * z,
        size_t node_number,
        size_t Nx, size_t Ny, size_t Nz,
        double hx, double hy, double hz
) {
    size_t i = node_number % Nx;
    size_t j = (node_number/Nx) % Ny;
    size_t k = node_number / (Nx*Ny);

    *x = i*hx;
    *y = j*hy;
    *z = k*hz;
    return;
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

//int job_t::importNodesandParticles(const char *nfilename, const char *pfilename){
int job_t::importNodesandParticles(std::string nfilename, std::string pfilename){
    //only reads particle file for 2 bodies

    if (this->use_3d == 1) {
        FILE *fp;

        size_t numNodes = 0;
        size_t numElements = 0;
        size_t numLinearNodesX = 0;
        size_t numLinearNodesY = 0;
        size_t numLinearNodesZ = 0;
        //double Lx;
        //double hx;
        size_t numParticles = 0;
        size_t numParticles1 = 0;
        size_t numParticles2 = 0;
        size_t numBodies = 0;

        char s[16384]; //#from mpm-2d-legacy

        //open grid file and parse
        std::ifstream fin(nfilename);
        if (fin.is_open()) {
            std::string line;
            if (std::getline(fin, line)) {
                std::stringstream ss(line);
                if (!(ss >> numLinearNodesX >> numLinearNodesY >> numLinearNodesZ)) {
                    std::cout << "Cannot parse grid file: " << nfilename << "\n";
                    return -1;
                }
                this->Nx = numLinearNodesX;
                this->Ny = numLinearNodesY;
                this->Nz = numLinearNodesZ;
                numElements = (numLinearNodesX - 1) * (numLinearNodesY - 1) *
                              (numLinearNodesZ - 1); //(number of linear nodes - 1) ^3
                numNodes = numLinearNodesX * numLinearNodesY * numLinearNodesZ; //3d domain
                this->num_nodes = numNodes;
                this->num_elements = numElements;

                this->u_dirichlet.resize(numNodes * NODAL_DOF);
                this->u_dirichlet_mask.resize(numNodes * NODAL_DOF);
                this->node_number_override.resize(numNodes * NODAL_DOF);
            }
            if (std::getline(fin, line)) {
                std::stringstream ss(line);
                if (!(ss >> this->Lx)) {
                    std::cout << "Cannot parse grid file: " << nfilename << "\n";
                    return -1;
                }
                this->Ly = this->Lx;
                this->Lz = this->Lx;
                this->hx = this->Lx / (numLinearNodesX - 1);
                this->hy = this->Ly / (numLinearNodesY - 1);
                this->hz = this->Lz / (numLinearNodesZ - 1);
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
        if (fin.is_open()) {
            std::string line;

            //parse header
            if (getline(fin, line)) {
                std::stringstream ss(line);
                if (!(ss >> numParticles >> numParticles1 >> numParticles2)) {
                    std::cout << "Cannot parse particle file: " << pfilename << "\n";
                    std::cout << "Expected 3 particle counts at file header." << "\n";
                    return -1;
                }
            }
            //create body objects
            numBodies = 0;
            if (numParticles != 0) {
                if (numParticles1 != 0) {
                    this->bodies.push_back(Body(numNodes, numParticles1, numElements, ++numBodies));
                }
                if (numParticles2 != 0) {
                    this->bodies.push_back(Body(numNodes, numParticles1, numElements, ++numBodies));
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
            while (getline(fin, line)) {
                std::stringstream ss(line);
                double b, m, v, x, y, z, x_t, y_t, z_t;
                size_t idOut;

                if (!(ss >> b >> m >> v >> x >> y >> z >> x_t >> y_t >> z_t)) {
                    std::cout << "Cannot parse particle file: " << pfilename << "\n";
                    return -1;
                }
                //currently assumes that bodies are 1-indexed
                if (b == 1) {
                    idOut = pb1;
                    pb1 += 1;
                } else {
                    idOut = pb2;
                    pb2 += 1;
                }
                //std::cout << "{" << x << " " << y << " " << z << "} -> " << idOut << "\n";
                this->bodies[b - 1].addParticle(m, v, x, y, z, x_t, y_t, z_t, idOut); //0-index particles
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
        for (size_t i = 0; i < numBodies; i++) {
            for (size_t nn = 0; nn < numNodes; nn++) {
                double x, y, z;
                job_t::node_number_to_coords(&x, &y, &z, nn, numLinearNodesX, numLinearNodesY, numLinearNodesZ, this->hx, this->hy, this->hz);
                this->bodies[i].addNode(x, y, z, nn);
            }

        }

        std::cout << "Nodes created (" << numNodes << ").\n";

        //assign elements to bodies and nodes to elements
        for (size_t i = 0; i < numBodies; i++) {
            for (size_t ne = 0; ne < numElements; ne++) {
                size_t nodeIDs[8];
                size_t c = ne % (numLinearNodesX - 1);
                size_t r = (ne / (numLinearNodesX - 1)) % (numLinearNodesY - 1);
                size_t l = ne / ((numLinearNodesX - 1) * (numLinearNodesY - 1));
                size_t n = ijkton_safe(c, r, l, numLinearNodesX, numLinearNodesY, numLinearNodesZ);

                nodeIDs[0] = n + ijkton_safe(0, 0, 0, numLinearNodesX, numLinearNodesY, numLinearNodesZ);
                nodeIDs[1] = n + ijkton_safe(1, 0, 0, numLinearNodesX, numLinearNodesY, numLinearNodesZ);
                nodeIDs[2] = n + ijkton_safe(0, 1, 0, numLinearNodesX, numLinearNodesY, numLinearNodesZ);
                nodeIDs[3] = n + ijkton_safe(1, 1, 0, numLinearNodesX, numLinearNodesY, numLinearNodesZ);
                nodeIDs[4] = n + ijkton_safe(0, 0, 1, numLinearNodesX, numLinearNodesY, numLinearNodesZ);
                nodeIDs[5] = n + ijkton_safe(1, 0, 1, numLinearNodesX, numLinearNodesY, numLinearNodesZ);
                nodeIDs[6] = n + ijkton_safe(0, 1, 1, numLinearNodesX, numLinearNodesY, numLinearNodesZ);
                nodeIDs[7] = n + ijkton_safe(1, 1, 1, numLinearNodesX, numLinearNodesY, numLinearNodesZ);

                //check uniqueness of element
                //for(size_t c1=0;c1<8;c1++){
                //    for(size_t c2;c2<8;c2++){
                //        if(c2!=c1 && nodeIDs[c1]==nodeIDs[c2]){
                //            std::cout << "Error: node " << c1 << " and node " << c2 << " identical in element " << ne << ".\n";
                //        }
                //    }
                //}

                this->bodies[i].addElement(nodeIDs, ne);
            }
        }

        std::cout << "Elements created (" << numElements << ").\n";
        return 1;
    } else {
        return this->importNodesandParticles2D(nfilename,pfilename);
    }
}

int job_t::importNodesandParticles2D(std::string nfilename, std::string pfilename){
    //only reads particle file for 2 bodies

    std::cout << "Setting up 2D simulation.\n";

    FILE *fp;

    size_t numNodes = 0;
    size_t numElements = 0;
    size_t numLinearNodesX = 0;
    size_t numLinearNodesY = 0;
    size_t numLinearNodesZ = 0;
    //double Lx;
    //double hx;
    size_t numParticles = 0;
    size_t numParticles1 = 0;
    size_t numParticles2 = 0;
    size_t numBodies = 0;

    char s[16384]; //#from mpm-2d-legacy

    //open grid file and parse
    std::ifstream fin(nfilename);
    if (fin.is_open()){
        std::string line;
        if (std::getline(fin,line)){
            std::stringstream ss(line);
            if (!(ss >> numLinearNodesX >> numLinearNodesY >> numLinearNodesZ)){
                std::cout << "Cannot parse grid file: " << nfilename << "\n";
                return -1;
            }
            this->Nx = numLinearNodesX;
            this->Ny = numLinearNodesY;
            this->Nz = 2;//numLinearNodes;
            numElements = (numLinearNodesX-1)*(numLinearNodesY-1); //(number of linear nodes - 1) ^2
            numNodes = numLinearNodesX*numLinearNodesY*2; //2d domain
            this->num_nodes = numNodes;
            this->num_elements = numElements;

            this->u_dirichlet.resize(numNodes*NODAL_DOF);
            this->u_dirichlet_mask.resize(numNodes*NODAL_DOF);
            this->node_number_override.resize(numNodes*NODAL_DOF);
        }
        if (std::getline(fin,line)) {
            std::stringstream ss(line);
            if (!(ss >> this->Lx >> this->Ly >> this->Lz)) {
                std::cout << "Cannot parse grid file: " << nfilename << "\n";
                return -1;
            }
            //this->Ly = this->Lx;
            //this->Lz = this->Lx;
            this->hx = this->Lx / (numLinearNodesX - 1);
            this->hy = this->Ly / (numLinearNodesY - 1);
            this->hz = this->Lz / 2;
            //this->Lz = this->hz; //1 element
        }
    } else {
        std::cout << "Cannot parse grid file: " << nfilename << "\n";
        return -1;
    }

    fin.close();

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
            //std::cout << "{" << x << " " << y << " " << z << "} -> " << idOut << "\n";
            this->bodies[b-1].addParticle(m,v,x,y,z,x_t,y_t,z_t,idOut); //0-index particles
        }
    } else {
        std::cout << "Cannot parse particle file: " << pfilename << "\n";
        return -1;
    }
    fin.close();

    std::cout << "Particles created (" << numParticles << ").\n";

    //assign nodes to bodies
    for (size_t i=0; i<numBodies; i++){
        for (size_t nn=0; nn<numNodes; nn++){
            double x, y, z;
            job_t::node_number_to_coords(&x,&y,&z,nn,numLinearNodesX,numLinearNodesY,numLinearNodesZ,this->hx,this->hy,this->hz); //ok for 2d as nodes count in x,y first
            this->bodies[i].addNode(x,y,z,nn);
        }

    }

    std::cout << "Nodes created (" << numNodes << ").\n";

    //assign elements to bodies and nodes to elements
    for (size_t i=0; i<numBodies;i++){
        for(size_t ne=0; ne<numElements; ne++){
            size_t nodeIDs[8];
            size_t c = ne % (this->Nx-1);
            size_t r = (ne/(this->Nx-1)) % (this->Ny-1);
            size_t l = ne / ((this->Nx-1)*(this->Ny-1));
            size_t n = ijkton_safe(c,r,l,this->Nx,this->Ny,this->Nz);

            nodeIDs[0] = n + ijkton_safe(0,0,0,this->Nx,this->Ny,this->Nz);
            nodeIDs[1] = n + ijkton_safe(1,0,0,this->Nx,this->Ny,this->Nz);
            nodeIDs[2] = n + ijkton_safe(0,1,0,this->Nx,this->Ny,this->Nz);
            nodeIDs[3] = n + ijkton_safe(1,1,0,this->Nx,this->Ny,this->Nz);
            nodeIDs[4] = n + ijkton_safe(0,0,1,this->Nx,this->Ny,this->Nz);
            nodeIDs[5] = n + ijkton_safe(1,0,1,this->Nx,this->Ny,this->Nz);
            nodeIDs[6] = n + ijkton_safe(0,1,1,this->Nx,this->Ny,this->Nz);
            nodeIDs[7] = n + ijkton_safe(1,1,1,this->Nx,this->Ny,this->Nz);

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
            this->bodies[i].material.calculate_stress_implicit = material1::calculate_stress_implicit;
        } else {
            this->bodies[i].material.calculate_stress = material2::calculate_stress;
            this->bodies[i].material.calculate_stress_threaded = material2::calculate_stress_threaded;
            this->bodies[i].material.material_init = material2::material_init;
            this->bodies[i].material.calculate_stress_implicit = material2::calculate_stress_implicit;
        }
    }
    std::cout << "Materials assigned (" << num_bodies << ").\n";
    return 1;
}

int job_t::assignMaterials(const char* matFile1, const char* matFile2){
    // for defining material based on .so files (not implemented as of 9/8/16
    return -1;
}


int job_t::assignBoundaryConditions(){
    this->boundary.bc_init = boundary::bc_init;
    this->boundary.bc_validate = boundary::bc_validate;
    this->boundary.bc_time_varying = boundary::bc_time_varying;
    this->boundary.generate_dirichlet_bcs = boundary::generate_dirichlet_bcs;
    this->boundary.generate_node_number_override = boundary::generate_node_number_override;
    this->boundary.bc_momentum = boundary::bc_momentum;
    this->boundary.bc_force = boundary::bc_force;

    std::cout << "Boundary Conditions assigned (1).\n";

    return 1;
}

int job_t::assignBoundaryConditions(const char* bcFile){
    // for defining boundary condition based on .so files (not implemented as of 10/26/16)
    return -1;
}

int job_t::createMappings() {
    //loop over bodies
    for (size_t b = 0; b < this->num_bodies; b++) {
        //reset previous maps
        this->bodies[b].Phi.setZero();
        this->bodies[b].gradPhiX.setZero();
        this->bodies[b].gradPhiY.setZero();
        this->bodies[b].gradPhiZ.setZero();
        //reset triplets
        this->bodies[b].PhiTriplets.clear();
        this->bodies[b].gradPhiXTriplets.clear();
        this->bodies[b].gradPhiYTriplets.clear();
        this->bodies[b].gradPhiZTriplets.clear();
        //loop over particles
        for (size_t p = 0; p < this->bodies[b].p; p++) {
            if((this->bodies[b].particles[p].updateActive(this))==1) {
                this->bodies[b].particles[p].updateCorners(this);
                //check that all corners are in domain
                if (this->use_cpdi==0 || std::count(this->bodies[b].particles[p].corner_elements,
                               std::end(this->bodies[b].particles[p].corner_elements), -1) > 0) {
                    //set corners to particle position
                    this->bodies[b].particles[p].resetCorners(this);
                }

                //use cpdi (or count center 8 times if corners reset)
                for (size_t c=0;c<8;c++) {
                    size_t e = this->bodies[b].particles[p].corner_elements[c];
                    if (e != this->bodies[0].elements[e].id){
                        std::cout << "Elemend ID Error! " << e << " != " << this->bodies[0].elements[e].id << "\n";
                    }
                    this->bodies[b].elements[e].calculatePhic(&(this->bodies[b]), &(this->bodies[b].particles[p]), c, this->use_cpdi);
                    //std::cout << "b:" << b << " p:" << p << " c:" << c << "\r";
                }
            }
        }
        //build Phic from triples
        //std::cout << "complete map\n";
        this->bodies[b].Phi.setFromTriplets(this->bodies[b].PhiTriplets.begin(),this->bodies[b].PhiTriplets.end());
        this->bodies[b].gradPhiX.setFromTriplets(this->bodies[b].gradPhiXTriplets.begin(),this->bodies[b].gradPhiXTriplets.end());
        this->bodies[b].gradPhiY.setFromTriplets(this->bodies[b].gradPhiYTriplets.begin(),this->bodies[b].gradPhiYTriplets.end());
        this->bodies[b].gradPhiZ.setFromTriplets(this->bodies[b].gradPhiZTriplets.begin(),this->bodies[b].gradPhiZTriplets.end());
    }
    return 1;
}


void job_t::mapParticles2Grid() {
    for (size_t b=0; b<this->num_bodies; b++) {
        //use Eigen Map to point to particle array
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;
        
        Eigen::MatrixXd p_m_x_t(numRowsP,numColsP);
        Eigen::MatrixXd p_m_y_t(numRowsP,numColsP);
        Eigen::MatrixXd p_m_z_t(numRowsP,numColsP);
        p_m_x_t = this->bodies[b].particle_m.array()*this->bodies[b].particle_x_t.array();
        p_m_y_t = this->bodies[b].particle_m.array()*this->bodies[b].particle_y_t.array();
        p_m_z_t = this->bodies[b].particle_m.array()*this->bodies[b].particle_z_t.array();

        Eigen::MatrixXd p_m_bx(numRowsP,numColsP);
        Eigen::MatrixXd p_m_by(numRowsP,numColsP);
        Eigen::MatrixXd p_m_bz(numRowsP,numColsP);
        p_m_bx = this->bodies[b].particle_m.array()*this->bodies[b].particle_bx.array();
        p_m_by = this->bodies[b].particle_m.array()*this->bodies[b].particle_by.array();
        p_m_bz = this->bodies[b].particle_m.array()*this->bodies[b].particle_bz.array();

        //use Eigen Map to point to node array
        size_t numRowsN = this->bodies[b].n;
        size_t numColsN = 1;
        
        //use to create dummy pvec and ones
        Eigen::VectorXd pvec(numRowsP);

        //use Phi to map particles to nodes
        this->bodies[b].node_m = this->bodies[b].Phi*this->bodies[b].particle_m;
        this->bodies[b].node_mx_t = this->bodies[b].Phi*p_m_x_t;
        this->bodies[b].node_my_t = this->bodies[b].Phi*p_m_y_t;
        this->bodies[b].node_mz_t = this->bodies[b].Phi*p_m_z_t;
        this->bodies[b].node_fx = this->bodies[b].Phi*p_m_bx; //need to add stress
        this->bodies[b].node_fy = this->bodies[b].Phi*p_m_by; //need to add stress
        this->bodies[b].node_fz = this->bodies[b].Phi*p_m_bz; //need to add stress

        //use gradPhi to map particles to nodes
        this->bodies[b].node_contact_normal_x = this->bodies[b].gradPhiX*Eigen::MatrixXd::Constant(numRowsP,numColsP,1);//this->bodies[b].particle_m;
        this->bodies[b].node_contact_normal_y = this->bodies[b].gradPhiY*Eigen::MatrixXd::Constant(numRowsP,numColsP,1);//this->bodies[b].particle_m;
        this->bodies[b].node_contact_normal_z = this->bodies[b].gradPhiZ*Eigen::MatrixXd::Constant(numRowsP,numColsP,1);//this->bodies[b].particle_m;

        /*Eigen::VectorXd normMag(numRowsN);
        normMag = this->bodies[b].node_contact_normal_x.array().square()
                + this->bodies[b].node_contact_normal_y.array().square()
                + this->bodies[b].node_contact_normal_z.array().square();

        this->bodies[b].node_contact_normal_x = this->bodies[b].node_contact_normal_x.array()/normMag.array().sqrt();
        this->bodies[b].node_contact_normal_y = this->bodies[b].node_contact_normal_y.array()/normMag.array().sqrt();
        this->bodies[b].node_contact_normal_z = this->bodies[b].node_contact_normal_z.array()/normMag.array().sqrt();
         */

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[XX];
        }
        this->bodies[b].node_fx -= this->bodies[b].gradPhiX*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[XY];
        }
        this->bodies[b].node_fx -= this->bodies[b].gradPhiY*pvec;
        this->bodies[b].node_fy -= this->bodies[b].gradPhiX*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[XZ];
        }
        this->bodies[b].node_fx -= this->bodies[b].gradPhiZ*pvec;
        this->bodies[b].node_fz -= this->bodies[b].gradPhiX*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[YY];
        }
        this->bodies[b].node_fy -= this->bodies[b].gradPhiY*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[YZ];
        }
        this->bodies[b].node_fy -= this->bodies[b].gradPhiZ*pvec;
        this->bodies[b].node_fz -= this->bodies[b].gradPhiY*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[ZZ];
        }
        this->bodies[b].node_fz -= this->bodies[b].gradPhiZ*pvec;

    }
    return;
}

void job_t::addContactForces(){
    //if (this->use_3d==1) {
    for (size_t b = 0; b < this->num_bodies; b++) {
        this->bodies[b].node_contact_mx_t = this->bodies[b].node_mx_t;
        this->bodies[b].node_contact_my_t = this->bodies[b].node_my_t;
        this->bodies[b].node_contact_mz_t = this->bodies[b].node_mz_t;

        //the following appear unused
        this->bodies[b].node_contact_x_t = this->bodies[b].node_mx_t.array()/this->bodies[b].node_m.array();
        this->bodies[b].node_contact_y_t = this->bodies[b].node_my_t.array()/this->bodies[b].node_m.array();
        this->bodies[b].node_contact_z_t = this->bodies[b].node_mz_t.array()/this->bodies[b].node_m.array();

        this->bodies[b].node_contact_fx = this->bodies[b].node_fx;
        this->bodies[b].node_contact_fy = this->bodies[b].node_fy;
        this->bodies[b].node_contact_fz = this->bodies[b].node_fz;
    }

    if (this->num_bodies > 1) {
        //look for contacts if there are two bodies
        for (size_t i = 0; i < this->num_nodes; i++) {
            //test every node for contact
            if (this->bodies[0].node_m[i] > TOL && this->bodies[1].node_m[i] > TOL) {
                //use normal from body 1
                Eigen::Vector3d n1i;
                n1i << this->bodies[0].node_contact_normal_x[i],
                        this->bodies[0].node_contact_normal_y[i],
                        this->bodies[0].node_contact_normal_z[i];
                //enforce unit length
                n1i /= sqrt(n1i.dot(n1i));

                //determine 'center of mass' velocity
                double m1 = this->bodies[0].node_m[i];
                double m2 = this->bodies[1].node_m[i];
                Eigen::Vector3d mv1i;
                Eigen::Vector3d mv2i;
                Eigen::Vector3d vCMi;
                mv1i << this->bodies[0].node_contact_mx_t[i],
                        this->bodies[0].node_contact_my_t[i],
                        this->bodies[0].node_contact_mz_t[i];
                mv2i << this->bodies[1].node_contact_mx_t[i],
                        this->bodies[1].node_contact_my_t[i],
                        this->bodies[1].node_contact_mz_t[i];
                vCMi = (mv1i + mv2i) / (m1 + m2);

                //determine normal force
                double fn1i;
                fn1i = m1 * m2 / (this->dt * (m1 + m2)) * (mv2i.dot(n1i) / m2 - mv1i.dot(n1i) / m1);

                //determine shear force and shear vector
                double ft1i;
                Eigen::Vector3d s1i;
                s1i = m1 / this->dt * (vCMi - mv1i / m1) - fn1i * n1i;
                ft1i = sqrt(s1i.dot(s1i));
                s1i /= ft1i;

                //add forces
                Eigen::Vector3d fcti;
                fcti = std::min(0.0, fn1i)*n1i + std::min(MU_F*std::abs(fn1i),std::abs(ft1i))*s1i;

                //set contact forces
                this->bodies[0].node_contact_fx[i] = fcti[0];
                this->bodies[0].node_contact_fy[i] = fcti[1];
                this->bodies[0].node_contact_fz[i] = fcti[2];

                this->bodies[1].node_contact_fx[i] = -fcti[0];
                this->bodies[1].node_contact_fy[i] = -fcti[1];
                this->bodies[1].node_contact_fz[i] = -fcti[2];

                //adjust nodal velocities for non-penetration
                mv1i = mv1i - n1i.dot(mv1i-m1*vCMi)*n1i;
                mv2i = mv2i - n1i.dot(mv2i-m2*vCMi)*n1i;

                this->bodies[0].node_contact_mx_t[i] = mv1i[0];
                this->bodies[0].node_contact_my_t[i] = mv1i[1];
                this->bodies[0].node_contact_mz_t[i] = mv1i[2];

                this->bodies[0].node_contact_x_t[i] = mv1i[0]/m1;
                this->bodies[0].node_contact_y_t[i] = mv1i[1]/m1;
                this->bodies[0].node_contact_z_t[i] = mv1i[2]/m1;

                this->bodies[1].node_contact_mx_t[i] = mv2i[0];
                this->bodies[1].node_contact_my_t[i] = mv2i[1];
                this->bodies[1].node_contact_mz_t[i] = mv2i[2];

                this->bodies[1].node_contact_x_t[i] = mv2i[0]/m2;
                this->bodies[1].node_contact_y_t[i] = mv2i[1]/m2;
                this->bodies[1].node_contact_z_t[i] = mv2i[2]/m2;
            }
        }
    }

    return;
    //} else {
    //    return this->addContactForces2D();
    //}
}

void job_t::addContactForces2D(){
    //resolve conflicts between grid velocites
    //implement later 10/21/16
    for (size_t b=0;b<this->num_bodies;b++) {
        this->bodies[b].node_contact_mx_t = this->bodies[b].node_mx_t;
        this->bodies[b].node_contact_my_t = this->bodies[b].node_my_t;
        this->bodies[b].node_contact_mz_t = this->bodies[b].node_mz_t;

        //the following appear unused
        this->bodies[b].node_contact_x_t = this->bodies[b].node_mx_t.array()/this->bodies[b].node_m.array();
        this->bodies[b].node_contact_y_t = this->bodies[b].node_my_t.array()/this->bodies[b].node_m.array();
        this->bodies[b].node_contact_z_t = this->bodies[b].node_mz_t.array()/this->bodies[b].node_m.array();

        this->bodies[b].node_contact_fx = this->bodies[b].node_fx;
        this->bodies[b].node_contact_fy = this->bodies[b].node_fy;
        this->bodies[b].node_contact_fz = this->bodies[b].node_fz;
    }

    if (this->num_bodies > 1) {
        //look for contacts if there are two bodies
        for (size_t i = 0; i < this->num_nodes; i++) {
            //test every node for contact
            if (this->bodies[0].node_m[i] > TOL && this->bodies[1].node_m[i] > TOL) {
                //use normal from body 1
                Eigen::Vector2d n1i;
                n1i << this->bodies[0].node_contact_normal_x[i],
                        this->bodies[0].node_contact_normal_y[i];
                //enforce unit length
                n1i /= sqrt(n1i.dot(n1i));

                //determine 'center of mass' velocity
                double m1 = this->bodies[0].node_m[i];
                double m2 = this->bodies[1].node_m[i];
                Eigen::Vector2d mv1i;
                Eigen::Vector2d mv2i;
                Eigen::Vector2d vCMi;
                mv1i << this->bodies[0].node_contact_mx_t[i],
                        this->bodies[0].node_contact_my_t[i];
                mv2i << this->bodies[1].node_contact_mx_t[i],
                        this->bodies[1].node_contact_my_t[i];
                vCMi = (mv1i + mv2i) / (m1 + m2);

                //determine normal force
                double fn1i;
                fn1i = m1 * m2 / (this->dt * (m1 + m2)) * (mv2i.dot(n1i) / m2 - mv1i.dot(n1i) / m1);

                //determine shear force and shear vector
                double ft1i;
                Eigen::Vector2d s1i;
                s1i = m1 / this->dt * (vCMi - mv1i / m1) - fn1i * n1i;
                ft1i = sqrt(s1i.dot(s1i));
                s1i /= ft1i;

                //add forces
                Eigen::Vector2d fcti;
                fcti = std::min(0.0, fn1i)*n1i + std::min(MU_F*std::abs(fn1i),std::abs(ft1i))*s1i;

                //set contact forces
                this->bodies[0].node_contact_fx[i] = fcti[0];
                this->bodies[0].node_contact_fy[i] = fcti[1];

                this->bodies[1].node_contact_fx[i] = -fcti[0];
                this->bodies[1].node_contact_fy[i] = -fcti[1];

                //adjust nodal velocities for non-penetration
                mv1i = mv1i - n1i.dot(mv1i-m1*vCMi)*n1i;
                mv2i = mv2i - n1i.dot(mv2i-m2*vCMi)*n1i;

                this->bodies[0].node_contact_mx_t[i] = mv1i[0];
                this->bodies[0].node_contact_my_t[i] = mv1i[1];

                this->bodies[0].node_contact_x_t[i] = mv1i[0]/m1;
                this->bodies[0].node_contact_y_t[i] = mv1i[1]/m1;

                this->bodies[1].node_contact_mx_t[i] = mv2i[0];
                this->bodies[1].node_contact_my_t[i] = mv2i[1];

                this->bodies[1].node_contact_x_t[i] = mv2i[0]/m2;
                this->bodies[1].node_contact_y_t[i] = mv2i[1]/m2;
            }
        }
    }

    return;
}

void job_t::addBoundaryConditions(){
    //initialize and calculate
    this->boundary.bc_time_varying(this);
    this->boundary.bc_momentum(this);
    this->boundary.bc_force(this);
    return;
}

void job_t::moveGridExplicit(){
    //if (this->use_3d==1) {
    for (size_t b = 0; b < this->num_bodies; b++) {
        for (size_t i = 0; i < this->num_nodes; i++) {
            double m = this->bodies[b].nodes[i].m[0];
            if (m > TOL) {
                this->bodies[b].nodes[i].contact_mx_t[0] += this->dt * this->bodies[b].nodes[i].contact_fx[0];
                this->bodies[b].nodes[i].contact_my_t[0] += this->dt * this->bodies[b].nodes[i].contact_fy[0];
                this->bodies[b].nodes[i].contact_mz_t[0] += this->dt * this->bodies[b].nodes[i].contact_fz[0];

                this->bodies[b].nodes[i].contact_x_t[0] = this->bodies[b].nodes[i].contact_mx_t[0] / m;
                this->bodies[b].nodes[i].contact_y_t[0] = this->bodies[b].nodes[i].contact_my_t[0] / m;
                this->bodies[b].nodes[i].contact_z_t[0] = this->bodies[b].nodes[i].contact_mz_t[0] / m;
            } else {
                this->bodies[b].nodes[i].contact_mx_t[0] = 0;
                this->bodies[b].nodes[i].contact_my_t[0] = 0;
                this->bodies[b].nodes[i].contact_mz_t[0] = 0;

                this->bodies[b].nodes[i].contact_x_t[0] = 0;
                this->bodies[b].nodes[i].contact_y_t[0] = 0;
                this->bodies[b].nodes[i].contact_z_t[0] = 0;
            }

        }
    }
    return;
    //} else {
    //    return this->moveGridExplicit2D();
    //}
}

void job_t::moveGridExplicit2D(){
    for (size_t b = 0; b < this->num_bodies; b++){
        for (size_t i = 0; i < this->num_nodes; i++){
            double m = this->bodies[b].nodes[i].m[0];
            if (m > TOL) {
                this->bodies[b].nodes[i].contact_mx_t[0] += this->dt * this->bodies[b].nodes[i].contact_fx[0];
                this->bodies[b].nodes[i].contact_my_t[0] += this->dt * this->bodies[b].nodes[i].contact_fy[0];
                this->bodies[b].nodes[i].contact_mz_t[0] = 0;//this->dt * this->bodies[b].nodes[i].contact_fz[0];

                this->bodies[b].nodes[i].contact_x_t[0] = this->bodies[b].nodes[i].contact_mx_t[0] / m;
                this->bodies[b].nodes[i].contact_y_t[0] = this->bodies[b].nodes[i].contact_my_t[0] / m;
                this->bodies[b].nodes[i].contact_z_t[0] = 0;//this->bodies[b].nodes[i].contact_mz_t[0] / m;
            } else {
                this->bodies[b].nodes[i].contact_mx_t[0] = 0;
                this->bodies[b].nodes[i].contact_my_t[0] = 0;
                this->bodies[b].nodes[i].contact_mz_t[0] = 0;

                this->bodies[b].nodes[i].contact_x_t[0] = 0;
                this->bodies[b].nodes[i].contact_y_t[0] = 0;
                this->bodies[b].nodes[i].contact_z_t[0] = 0;
            }

        }
    }
    return;
}

void job_t::moveParticlesExplicit(){
    //if (this->use_3d==1) {
    for (size_t b = 0; b < this->num_bodies; b++) {
        //use Eigen Map to point to particle array
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;

        //use Eigen Map to point to node array
        size_t numRowsN = this->bodies[b].n;
        size_t numColsN = 1;

        //use to create dummy pvec
        Eigen::VectorXd pvec(numRowsP);

        for (size_t i = 0; i < this->bodies[b].n; i++) {
            double m = this->bodies[b].nodes[i].m[0];
            if (m != 0) {
                this->bodies[b].nodes[i].ux[0] = this->dt * this->bodies[b].nodes[i].contact_mx_t[0] / m;
                this->bodies[b].nodes[i].uy[0] = this->dt * this->bodies[b].nodes[i].contact_my_t[0] / m;
                this->bodies[b].nodes[i].uz[0] = this->dt * this->bodies[b].nodes[i].contact_mz_t[0] / m;

                this->bodies[b].nodes[i].diff_x_t[0] = this->dt * this->bodies[b].nodes[i].contact_fx[0] / m;
                this->bodies[b].nodes[i].diff_y_t[0] = this->dt * this->bodies[b].nodes[i].contact_fy[0] / m;
                this->bodies[b].nodes[i].diff_z_t[0] = this->dt * this->bodies[b].nodes[i].contact_fz[0] / m;
            } else {
                this->bodies[b].nodes[i].ux[0] = 0;
                this->bodies[b].nodes[i].uy[0] = 0;
                this->bodies[b].nodes[i].uz[0] = 0;

                this->bodies[b].nodes[i].diff_x_t[0] = 0;
                this->bodies[b].nodes[i].diff_y_t[0] = 0;
                this->bodies[b].nodes[i].diff_z_t[0] = 0;
            }
        }

        //map back to particles using S transpose
        this->bodies[b].particle_x += this->bodies[b].Phi.transpose() * this->bodies[b].node_ux;
        this->bodies[b].particle_y += this->bodies[b].Phi.transpose() * this->bodies[b].node_uy;
        this->bodies[b].particle_z += this->bodies[b].Phi.transpose() * this->bodies[b].node_uz;

        this->bodies[b].particle_ux += this->bodies[b].Phi.transpose() * this->bodies[b].node_ux;
        this->bodies[b].particle_uy += this->bodies[b].Phi.transpose() * this->bodies[b].node_uy;
        this->bodies[b].particle_uz += this->bodies[b].Phi.transpose() * this->bodies[b].node_uz;

        this->bodies[b].particle_x_t += this->bodies[b].Phi.transpose() * this->bodies[b].node_diff_x_t;
        this->bodies[b].particle_y_t += this->bodies[b].Phi.transpose() * this->bodies[b].node_diff_y_t;
        this->bodies[b].particle_z_t += this->bodies[b].Phi.transpose() * this->bodies[b].node_diff_z_t;

    }
    return;
    //} else {
    //    return this->moveParticlesExplicit2D();
    //}
}

void job_t::moveParticlesImplicit(){
    //if (this->use_3d==1) {
    //REQUIRES NODE mx_t TO BE STORED AT START OF TIMESTEP
    for (size_t b = 0; b < this->num_bodies; b++) {
        //use Eigen Map to point to particle array
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;

        //use Eigen Map to point to node array
        size_t numRowsN = this->bodies[b].n;
        size_t numColsN = 1;

        //use to create dummy pvec
        Eigen::VectorXd pvec(numRowsP);

        for (size_t i = 0; i < this->bodies[b].n; i++) {
            double m = this->bodies[b].nodes[i].m[0];
            if (m != 0) {
                this->bodies[b].nodes[i].ux[0] = this->dt * (this->bodies[b].node_mx_t_k[i] + this->bodies[b].nodes[i].contact_mx_t[0]) / (2*m);
                this->bodies[b].nodes[i].uy[0] = this->dt * (this->bodies[b].node_my_t_k[i] + this->bodies[b].nodes[i].contact_my_t[0]) / (2*m);
                this->bodies[b].nodes[i].uz[0] = this->dt * (this->bodies[b].node_mz_t_k[i] + this->bodies[b].nodes[i].contact_mz_t[0]) / (2*m);

                //this->bodies[b].nodes[i].uy[0] = this->dt * this->bodies[b].nodes[i].contact_my_t[0] / m;
                //this->bodies[b].nodes[i].uz[0] = this->dt * this->bodies[b].nodes[i].contact_mz_t[0] / m;

                this->bodies[b].nodes[i].diff_x_t[0] = this->dt * this->bodies[b].nodes[i].contact_fx[0] / m;
                this->bodies[b].nodes[i].diff_y_t[0] = this->dt * this->bodies[b].nodes[i].contact_fy[0] / m;
                this->bodies[b].nodes[i].diff_z_t[0] = this->dt * this->bodies[b].nodes[i].contact_fz[0] / m;
            } else {
                this->bodies[b].nodes[i].ux[0] = 0;
                this->bodies[b].nodes[i].uy[0] = 0;
                this->bodies[b].nodes[i].uz[0] = 0;

                this->bodies[b].nodes[i].diff_x_t[0] = 0;
                this->bodies[b].nodes[i].diff_y_t[0] = 0;
                this->bodies[b].nodes[i].diff_z_t[0] = 0;
            }
        }

        //map back to particles using S transpose
        this->bodies[b].particle_x += this->bodies[b].Phi.transpose() * this->bodies[b].node_ux;
        this->bodies[b].particle_y += this->bodies[b].Phi.transpose() * this->bodies[b].node_uy;
        this->bodies[b].particle_z += this->bodies[b].Phi.transpose() * this->bodies[b].node_uz;

        this->bodies[b].particle_ux += this->bodies[b].Phi.transpose() * this->bodies[b].node_ux;
        this->bodies[b].particle_uy += this->bodies[b].Phi.transpose() * this->bodies[b].node_uy;
        this->bodies[b].particle_uz += this->bodies[b].Phi.transpose() * this->bodies[b].node_uz;

        this->bodies[b].particle_x_t += this->bodies[b].Phi.transpose() * this->bodies[b].node_diff_x_t;
        this->bodies[b].particle_y_t += this->bodies[b].Phi.transpose() * this->bodies[b].node_diff_y_t;
        this->bodies[b].particle_z_t += this->bodies[b].Phi.transpose() * this->bodies[b].node_diff_z_t;

    }
    return;
    //} else {
    //    return this->moveParticlesExplicit2D();
    //}
}

void job_t::moveParticlesExplicit2D(){
    for (size_t b=0; b<this->num_bodies; b++) {
        //use Eigen Map to point to particle array
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;

        //use Eigen Map to point to node array
        size_t numRowsN = this->bodies[b].n;
        size_t numColsN = 1;

        //use to create dummy pvec
        //Eigen::VectorXd pvec(numRowsP);
        Eigen::VectorXd pvec(numRowsP);

        for (size_t i=0;i<this->bodies[b].n;i++){
            double m = this->bodies[b].nodes[i].m[0];
            if (m!=0){
                this->bodies[b].nodes[i].ux[0] = this->dt * this->bodies[b].nodes[i].contact_mx_t[0] / m;
                this->bodies[b].nodes[i].uy[0] = this->dt * this->bodies[b].nodes[i].contact_my_t[0] / m;
                //this->bodies[b].nodes[i].uz[0] = this->dt * this->bodies[b].nodes[i].contact_mz_t[0] / m;

                this->bodies[b].nodes[i].diff_x_t[0] = this->dt * this->bodies[b].nodes[i].contact_fx[0] / m;
                this->bodies[b].nodes[i].diff_y_t[0] = this->dt * this->bodies[b].nodes[i].contact_fy[0] / m;
                //this->bodies[b].nodes[i].diff_z_t[0] = this->dt * this->bodies[b].nodes[i].contact_fz[0] / m;
            } else {
                this->bodies[b].nodes[i].ux[0] = 0;
                this->bodies[b].nodes[i].uy[0] = 0;
                this->bodies[b].nodes[i].uz[0] = 0;

                this->bodies[b].nodes[i].diff_x_t[0] = 0;
                this->bodies[b].nodes[i].diff_y_t[0] = 0;
                this->bodies[b].nodes[i].diff_z_t[0] = 0;
            }
        }

        //map back to particles using S transpose
        this->bodies[b].particle_x += this->bodies[b].Phi.transpose()*this->bodies[b].node_ux;
        this->bodies[b].particle_y += this->bodies[b].Phi.transpose()*this->bodies[b].node_uy;
        //this->bodies[b].particle_z += this->bodies[b].Phi.transpose()*this->bodies[b].node_uz;

        this->bodies[b].particle_ux += this->bodies[b].Phi.transpose()*this->bodies[b].node_ux;
        this->bodies[b].particle_uy += this->bodies[b].Phi.transpose()*this->bodies[b].node_uy;
        //this->bodies[b].particle_uz += this->bodies[b].Phi.transpose()*this->bodies[b].node_uz;

        this->bodies[b].particle_x_t += this->bodies[b].Phi.transpose()*this->bodies[b].node_diff_x_t;
        this->bodies[b].particle_y_t += this->bodies[b].Phi.transpose()*this->bodies[b].node_diff_y_t;
        //this->bodies[b].particle_z_t += this->bodies[b].Phi.transpose()*this->bodies[b].node_diff_z_t;

    }
    return;
}

void job_t::calculateStrainRate() {
    //if (this->use_3d==1) {
    for (size_t b = 0; b < this->num_bodies; b++) {
        //map nodal velocities with Eigen
        size_t numRowsN = this->bodies[b].n;
        size_t numColsN = 1;

        //nodal velocities
        this->bodies[b].node_contact_x_t =
                this->bodies[b].node_contact_mx_t.array() / this->bodies[b].node_m.array();
        this->bodies[b].node_contact_y_t =
                this->bodies[b].node_contact_my_t.array() / this->bodies[b].node_m.array();
        this->bodies[b].node_contact_z_t =
                this->bodies[b].node_contact_mz_t.array() / this->bodies[b].node_m.array();

        //use to create dummy pvec
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;
        Eigen::VectorXd pvec(numRowsP);

        //calculate particle[i].L[9]
        pvec = this->bodies[b].gradPhiX.transpose() * this->bodies[b].node_contact_x_t;
        for (size_t i = 0; i < this->bodies[b].p; i++) {
            this->bodies[b].particles[i].L[XX] = pvec[i];
        }

        pvec = this->bodies[b].gradPhiY.transpose() * this->bodies[b].node_contact_x_t;
        for (size_t i = 0; i < this->bodies[b].p; i++) {
            this->bodies[b].particles[i].L[XY] = pvec[i];
        }

        pvec = this->bodies[b].gradPhiZ.transpose() * this->bodies[b].node_contact_x_t;
        for (size_t i = 0; i < this->bodies[b].p; i++) {
            this->bodies[b].particles[i].L[XZ] = pvec[i];
        }

        pvec = this->bodies[b].gradPhiX.transpose() * this->bodies[b].node_contact_y_t;
        for (size_t i = 0; i < this->bodies[b].p; i++) {
            this->bodies[b].particles[i].L[YX] = pvec[i];
        }

        pvec = this->bodies[b].gradPhiY.transpose() * this->bodies[b].node_contact_y_t;
        for (size_t i = 0; i < this->bodies[b].p; i++) {
            this->bodies[b].particles[i].L[YY] = pvec[i];
        }

        pvec = this->bodies[b].gradPhiZ.transpose() * this->bodies[b].node_contact_y_t;
        for (size_t i = 0; i < this->bodies[b].p; i++) {
            this->bodies[b].particles[i].L[YZ] = pvec[i];
        }

        pvec = this->bodies[b].gradPhiX.transpose() * this->bodies[b].node_contact_z_t;
        for (size_t i = 0; i < this->bodies[b].p; i++) {
            this->bodies[b].particles[i].L[ZX] = pvec[i];
        }

        pvec = this->bodies[b].gradPhiY.transpose() * this->bodies[b].node_contact_z_t;
        for (size_t i = 0; i < this->bodies[b].p; i++) {
            this->bodies[b].particles[i].L[ZY] = pvec[i];
        }

        pvec = this->bodies[b].gradPhiZ.transpose() * this->bodies[b].node_contact_z_t;
        for (size_t i = 0; i < this->bodies[b].p; i++) {
            this->bodies[b].particles[i].L[ZZ] = pvec[i];
        }
    }
    return;
    //} else {
    //    return this->calculateStrainRate2D();
    //}
}

void job_t::calculateStrainRate2D() {
    for (size_t b=0; b<this->num_bodies; b++) {
        //map nodal velocities with Eigen
        size_t numRowsN = this->bodies[b].n;
        size_t numColsN = 1;

        //nodal velocities
        this->bodies[b].node_contact_x_t = this->bodies[b].node_contact_mx_t.array()/this->bodies[b].node_m.array();
        this->bodies[b].node_contact_y_t = this->bodies[b].node_contact_my_t.array()/this->bodies[b].node_m.array();
        //this->bodies[b].node_contact_z_t = this->bodies[b].node_contact_mz_t.array()/this->bodies[b].node_m.array();

        //use to create dummy pvec
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;
        //Eigen::MatrixXd pvec(numRowsP,numColsP);
        Eigen::VectorXd pvec(numRowsP);

        //calculate particle[i].L[9]
        pvec = this->bodies[b].gradPhiX.transpose()*this->bodies[b].node_contact_x_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[XX] = pvec[i];
        }

        pvec = this->bodies[b].gradPhiY.transpose()*this->bodies[b].node_contact_x_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[XY] = pvec[i];
        }

        //pvec = this->bodies[b].gradPhiZ.transpose()*this->bodies[b].node_contact_x_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[XZ] = 0;
        }

        pvec = this->bodies[b].gradPhiX.transpose()*this->bodies[b].node_contact_y_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[YX] = pvec[i];
        }

        pvec = this->bodies[b].gradPhiY.transpose()*this->bodies[b].node_contact_y_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[YY] = pvec[i];
        }

        //pvec = this->bodies[b].gradPhiZ.transpose()*this->bodies[b].node_contact_y_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[YZ] = 0;
        }

        //pvec = this->bodies[b].gradPhiX.transpose()*this->bodies[b].node_contact_z_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[ZX] = 0;
        }

        //pvec = this->bodies[b].gradPhiY.transpose()*this->bodies[b].node_contact_z_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[ZY] = 0;
        }

        //pvec = this->bodies[b].gradPhiZ.transpose()*this->bodies[b].node_contact_z_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[ZZ] = 0;
        }
    }
    return;
}

void job_t::updateDensity(){
    //update density of particles per sachiths code
    for (size_t b=0;b<this->num_bodies;b++){
        for (size_t i=0;i<this->bodies[b].p;i++){
            double trL = 0;
            tensor_trace3(&trL,this->bodies[b].particles[i].L);
            this->bodies[b].particles[i].v[0] *= exp(this->dt * trL);
        }
    }
    return;
}

void job_t::updateTrialDensity(){
    //update density of particles per sachiths code
    for (size_t b=0;b<this->num_bodies;b++){
        for (size_t i=0;i<this->bodies[b].p;i++){
            double trL = 0;
            tensor_trace3(&trL,this->bodies[b].particles[i].L);
            this->bodies[b].particles[i].v_trial[0] = this->bodies[b].particles[i].v[0] * exp(this->dt * trL);
        }
    }
    return;
}

void job_t::updateStress(){
    //calculate stress
    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].material.calculate_stress(&(this->bodies[b]),this->dt);
    }
    return;
}

void job_t::updateTrialStress(){
    //calculate stress
    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].material.calculate_stress_implicit(&(this->bodies[b]),this->dt);
    }
    return;
}

void job_t::mapTrialStress2Grid() {
    for (size_t b=0; b<this->num_bodies; b++) {
        //use Eigen Map to point to particle array
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;

        Eigen::MatrixXd p_m_bx(numRowsP, numColsP);
        Eigen::MatrixXd p_m_by(numRowsP, numColsP);
        Eigen::MatrixXd p_m_bz(numRowsP, numColsP);
        p_m_bx = this->bodies[b].particle_m.array() * this->bodies[b].particle_bx.array();
        p_m_by = this->bodies[b].particle_m.array() * this->bodies[b].particle_by.array();
        p_m_bz = this->bodies[b].particle_m.array() * this->bodies[b].particle_bz.array();

        this->bodies[b].node_fx = this->bodies[b].Phi * p_m_bx; //need to add stress
        this->bodies[b].node_fy = this->bodies[b].Phi * p_m_by; //need to add stress
        this->bodies[b].node_fz = this->bodies[b].Phi * p_m_bz; //need to add stress

        //use to create dummy pvec and ones
        Eigen::VectorXd pvec(numRowsP);

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].Ttrial[XX];
        }
        this->bodies[b].node_fx -= this->bodies[b].gradPhiX*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].Ttrial[XY];
        }
        this->bodies[b].node_fx -= this->bodies[b].gradPhiY*pvec;
        this->bodies[b].node_fy -= this->bodies[b].gradPhiX*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].Ttrial[XZ];
        }
        this->bodies[b].node_fx -= this->bodies[b].gradPhiZ*pvec;
        this->bodies[b].node_fz -= this->bodies[b].gradPhiX*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].Ttrial[YY];
        }
        this->bodies[b].node_fy -= this->bodies[b].gradPhiY*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].Ttrial[YZ];
        }
        this->bodies[b].node_fy -= this->bodies[b].gradPhiZ*pvec;
        this->bodies[b].node_fz -= this->bodies[b].gradPhiY*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec[i] = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].Ttrial[ZZ];
        }
        this->bodies[b].node_fz -= this->bodies[b].gradPhiZ*pvec;
    }
    return;
}

void job_t::calculateImplicitResidual() {
    for (size_t b=0;b<this->num_bodies;b++) {
        /*this->bodies[b].Rx = this->bodies[b].node_x_t_trial.array()*this->bodies[b].node_m.array()
                                  - this->dt*this->bodies[b].node_fx_L.array() + this->dt*this->bodies[b].node_fx_k.array()
                                  - this->bodies[b].node_x_t_explicit.array()*this->bodies[b].node_m.array();

        this->bodies[b].Ry = this->bodies[b].node_y_t_trial.array()*this->bodies[b].node_m.array()
                                  - this->dt*this->bodies[b].node_fy_L.array() + this->dt*this->bodies[b].node_fy_k.array()
                                  - this->bodies[b].node_y_t_explicit.array()*this->bodies[b].node_m.array();

        this->bodies[b].Rz = this->bodies[b].node_z_t_trial.array()*this->bodies[b].node_m.array()
                                  - this->dt*this->bodies[b].node_fz_L.array() + this->dt*this->bodies[b].node_fz_k.array()
                                  - this->bodies[b].node_z_t_explicit.array()*this->bodies[b].node_m.array();
        */

        this->bodies[b].Rx = this->bodies[b].node_x_t_trial.array()*this->bodies[b].node_m.array();
        this->bodies[b].Rx -= this->dt*this->bodies[b].node_fx_L;
        this->bodies[b].Rx -= this->bodies[b].node_mx_t_k;

        this->bodies[b].Ry = this->bodies[b].node_y_t_trial.array()*this->bodies[b].node_m.array();
        this->bodies[b].Ry -= this->dt*this->bodies[b].node_fy_L;
        this->bodies[b].Ry -= this->bodies[b].node_my_t_k;

        this->bodies[b].Rz = this->bodies[b].node_z_t_trial.array()*this->bodies[b].node_m.array();
        this->bodies[b].Rz -= this->dt*this->bodies[b].node_fz_L;
        this->bodies[b].Rz -= this->bodies[b].node_mz_t_k;

        for (size_t i=0;i<this->num_nodes;i++){
            if (this->u_dirichlet_mask[i] != 0 || this->bodies[b].node_m[i] <= TOL){
                //nodal boundary
                this->bodies[b].Rx[i] = 0;
                this->bodies[b].Ry[i] = 0;
                this->bodies[b].Rz[i] = 0;
            } else {
                //this->bodies[b].Rx[i] /= this->bodies[b].node_m[i];
                //this->bodies[b].Ry[i] /= this->bodies[b].node_m[i];
                //this->bodies[b].Rz[i] /= this->bodies[b].node_m[i];
            }

        }
    }
}

void job_t::moveGridImplicitCG() {
    //full conjugate gradient mathod for solving Js = -F(v)

    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].node_mx_t_k = this->bodies[b].node_contact_mx_t;
        this->bodies[b].node_my_t_k = this->bodies[b].node_contact_my_t;
        this->bodies[b].node_mz_t_k = this->bodies[b].node_contact_mz_t;
    }

    //move grid
    //this->moveGridExplicit();

    //save forces and initial trial velocity
    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].node_fx_k = this->bodies[b].node_contact_fx;
        this->bodies[b].node_fy_k = this->bodies[b].node_contact_fy;
        this->bodies[b].node_fz_k = this->bodies[b].node_contact_fz;

        /*this->bodies[b].node_x_t_explicit.setZero();
        this->bodies[b].node_y_t_explicit.setZero();
        this->bodies[b].node_z_t_explicit.setZero();

        for (size_t i=0;i<this->num_nodes;i++){
            if (this->bodies[b].node_m[i] > TOL) {
                this->bodies[b].node_x_t_explicit[i] = this->bodies[b].node_contact_mx_t[i] / this->bodies[b].node_m[i];
                this->bodies[b].node_y_t_explicit[i] = this->bodies[b].node_contact_my_t[i] / this->bodies[b].node_m[i];
                this->bodies[b].node_z_t_explicit[i] = this->bodies[b].node_contact_mz_t[i] / this->bodies[b].node_m[i];
            }
        }

        this->bodies[b].node_x_t_n = this->bodies[b].node_x_t_explicit;
        this->bodies[b].node_y_t_n = this->bodies[b].node_y_t_explicit;
        this->bodies[b].node_z_t_n = this->bodies[b].node_z_t_explicit;
         */

        this->bodies[b].node_x_t_n.setZero();
        this->bodies[b].node_y_t_n.setZero();
        this->bodies[b].node_z_t_n.setZero();

        for (size_t i=0;i<this->num_nodes;i++){
            if (this->bodies[b].node_m[i] > TOL) {
                this->bodies[b].node_x_t_n[i] = this->bodies[b].node_contact_mx_t[i] / this->bodies[b].node_m[i];
                this->bodies[b].node_y_t_n[i] = this->bodies[b].node_contact_my_t[i] / this->bodies[b].node_m[i];
                this->bodies[b].node_z_t_n[i] = this->bodies[b].node_contact_mz_t[i] / this->bodies[b].node_m[i];
            }
        }

        this->bodies[b].node_x_t_trial = this->bodies[b].node_x_t_n;
        this->bodies[b].node_y_t_trial = this->bodies[b].node_y_t_n;
        this->bodies[b].node_z_t_trial = this->bodies[b].node_z_t_n;
    }

    //calculate L on particles
    this->calculateStrainRate();

    //update particle densities
    this->updateTrialDensity();

    //material stress update
    this->updateTrialStress();

    //add body forces
    time_varying_loads(this);

    //map particle stress back to nodes
    this->mapTrialStress2Grid();

    //add contact forces
    this->addContactForces();

    //enforce boundary conditions
    this->addBoundaryConditions();

    //save explicit solution
    for (size_t b=0;b<this->num_bodies;b++){
        //this->bodies[b].node_x_t_explicit = this->bodies[b].node_contact_mx_t.array()/this->bodies[b].node_m.array();
        //this->bodies[b].node_y_t_explicit = this->bodies[b].node_contact_my_t.array()/this->bodies[b].node_m.array();
        //this->bodies[b].node_z_t_explicit = this->bodies[b].node_contact_mz_t.array()/this->bodies[b].node_m.array();

        this->bodies[b].node_fx_L = this->bodies[b].node_contact_fx;
        this->bodies[b].node_fy_L = this->bodies[b].node_contact_fy;
        this->bodies[b].node_fz_L = this->bodies[b].node_contact_fz;
    }

    //calculate residual for expicit step
    this->calculateImplicitResidual();

    double rhoSum = 0;
    double rhoTOL = 0;
    rhoTOL = this->newtonTOL;//R_TOL;//*this->num_bodies*this->num_nodes;

    size_t nIter = 0;
    do {
        //calculate norm from residual on both bodies and save residual of trial velocity
        //setup iteration for s
        rhoSum = 0;
        for (size_t b = 0; b < this->num_bodies; b++) {
            this->bodies[b].Rvx = this->bodies[b].Rx;
            this->bodies[b].Rvy = this->bodies[b].Ry;
            this->bodies[b].Rvz = this->bodies[b].Rz;

            this->bodies[b].rk.setZero();

            this->bodies[b].rk << -this->bodies[b].Rvx, -this->bodies[b].Rvy, -this->bodies[b].Rvz;

            this->bodies[b].rhok = this->bodies[b].rk.squaredNorm();

            rhoSum += this->bodies[b].rhok;

            this->bodies[b].pk = this->bodies[b].rk / std::sqrt(this->bodies[b].rhok);

            this->bodies[b].sk.setZero();

            /*if (this->bodies[b].rhok > 2){
                this->bodies[b].rhok = 1;
            }*/
        }

        std::cout << "RHO: " << rhoSum << " ?< " << rhoTOL <<  " vnorm: " << this->bodies[0].node_x_t_n.lpNorm<Eigen::Infinity>() <<
        "," << this->bodies[0].node_y_t_n.lpNorm<Eigen::Infinity>() <<
        "," << this->bodies[0].node_z_t_n.lpNorm<Eigen::Infinity>() << std::endl;

        //solve for 's' to iterate 'v' [Sulsky 2003]
        size_t k = 0;
        while (/*k < 3*this->num_nodes &&*/ rhoSum > rhoTOL && rhoSum < R_MAX) { //(this->num_bodies * this->num_nodes * R_TOL)) {
            double h = this->linearStepSize;//-6; //per paper suggestion

            for (size_t b = 0; b < this->num_bodies; b++) {
                double vNorm = std::sqrt(
                        this->bodies[b].node_x_t_n.squaredNorm() + this->bodies[b].node_y_t_n.squaredNorm()
                        + this->bodies[b].node_z_t_n.squaredNorm());
                if (vNorm <= TOL) {
                    vNorm = 1.0;
                }
                double sNorm = this->bodies[b].pk.norm();
                if (sNorm>TOL) {
                    this->bodies[b].node_x_t_trial =
                            this->bodies[b].node_x_t_n + h * vNorm * this->bodies[b].pk.segment(0,this->bodies[b].n) / sNorm;
                    this->bodies[b].node_y_t_trial =
                            this->bodies[b].node_y_t_n + h * vNorm * this->bodies[b].pk.segment(this->bodies[b].n,this->bodies[b].n)/ sNorm;
                    this->bodies[b].node_z_t_trial =
                            this->bodies[b].node_z_t_n + h * vNorm * this->bodies[b].pk.segment(2*this->bodies[b].n,this->bodies[b].n) / sNorm;
                } else {
                    this->bodies[b].node_x_t_trial = this->bodies[b].node_x_t_n;
                    this->bodies[b].node_y_t_trial = this->bodies[b].node_y_t_n;
                    this->bodies[b].node_z_t_trial = this->bodies[b].node_z_t_n;
                }

                this->bodies[b].node_mx_t = this->bodies[b].node_x_t_trial.array() * this->bodies[b].node_m.array();
                this->bodies[b].node_my_t = this->bodies[b].node_y_t_trial.array() * this->bodies[b].node_m.array();
                this->bodies[b].node_mz_t = this->bodies[b].node_z_t_trial.array() * this->bodies[b].node_m.array();

                this->bodies[b].node_x_t = this->bodies[b].node_x_t_trial;
                this->bodies[b].node_y_t = this->bodies[b].node_y_t_trial;
                this->bodies[b].node_z_t = this->bodies[b].node_z_t_trial;
            }

            //add contact forces
            this->addContactForces();

            //enforce boundary conditions
            this->addBoundaryConditions();

            //calculate L on particles
            this->calculateStrainRate();

            //update particle densities
            this->updateTrialDensity();

            //material stress update
            this->updateTrialStress();

            this->mapTrialStress2Grid();

            //add contact forces
            this->addContactForces();

            //enforce boundary conditions
            this->addBoundaryConditions();

            for (size_t b=0;b<this->num_bodies;b++) {
                this->bodies[b].node_fx_L = this->bodies[b].node_contact_fx;
                this->bodies[b].node_fy_L = this->bodies[b].node_contact_fy;
                this->bodies[b].node_fz_L = this->bodies[b].node_contact_fz;
            }

            this->calculateImplicitResidual();

            for (size_t b = 0; b < this->num_bodies; b++) {
                double vNorm = std::sqrt(
                        this->bodies[b].node_x_t_n.squaredNorm() + this->bodies[b].node_y_t_n.squaredNorm() +
                        this->bodies[b].node_z_t_n.squaredNorm());
                if (vNorm <= TOL) {
                    vNorm = 1.0;
                }
                double sNorm = this->bodies[b].pk.norm();
                this->bodies[b].DhRx = sNorm / (h * vNorm) * (this->bodies[b].Rx - this->bodies[b].Rvx);
                this->bodies[b].DhRy = sNorm / (h * vNorm) * (this->bodies[b].Ry - this->bodies[b].Rvy);
                this->bodies[b].DhRz = sNorm / (h * vNorm) * (this->bodies[b].Rz - this->bodies[b].Rvz);

                this->bodies[b].wk << this->bodies[b].DhRx, this->bodies[b].DhRy, this->bodies[b].DhRz;
                this->bodies[b].ak = this->bodies[b].rhok / (this->bodies[b].pk.transpose() * this->bodies[b].wk);

                this->bodies[b].sk += this->bodies[b].ak * this->bodies[b].pk;

                this->bodies[b].rk -= this->bodies[b].ak * this->bodies[b].wk;

                this->bodies[b].bk = this->bodies[b].rk.squaredNorm() / this->bodies[b].rhok;

                this->bodies[b].rhok = this->bodies[b].rk.squaredNorm();

                this->bodies[b].pk = this->bodies[b].rk + this->bodies[b].bk * this->bodies[b].pk;
            }
            k += 1;
            rhoSum = 0;
            for (size_t b = 0; b < this->num_bodies; b++) {
                rhoSum += this->bodies[b].rhok;
            }
            std::cout << "\rn: " << nIter << " k: " << k <<  " r: " << rhoSum << "      \r" << std::flush;
        }

        for (size_t b = 0; b < this->num_bodies; b++) {
            this->bodies[b].node_x_t_n += this->bodies[b].sk.segment(0,this->bodies[b].n);
            this->bodies[b].node_y_t_n += this->bodies[b].sk.segment(this->bodies[b].n,this->bodies[b].n);
            this->bodies[b].node_z_t_n += this->bodies[b].sk.segment(2*this->bodies[b].n,this->bodies[b].n);

            this->bodies[b].node_x_t_trial = this->bodies[b].node_x_t_n;
            this->bodies[b].node_y_t_trial = this->bodies[b].node_y_t_n;
            this->bodies[b].node_z_t_trial = this->bodies[b].node_z_t_n;

            this->bodies[b].node_mx_t = this->bodies[b].node_x_t_trial.array() * this->bodies[b].node_m.array();
            this->bodies[b].node_my_t = this->bodies[b].node_y_t_trial.array() * this->bodies[b].node_m.array();
            this->bodies[b].node_mz_t = this->bodies[b].node_z_t_trial.array() * this->bodies[b].node_m.array();

            this->bodies[b].node_x_t = this->bodies[b].node_x_t_trial;
            this->bodies[b].node_y_t = this->bodies[b].node_y_t_trial;
            this->bodies[b].node_z_t = this->bodies[b].node_z_t_trial;
        }

        std::cout << "n: " << nIter << " k: " << k <<  " r: " << rhoSum << std::endl;
        nIter += 1;

        //add contact forces
        this->addContactForces();

        //enforce boundary conditions
        this->addBoundaryConditions();
        //calculate L on particles
        this->calculateStrainRate();

        //update particle densities
        this->updateTrialDensity();

        //material stress update
        this->updateTrialStress();

        this->mapTrialStress2Grid();
        //calculate residual for expicit step

        //add contact forces
        this->addContactForces();

        //enforce boundary conditions
        this->addBoundaryConditions();

        for (size_t b=0;b<this->num_bodies;b++) {
            this->bodies[b].node_fx_L = this->bodies[b].node_contact_fx;
            this->bodies[b].node_fy_L = this->bodies[b].node_contact_fy;
            this->bodies[b].node_fz_L = this->bodies[b].node_contact_fz;
        }

        this->calculateImplicitResidual();

        //setup for next iteration
        rhoSum = 0;
        for (size_t b = 0; b < this->num_bodies; b++) {
            //this->bodies[b].Rvx = this->bodies[b].Rx;
            //this->bodies[b].Rvy = this->bodies[b].Ry;
            //this->bodies[b].Rvz = this->bodies[b].Rz;

            this->bodies[b].rk << -this->bodies[b].Rx, -this->bodies[b].Ry, -this->bodies[b].Rz;

            this->bodies[b].rhok = this->bodies[b].rk.squaredNorm();

            rhoSum += this->bodies[b].rhok;

            //this->bodies[b].pk = this->bodies[b].rk / std::sqrt(this->bodies[b].rhok);

            //this->bodies[b].sk.setZero();
        }
        //std::cout << "n: " << nIter << " k: " << k <<  " r: " << rhoSum << std::endl;
    } while (rhoSum>rhoTOL && rhoSum<R_MAX); //(this->num_bodies * this->num_nodes * R_TOL));

    //catch error where rhoSum is nan
    if (!std::isfinite(rhoSum) || rhoSum>=R_MAX){
        std::cout << "Error: Residual in mpmStepUSLImplicit() is infinite." << std::endl;
        /*if (this->dt*0.5 < this->dt_minimum){
            std::cout << "dt below threshhold. Exiting." << std::endl;
            exit(0);
        } else {
            //reset step
            this->t -= dt;
            this->stepcount -= 1;

            //adjust timestep
            this->dt *= 0.5;
            std::cout << "Reducing dt: " << this->dt << std::endl;

            //restart step
            this->mpmStepUSLImplicit();

            //reset timestep
            this->dt = this->dt_base;

            return 1;
        }*/
        exit(0);
    }

    return;
}

void job_t::moveGridImplicitBiCGSTAB() {
    //full biconjugate gradient stabilized mathod for solving Js = -F(v)

    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].node_mx_t_k = this->bodies[b].node_contact_mx_t;
        this->bodies[b].node_my_t_k = this->bodies[b].node_contact_my_t;
        this->bodies[b].node_mz_t_k = this->bodies[b].node_contact_mz_t;
    }

    //move grid
    //this->moveGridExplicit();

    //save forces and initial trial velocity
    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].node_fx_k = this->bodies[b].node_contact_fx;
        this->bodies[b].node_fy_k = this->bodies[b].node_contact_fy;
        this->bodies[b].node_fz_k = this->bodies[b].node_contact_fz;

        this->bodies[b].node_x_t_n.setZero();
        this->bodies[b].node_y_t_n.setZero();
        this->bodies[b].node_z_t_n.setZero();

        for (size_t i=0;i<this->num_nodes;i++){
            if (this->bodies[b].node_m[i] > TOL) {
                this->bodies[b].node_x_t_n[i] = this->bodies[b].node_contact_mx_t[i] / this->bodies[b].node_m[i];
                this->bodies[b].node_y_t_n[i] = this->bodies[b].node_contact_my_t[i] / this->bodies[b].node_m[i];
                this->bodies[b].node_z_t_n[i] = this->bodies[b].node_contact_mz_t[i] / this->bodies[b].node_m[i];
            }
        }

        this->bodies[b].node_x_t_trial = this->bodies[b].node_x_t_n;
        this->bodies[b].node_y_t_trial = this->bodies[b].node_y_t_n;
        this->bodies[b].node_z_t_trial = this->bodies[b].node_z_t_n;
    }

    //calculate L on particles
    this->calculateStrainRate();

    //update particle densities
    this->updateTrialDensity();

    //material stress update
    this->updateTrialStress();

    //add body forces
    time_varying_loads(this);

    //map particle stress back to nodes
    this->mapTrialStress2Grid();

    //add contact forces
    this->addContactForces();

    //enforce boundary conditions
    this->addBoundaryConditions();

    //save final forces on nodes
    for (size_t b=0;b<this->num_bodies;b++){

        this->bodies[b].node_fx_L = this->bodies[b].node_contact_fx;
        this->bodies[b].node_fy_L = this->bodies[b].node_contact_fy;
        this->bodies[b].node_fz_L = this->bodies[b].node_contact_fz;
    }

    //calculate initial residual
    this->calculateImplicitResidual();

    double rhoSum = 0;
    double rhoTOL = 0;
    rhoTOL = this->newtonTOL;//R_TOL;//*this->num_bodies*this->num_nodes;

    size_t nIter = 0;
    do {
        //calculate norm from residual on both bodies and save residual of trial velocity
        //setup iteration for s
        rhoSum = 0;
        for (size_t b = 0; b < this->num_bodies; b++) {
            this->bodies[b].Rvx = this->bodies[b].Rx;
            this->bodies[b].Rvy = this->bodies[b].Ry;
            this->bodies[b].Rvz = this->bodies[b].Rz;

            //initialize bicgstab variables
            this->bodies[b].sk.setZero();
            this->bodies[b].wk.setZero();
            this->bodies[b].pk.setZero();
            this->bodies[b].ak = 1;
            this->bodies[b].rhok = 1;
            this->bodies[b].ok = 1;
            this->bodies[b].rk.setZero();

            this->bodies[b].rk << -this->bodies[b].Rvx, -this->bodies[b].Rvy, -this->bodies[b].Rvz;
            this->bodies[b].r0 = this->bodies[b].rk;

            rhoSum += this->bodies[b].rk.squaredNorm();
        }

        std::cout << "RHO: " << rhoSum << " ?< " << rhoTOL <<  " vnorm: " << this->bodies[0].node_x_t_n.lpNorm<Eigen::Infinity>() <<
        "," << this->bodies[0].node_y_t_n.lpNorm<Eigen::Infinity>() <<
        "," << this->bodies[0].node_z_t_n.lpNorm<Eigen::Infinity>() << std::endl;

        //solve for 's' to iterate 'v' [Sulsky 2003]
        size_t k = 0;
        while (/*k < 3*this->num_nodes &&*/ rhoSum > rhoTOL && rhoSum < R_MAX) {
            double h = this->linearStepSize;

            for (size_t b = 0; b < this->num_bodies; b++) {

                this->bodies[b].bk = (this->bodies[b].r0.dot(this->bodies[b].rk))/(this->bodies[b].rhok)*this->bodies[b].ak/this->bodies[b].ok;
                this->bodies[b].rhok = this->bodies[b].r0.dot(this->bodies[b].rk);
                this->bodies[b].pk = this->bodies[b].rk + this->bodies[b].bk*(this->bodies[b].pk - this->bodies[b].ok*this->bodies[b].wk);

                //wk = DhDF(vn,pk)
                double vNorm = std::sqrt(
                        this->bodies[b].node_x_t_n.squaredNorm() + this->bodies[b].node_y_t_n.squaredNorm()
                        + this->bodies[b].node_z_t_n.squaredNorm());
                if (vNorm <= TOL) {
                    vNorm = 1.0;
                }
                double sNorm = this->bodies[b].pk.norm();
                if (sNorm>TOL) {
                    this->bodies[b].node_x_t_trial =
                            this->bodies[b].node_x_t_n + h * vNorm * this->bodies[b].pk.segment(0,this->bodies[b].n) / sNorm;
                    this->bodies[b].node_y_t_trial =
                            this->bodies[b].node_y_t_n + h * vNorm * this->bodies[b].pk.segment(this->bodies[b].n,this->bodies[b].n)/ sNorm;
                    this->bodies[b].node_z_t_trial =
                            this->bodies[b].node_z_t_n + h * vNorm * this->bodies[b].pk.segment(2*this->bodies[b].n,this->bodies[b].n) / sNorm;
                } else {
                    this->bodies[b].node_x_t_trial = this->bodies[b].node_x_t_n;
                    this->bodies[b].node_y_t_trial = this->bodies[b].node_y_t_n;
                    this->bodies[b].node_z_t_trial = this->bodies[b].node_z_t_n;
                }

                this->bodies[b].node_mx_t = this->bodies[b].node_x_t_trial.array() * this->bodies[b].node_m.array();
                this->bodies[b].node_my_t = this->bodies[b].node_y_t_trial.array() * this->bodies[b].node_m.array();
                this->bodies[b].node_mz_t = this->bodies[b].node_z_t_trial.array() * this->bodies[b].node_m.array();

                this->bodies[b].node_x_t = this->bodies[b].node_x_t_trial;
                this->bodies[b].node_y_t = this->bodies[b].node_y_t_trial;
                this->bodies[b].node_z_t = this->bodies[b].node_z_t_trial;
            }

            //add contact forces
            this->addContactForces();

            //enforce boundary conditions
            this->addBoundaryConditions();

            //calculate L on particles
            this->calculateStrainRate();

            //update particle densities
            this->updateTrialDensity();

            //material stress update
            this->updateTrialStress();

            this->mapTrialStress2Grid();

            //add contact forces
            this->addContactForces();

            //enforce boundary conditions
            this->addBoundaryConditions();

            for (size_t b=0;b<this->num_bodies;b++) {
                this->bodies[b].node_fx_L = this->bodies[b].node_contact_fx;
                this->bodies[b].node_fy_L = this->bodies[b].node_contact_fy;
                this->bodies[b].node_fz_L = this->bodies[b].node_contact_fz;
            }

            this->calculateImplicitResidual();

            for (size_t b = 0; b < this->num_bodies; b++) {
                double vNorm = std::sqrt(
                        this->bodies[b].node_x_t_n.squaredNorm() + this->bodies[b].node_y_t_n.squaredNorm() +
                        this->bodies[b].node_z_t_n.squaredNorm());
                if (vNorm <= TOL) {
                    vNorm = 1.0;
                }
                double sNorm = this->bodies[b].pk.norm();
                this->bodies[b].DhRx = sNorm / (h * vNorm) * (this->bodies[b].Rx - this->bodies[b].Rvx);
                this->bodies[b].DhRy = sNorm / (h * vNorm) * (this->bodies[b].Ry - this->bodies[b].Rvy);
                this->bodies[b].DhRz = sNorm / (h * vNorm) * (this->bodies[b].Rz - this->bodies[b].Rvz);

                this->bodies[b].wk << this->bodies[b].DhRx, this->bodies[b].DhRy, this->bodies[b].DhRz;
                this->bodies[b].ak = this->bodies[b].rhok / (this->bodies[b].r0.dot(this->bodies[b].wk));
                this->bodies[b].hk = this->bodies[b].sk + this->bodies[b].ak * this->bodies[b].pk;

                if(!std::isfinite(this->bodies[b].ak)){
                    std::cout << "ak is infinite\n";
                }

                //check hk for convergence?

                this->bodies[b].qk = this->bodies[b].rk - this->bodies[b].ak * this->bodies[b].wk;

                //calculate t = DhDF(vn,q)
                sNorm = this->bodies[b].qk.norm();
                if (sNorm>TOL) {
                    this->bodies[b].node_x_t_trial =
                            this->bodies[b].node_x_t_n + h * vNorm * this->bodies[b].qk.segment(0,this->bodies[b].n) / sNorm;
                    this->bodies[b].node_y_t_trial =
                            this->bodies[b].node_y_t_n + h * vNorm * this->bodies[b].qk.segment(this->bodies[b].n,this->bodies[b].n)/ sNorm;
                    this->bodies[b].node_z_t_trial =
                            this->bodies[b].node_z_t_n + h * vNorm * this->bodies[b].qk.segment(2*this->bodies[b].n,this->bodies[b].n) / sNorm;
                } else {
                    this->bodies[b].node_x_t_trial = this->bodies[b].node_x_t_n;
                    this->bodies[b].node_y_t_trial = this->bodies[b].node_y_t_n;
                    this->bodies[b].node_z_t_trial = this->bodies[b].node_z_t_n;
                }

                this->bodies[b].node_mx_t = this->bodies[b].node_x_t_trial.array() * this->bodies[b].node_m.array();
                this->bodies[b].node_my_t = this->bodies[b].node_y_t_trial.array() * this->bodies[b].node_m.array();
                this->bodies[b].node_mz_t = this->bodies[b].node_z_t_trial.array() * this->bodies[b].node_m.array();

                this->bodies[b].node_x_t = this->bodies[b].node_x_t_trial;
                this->bodies[b].node_y_t = this->bodies[b].node_y_t_trial;
                this->bodies[b].node_z_t = this->bodies[b].node_z_t_trial;
            }

            //add contact forces
            this->addContactForces();

            //enforce boundary conditions
            this->addBoundaryConditions();

            //calculate L on particles
            this->calculateStrainRate();

            //update particle densities
            this->updateTrialDensity();

            //material stress update
            this->updateTrialStress();

            this->mapTrialStress2Grid();

            //add contact forces
            this->addContactForces();

            //enforce boundary conditions
            this->addBoundaryConditions();

            for (size_t b=0;b<this->num_bodies;b++) {
                this->bodies[b].node_fx_L = this->bodies[b].node_contact_fx;
                this->bodies[b].node_fy_L = this->bodies[b].node_contact_fy;
                this->bodies[b].node_fz_L = this->bodies[b].node_contact_fz;
            }

            this->calculateImplicitResidual();

            for (size_t b = 0; b < this->num_bodies; b++) {
                double vNorm = std::sqrt(
                        this->bodies[b].node_x_t_n.squaredNorm() + this->bodies[b].node_y_t_n.squaredNorm() +
                        this->bodies[b].node_z_t_n.squaredNorm());
                if (vNorm <= TOL) {
                    vNorm = 1.0;
                }
                double sNorm = this->bodies[b].qk.norm();
                if (sNorm > TOL) {
                    this->bodies[b].DhRx = sNorm / (h * vNorm) * (this->bodies[b].Rx - this->bodies[b].Rvx);
                    this->bodies[b].DhRy = sNorm / (h * vNorm) * (this->bodies[b].Ry - this->bodies[b].Rvy);
                    this->bodies[b].DhRz = sNorm / (h * vNorm) * (this->bodies[b].Rz - this->bodies[b].Rvz);

                    this->bodies[b].tk << this->bodies[b].DhRx, this->bodies[b].DhRy, this->bodies[b].DhRz;
                    this->bodies[b].ok =
                            this->bodies[b].tk.dot(this->bodies[b].qk) / (this->bodies[b].tk.dot(this->bodies[b].tk));
                    this->bodies[b].sk = this->bodies[b].hk + this->bodies[b].ok * this->bodies[b].qk;

                    //check convergence of sk?

                    if (!std::isfinite(this->bodies[b].ok)) {
                        std::cout << "ok is infinite\n";
                    }

                    this->bodies[b].rk = this->bodies[b].qk - this->bodies[b].ok * this->bodies[b].tk;
                } else {
                    this->bodies[b].sk = this->bodies[b].hk;
                    this->bodies[b].rk = this->bodies[b].qk;
                }
            }
            k += 1;
            rhoSum = 0;
            for (size_t b = 0; b < this->num_bodies; b++) {
                rhoSum += this->bodies[b].rk.squaredNorm();
            }
            std::cout << "\rn: " << nIter << " k: " << k <<  " r: " << rhoSum << "      \r" << std::flush;
        }

        for (size_t b = 0; b < this->num_bodies; b++) {
            this->bodies[b].node_x_t_n += this->bodies[b].sk.segment(0,this->bodies[b].n);
            this->bodies[b].node_y_t_n += this->bodies[b].sk.segment(this->bodies[b].n,this->bodies[b].n);
            this->bodies[b].node_z_t_n += this->bodies[b].sk.segment(2*this->bodies[b].n,this->bodies[b].n);

            this->bodies[b].node_x_t_trial = this->bodies[b].node_x_t_n;
            this->bodies[b].node_y_t_trial = this->bodies[b].node_y_t_n;
            this->bodies[b].node_z_t_trial = this->bodies[b].node_z_t_n;

            this->bodies[b].node_mx_t = this->bodies[b].node_x_t_trial.array() * this->bodies[b].node_m.array();
            this->bodies[b].node_my_t = this->bodies[b].node_y_t_trial.array() * this->bodies[b].node_m.array();
            this->bodies[b].node_mz_t = this->bodies[b].node_z_t_trial.array() * this->bodies[b].node_m.array();

            this->bodies[b].node_x_t = this->bodies[b].node_x_t_trial;
            this->bodies[b].node_y_t = this->bodies[b].node_y_t_trial;
            this->bodies[b].node_z_t = this->bodies[b].node_z_t_trial;
        }

        std::cout << "n: " << nIter << " k: " << k <<  " r: " << rhoSum << std::endl;
        nIter += 1;

        //add contact forces
        this->addContactForces();

        //enforce boundary conditions
        this->addBoundaryConditions();
        //calculate L on particles
        this->calculateStrainRate();

        //update particle densities
        this->updateTrialDensity();

        //material stress update
        this->updateTrialStress();

        this->mapTrialStress2Grid();
        //calculate residual for expicit step

        //add contact forces
        this->addContactForces();

        //enforce boundary conditions
        this->addBoundaryConditions();

        for (size_t b=0;b<this->num_bodies;b++) {
            this->bodies[b].node_fx_L = this->bodies[b].node_contact_fx;
            this->bodies[b].node_fy_L = this->bodies[b].node_contact_fy;
            this->bodies[b].node_fz_L = this->bodies[b].node_contact_fz;
        }

        this->calculateImplicitResidual();

        //setup for next iteration
        rhoSum = 0;
        for (size_t b = 0; b < this->num_bodies; b++) {
            this->bodies[b].rk << -this->bodies[b].Rx, -this->bodies[b].Ry, -this->bodies[b].Rz;
            rhoSum += this->bodies[b].rk.squaredNorm();
        }
        //std::cout << "n: " << nIter << " k: " << k <<  " r: " << rhoSum << std::endl;
    } while (rhoSum>rhoTOL && rhoSum<R_MAX); //(this->num_bodies * this->num_nodes * R_TOL));

    //catch error where rhoSum is nan
    if (!std::isfinite(rhoSum) || rhoSum>=R_MAX){
        std::cout << "Error: Residual in mpmStepUSLImplicit() is infinite." << std::endl;
        exit(0);
    }

    return;
}

//*******************************************************************//
//**************************MPM STEP*********************************//
//*******************************************************************//

int job_t::mpmStepUSLExplicit() {
    //forward step
    this->t += this->dt;
    this->stepcount += 1;

    //create particle map
    this->createMappings();

    //map particles to grid
    this->mapParticles2Grid();

    //add contact forces
    this->addContactForces();

    //enforce boundary conditions
    this->addBoundaryConditions();

    //move grid
    this->moveGridExplicit();

    //move particles
    this->moveParticlesExplicit();

    //calculate L on particles
    this->calculateStrainRate();

    //update particle densities
    this->updateDensity();

    //material stress update
    this->updateStress();

    //add boddy forces
    time_varying_loads(this);

    return 1;
}

int job_t::mpmStepUSLExplicit2D() {
    //forward step
    this->t += this->dt;
    this->stepcount += 1;

    //create particle map
    this->createMappings();

    //map particles to grid
    this->mapParticles2Grid();

    //add contact forces
    this->addContactForces2D();

    //enforce boundary conditions
    this->addBoundaryConditions();

    //move grid
    this->moveGridExplicit2D();

    //move particles
    this->moveParticlesExplicit2D();

    //calculate L on particles
    this->calculateStrainRate2D();

    //update particle densities
    this->updateDensity();

    //material stress update
    this->updateStress();

    //add boddy forces
    time_varying_loads(this);

    return 1;
}

int job_t::mpmStepUSLImplicit() {
    //forward step
    this->t += this->dt;
    this->stepcount += 1;

    //create particle map
    this->createMappings();

    //map particles to grid
    this->mapParticles2Grid();

    //add contact forces
    this->addContactForces();

    //enforce boundary conditions
    this->addBoundaryConditions();

    //move grid
    //this->moveGridExplicit();
    //this->moveGridImplicitCG();
    this->moveGridImplicitBiCGSTAB();

    //add contact forces
    this->addContactForces();

    //enforce boundary conditions
    this->addBoundaryConditions();

    //move particles
    //this->moveParticlesExplicit();
    this->moveParticlesImplicit();

    //calculate L on particles
    this->calculateStrainRate();

    //update particle densities
    this->updateDensity();

    //material stress update
    this->updateStress();

    //add boddy forces
    time_varying_loads(this);

    return 1;
}
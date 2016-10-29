//
// Created by aaron on 8/27/16.
// process.cpp
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <algorithm>

#include "tensor.hpp"
#include "process.hpp"
#include "body.hpp"
#include "loading.hpp"

//mass tolerance
#define TOL 5e-11

//hard coded for now. need to change
job_t::job_t():
        use_cpdi(1),
        dt(1e-6),
        t(0.0),
        step_start_time(0.0)
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
        size_t Nx, double hx
) {
    size_t i = node_number % Nx;
    size_t j = (node_number/Nx) % Nx;
    size_t k = node_number / (Nx*Nx);

    *x = i*hx;
    *y = j*hx;
    *z = k*hx;
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

            this->u_dirichlet.resize(numNodes*NODAL_DOF);
            this->u_dirichlet_mask.resize(numNodes*NODAL_DOF);
            this->node_number_override.resize(numNodes*NODAL_DOF);
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
            //std::cout << "{" << x << " " << y << " " << z << "} -> " << idOut << "\n";
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
        this->bodies[b].Sip.setZero();
        this->bodies[b].gradSipX.setZero();
        this->bodies[b].gradSipY.setZero();
        this->bodies[b].gradSipZ.setZero();
        //reset triplets
        this->bodies[b].SipTriplets.clear();
        this->bodies[b].gradSipXTriplets.clear();
        this->bodies[b].gradSipYTriplets.clear();
        this->bodies[b].gradSipZTriplets.clear();
        //loop over particles
        for (size_t p = 0; p < this->bodies[b].p; p++) {
            if((this->bodies[b].particles[p].updateActive(this))==1) {
                this->bodies[b].particles[p].updateCorners(this);
                //check that all corners are in domain
                if (std::count(this->bodies[b].particles[p].corner_elements,
                               std::end(this->bodies[b].particles[p].corner_elements), -1) > 0) {
                    //set corners to particle position
                    this->bodies[b].particles[p].resetCorners(this);
                }

                //use cpdi (or count center 8 times if corners reset)
                for (size_t c=0;c<8;c++) {
                    size_t e = this->bodies[b].particles[p].corner_elements[c];
                    this->bodies[b].elements[e].calculateSipc(&(this->bodies[b]), &(this->bodies[b].particles[p]), c);
                    //std::cout << "b:" << b << " p:" << p << " c:" << c << "\r";
                }
            }
        }
        //build Sipc from triples
        //std::cout << "complete map\n";
        this->bodies[b].Sip.setFromTriplets(this->bodies[b].SipTriplets.begin(),this->bodies[b].SipTriplets.end());
        this->bodies[b].gradSipX.setFromTriplets(this->bodies[b].gradSipXTriplets.begin(),this->bodies[b].gradSipXTriplets.end());
        this->bodies[b].gradSipY.setFromTriplets(this->bodies[b].gradSipYTriplets.begin(),this->bodies[b].gradSipYTriplets.end());
        this->bodies[b].gradSipZ.setFromTriplets(this->bodies[b].gradSipZTriplets.begin(),this->bodies[b].gradSipZTriplets.end());
    }
    return 1;
}


void job_t::mapParticles2Grid() {
    for (size_t b=0; b<this->num_bodies; b++) {
            //use Eigen Map to point to particle array
            size_t numRowsP = this->bodies[b].p;
            size_t numColsP = 1;
            Eigen::Map<Eigen::MatrixXd> p_m(this->bodies[b].particle_m.data(),numRowsP,numColsP);
            Eigen::Map<Eigen::MatrixXd> p_x_t(this->bodies[b].particle_x_t.data(),numRowsP,numColsP);
            Eigen::Map<Eigen::MatrixXd> p_y_t(this->bodies[b].particle_y_t.data(),numRowsP,numColsP);
            Eigen::Map<Eigen::MatrixXd> p_z_t(this->bodies[b].particle_z_t.data(),numRowsP,numColsP);
            Eigen::Map<Eigen::MatrixXd> p_bx(this->bodies[b].particle_bx.data(),numRowsP,numColsP);
            Eigen::Map<Eigen::MatrixXd> p_by(this->bodies[b].particle_by.data(),numRowsP,numColsP);
            Eigen::Map<Eigen::MatrixXd> p_bz(this->bodies[b].particle_bz.data(),numRowsP,numColsP);

            Eigen::MatrixXd p_m_x_t(numRowsP,numColsP);
            Eigen::MatrixXd p_m_y_t(numRowsP,numColsP);
            Eigen::MatrixXd p_m_z_t(numRowsP,numColsP);
            p_m_x_t = p_m.array()*p_x_t.array();
            p_m_y_t = p_m.array()*p_y_t.array();
            p_m_z_t = p_m.array()*p_z_t.array();

            Eigen::MatrixXd p_m_bx(numRowsP,numColsP);
            Eigen::MatrixXd p_m_by(numRowsP,numColsP);
            Eigen::MatrixXd p_m_bz(numRowsP,numColsP);
            p_m_bx = p_m.array()*p_bx.array();
            p_m_by = p_m.array()*p_by.array();
            p_m_bz = p_m.array()*p_bz.array();

            //use Eigen Map to point to node array
            size_t numRowsN = this->bodies[b].n;
            size_t numColsN = 1;
            Eigen::Map<Eigen::MatrixXd> n_m(this->bodies[b].node_m.data(),numRowsN,numColsN);
            Eigen::Map<Eigen::MatrixXd> n_m_x_t(this->bodies[b].node_mx_t.data(),numRowsN,numColsN);
            Eigen::Map<Eigen::MatrixXd> n_m_y_t(this->bodies[b].node_my_t.data(),numRowsN,numColsN);
            Eigen::Map<Eigen::MatrixXd> n_m_z_t(this->bodies[b].node_mz_t.data(),numRowsN,numColsN);
            Eigen::Map<Eigen::MatrixXd> n_fx(this->bodies[b].node_contact_fx.data(),numRowsN,numColsN);
            Eigen::Map<Eigen::MatrixXd> n_fy(this->bodies[b].node_contact_fy.data(),numRowsN,numColsN);
            Eigen::Map<Eigen::MatrixXd> n_fz(this->bodies[b].node_contact_fz.data(),numRowsN,numColsN);
            Eigen::Map<Eigen::MatrixXd> n_nx(this->bodies[b].node_contact_normal_x.data(),numRowsN,numColsN);
            Eigen::Map<Eigen::MatrixXd> n_ny(this->bodies[b].node_contact_normal_y.data(),numRowsN,numColsN);
            Eigen::Map<Eigen::MatrixXd> n_nz(this->bodies[b].node_contact_normal_z.data(),numRowsN,numColsN);

            //use to create dummy pvec
            Eigen::MatrixXd pvec(numRowsP,numColsP);

        //use Sip to map particles to nodes
        n_m = this->bodies[b].Sip*p_m;
        n_m_x_t = this->bodies[b].Sip*p_m_x_t;
        n_m_y_t = this->bodies[b].Sip*p_m_y_t;
        n_m_z_t = this->bodies[b].Sip*p_m_z_t;
        n_fx = this->bodies[b].Sip*p_m_bx; //need to add stress
        n_fy = this->bodies[b].Sip*p_m_by; //need to add stress
        n_fz = this->bodies[b].Sip*p_m_bz; //need to add stress

        //use gradSip to map particles to nodes
        n_nx = this->bodies[b].gradSipX*p_m;
        n_ny = this->bodies[b].gradSipY*p_m;
        n_nz = this->bodies[b].gradSipZ*p_m;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec(i,0) = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[XX];
        }
        n_fx -= this->bodies[b].gradSipX*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec(i,0) = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[XY];
        }
        n_fx -= this->bodies[b].gradSipY*pvec;
        n_fy -= this->bodies[b].gradSipX*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec(i,0) = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[XZ];
        }
        n_fx -= this->bodies[b].gradSipZ*pvec;
        n_fz -= this->bodies[b].gradSipX*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec(i,0) = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[YY];
        }
        n_fy -= this->bodies[b].gradSipY*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec(i,0) = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[YZ];
        }
        n_fy -= this->bodies[b].gradSipZ*pvec;
        n_fz -= this->bodies[b].gradSipY*pvec;

        for (size_t i=0;i<this->bodies[b].p;i++){
            pvec(i,0) = this->bodies[b].particle_v[i] * this->bodies[b].particles[i].T[ZZ];
        }
        n_fz -= this->bodies[b].gradSipZ*pvec;

    }
    return;
}

void job_t::addContactForces(){
    //resolve conflicts between grid velocites
    //implement later 10/21/16
    for (size_t b=0;b<this->num_bodies;b++){
        for (size_t i=0;i<this->bodies[b].n;i++){
            this->bodies[b].node_contact_mx_t[i] = this->bodies[b].node_mx_t[i];
            this->bodies[b].node_contact_my_t[i] = this->bodies[b].node_my_t[i];
            this->bodies[b].node_contact_mz_t[i] = this->bodies[b].node_mz_t[i];

            //the following appear unused
            this->bodies[b].node_contact_x_t[i] = this->bodies[b].node_x_t[i];
            this->bodies[b].node_contact_y_t[i] = this->bodies[b].node_y_t[i];
            this->bodies[b].node_contact_z_t[i] = this->bodies[b].node_z_t[i];

            this->bodies[b].node_contact_fx[i] = this->bodies[b].node_fx[i];
            this->bodies[b].node_contact_fy[i] = this->bodies[b].node_fy[i];
            this->bodies[b].node_contact_fz[i] = this->bodies[b].node_fz[i];
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
    for (size_t b = 0; b < this->num_bodies; b++){
        for (size_t i = 0; i < this->num_nodes; i++){
            double m = this->bodies[b].nodes[i].m[0];
            if (m > TOL) {
                this->bodies[b].nodes[i].contact_mx_t[0] += this->dt * this->bodies[b].nodes[i].contact_fx[0];
                this->bodies[b].nodes[i].contact_my_t[0] += this->dt * this->bodies[b].nodes[i].contact_fy[0];
                this->bodies[b].nodes[i].contact_mz_t[0] += this->dt * this->bodies[b].nodes[i].contact_fz[0];

                this->bodies[b].nodes[i].contact_x_t[0] = this->bodies[b].nodes[i].contact_mx_t[0] / m;
                this->bodies[b].nodes[i].contact_y_t[0] = this->bodies[b].nodes[i].contact_my_t[0] / m;
                this->bodies[b].nodes[i].contact_z_t[0] = this->bodies[b].nodes[i].contact_mz_t[0] / m;
            } else {
                this->bodies[b].nodes[i].contact_x_t[0] = 0;
                this->bodies[b].nodes[i].contact_y_t[0] = 0;
                this->bodies[b].nodes[i].contact_z_t[0] = 0;
            }

        }
    }
    return;
}

void job_t::moveParticlesExplicit(){
    for (size_t b=0; b<this->num_bodies; b++) {
        //use Eigen Map to point to particle array
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;
        Eigen::Map<Eigen::MatrixXd> p_x(this->bodies[b].particle_x_t.data(), numRowsP, numColsP);
        Eigen::Map<Eigen::MatrixXd> p_y(this->bodies[b].particle_y_t.data(), numRowsP, numColsP);
        Eigen::Map<Eigen::MatrixXd> p_z(this->bodies[b].particle_z_t.data(), numRowsP, numColsP);
        Eigen::Map<Eigen::MatrixXd> p_x_t(this->bodies[b].particle_x_t.data(), numRowsP, numColsP);
        Eigen::Map<Eigen::MatrixXd> p_y_t(this->bodies[b].particle_y_t.data(), numRowsP, numColsP);
        Eigen::Map<Eigen::MatrixXd> p_z_t(this->bodies[b].particle_z_t.data(), numRowsP, numColsP);
        Eigen::Map<Eigen::MatrixXd> p_ux(this->bodies[b].particle_ux.data(), numRowsP, numColsP);
        Eigen::Map<Eigen::MatrixXd> p_uy(this->bodies[b].particle_uy.data(), numRowsP, numColsP);
        Eigen::Map<Eigen::MatrixXd> p_uz(this->bodies[b].particle_uz.data(), numRowsP, numColsP);

        //use Eigen Map to point to node array
        size_t numRowsN = this->bodies[b].n;
        size_t numColsN = 1;
        Eigen::Map<Eigen::MatrixXd> n_ux(this->bodies[b].node_ux.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_uy(this->bodies[b].node_uy.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_uz(this->bodies[b].node_uz.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_dx_t(this->bodies[b].node_diff_x_t.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_dy_t(this->bodies[b].node_diff_y_t.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_dz_t(this->bodies[b].node_diff_z_t.data(), numRowsN, numColsN);

        //use to create dummy pvec
        Eigen::MatrixXd pvec(numRowsP, numColsP);

        for (size_t i=0;i<this->bodies[b].n;i++){
            double m = this->bodies[b].nodes[i].m[0];
            if (m!=0){
                this->bodies[b].nodes[i].ux[0] = this->dt * this->bodies[b].nodes[i].contact_mx_t[0] / m;
                this->bodies[b].nodes[i].uy[0] = this->dt * this->bodies[b].nodes[i].contact_my_t[0] / m;
                this->bodies[b].nodes[i].uz[0] = this->dt * this->bodies[b].nodes[i].contact_mz_t[0] / m;

                this->bodies[b].nodes[i].diff_x_t[0] = this->dt * this->bodies[b].nodes[i].contact_fx[0] / m;
                this->bodies[b].nodes[i].diff_y_t[0] = this->dt * this->bodies[b].nodes[i].contact_fy[0] / m;
                this->bodies[b].nodes[i].diff_y_t[0] = this->dt * this->bodies[b].nodes[i].contact_fz[0] / m;
            } else {
                this->bodies[b].nodes[i].ux[0] = 0;
                this->bodies[b].nodes[i].uy[0] = 0;
                this->bodies[b].nodes[i].uz[0] = 0;

                this->bodies[b].nodes[i].diff_x_t[0] = 0;
                this->bodies[b].nodes[i].diff_y_t[0] = 0;
                this->bodies[b].nodes[i].diff_y_t[0] = 0;
            }
        }

        //map back to particles using S transpose
        p_x += this->bodies[b].Sip.transpose()*n_ux;
        p_y += this->bodies[b].Sip.transpose()*n_uy;
        p_z += this->bodies[b].Sip.transpose()*n_uz;

        p_ux += this->bodies[b].Sip.transpose()*n_ux;
        p_uy += this->bodies[b].Sip.transpose()*n_uy;
        p_uz += this->bodies[b].Sip.transpose()*n_uz;

        p_x_t += this->bodies[b].Sip.transpose()*n_dx_t;
        p_y_t += this->bodies[b].Sip.transpose()*n_dy_t;
        p_z_t += this->bodies[b].Sip.transpose()*n_dz_t;

    }
    return;
}

void job_t::calculateStrainRate() {
    for (size_t b=0; b<this->num_bodies; b++) {
        //map nodal velocities with Eigen
        size_t numRowsN = this->bodies[b].n;
        size_t numColsN = 1;
        Eigen::Map<Eigen::MatrixXd> n_m(this->bodies[b].node_m.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_c_mx_t(this->bodies[b].node_contact_mx_t.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_c_my_t(this->bodies[b].node_contact_my_t.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_c_mz_t(this->bodies[b].node_contact_mz_t.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_c_x_t(this->bodies[b].node_contact_x_t.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_c_y_t(this->bodies[b].node_contact_y_t.data(), numRowsN, numColsN);
        Eigen::Map<Eigen::MatrixXd> n_c_z_t(this->bodies[b].node_contact_z_t.data(), numRowsN, numColsN);

        //nodal velocities
        n_c_x_t = n_c_mx_t.array()/n_m.array();
        n_c_y_t = n_c_my_t.array()/n_m.array();
        n_c_z_t = n_c_mz_t.array()/n_m.array();

        //use to create dummy pvec
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;
        Eigen::MatrixXd pvec(numRowsP,numColsP);

        //calculate particle[i].L[9]
        pvec = this->bodies[b].gradSipX.transpose()*n_c_x_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[XX] = pvec(i,0);
        }

        pvec = this->bodies[b].gradSipY.transpose()*n_c_x_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[XY] = pvec(i,0);
        }

        pvec = this->bodies[b].gradSipZ.transpose()*n_c_x_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[XZ] = pvec(i,0);
        }

        pvec = this->bodies[b].gradSipX.transpose()*n_c_y_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[YX] = pvec(i,0);
        }

        pvec = this->bodies[b].gradSipY.transpose()*n_c_y_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[YY] = pvec(i,0);
        }

        pvec = this->bodies[b].gradSipZ.transpose()*n_c_y_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[YZ] = pvec(i,0);
        }

        pvec = this->bodies[b].gradSipX.transpose()*n_c_z_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[ZX] = pvec(i,0);
        }

        pvec = this->bodies[b].gradSipY.transpose()*n_c_z_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[ZY] = pvec(i,0);
        }

        pvec = this->bodies[b].gradSipZ.transpose()*n_c_z_t;
        for (size_t i=0;i<this->bodies[b].p;i++){
            this->bodies[b].particles[i].L[ZZ] = pvec(i,0);
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

void job_t::updateStress(){
    //calculate stress
    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].material.calculate_stress(&(this->bodies[b]),this->dt);
    }
    return;
}

//*******************************************************************//
//**************************MPM STEP*********************************//
//*******************************************************************//

int job_t::mpmStepUSLExplicit() {
    //forward step
    this->t += this->dt;

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
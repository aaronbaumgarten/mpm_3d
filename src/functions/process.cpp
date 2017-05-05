//
// Created by aaron on 8/27/16.
// process.cpp
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <Eigen/Core>

#include "tensor.hpp"
#include "process.hpp"
#include "body.hpp"
#include "loading.hpp"
#include "node.hpp"
#include "element.hpp"

//mass tolerance
//#define TOL 0//5e-11
//contact friction
#define MU_F 0.4
//squared norm error tolerance
//#define R_TOL 1e-5
#define R_MAX 1e32

//hard coded for now. need to change
job_t::job_t():
        use_cpdi(1), //default cpdi
        use_3d(1), //default 3d
        use_implicit(0), //default explicit
        use_smoothing(0),
        dt(1e-3),
        dt_base(dt),
        dt_minimum(1e-6),
        t(0.0),
        step_start_time(0.0),
        stepcount(0),
        newtonTOL(1e-5),
        linearStepSize(1e-6),
        num_contacts(0),
        num_bodies(0),
        num_particles(0),
        num_nodes(0),
        num_elements(0)
{
    std::cout << "Job created.\n";
    boundary = Boundary();
    elements = Elements();
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
        return 0;
    }
    if (j>jmax || j<0){
        return 0;
    }
    if (k>kmax || k<0){
        return 0;
    }
    return i + j*imax + k*jmax*imax;
}

int job_t::importNodesandParticles(std::string nfilename, std::string pfilename){
    //only reads particle file for n bodies
    this->particleFile = pfilename;
    this->nodeFile = nfilename;

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
        size_t numBodies = 0;

        //open grid file and parse
        std::ifstream fin(nfilename);
        if (fin.is_open()) {
            std::string line;
            if (std::getline(fin, line)) {
                std::stringstream ss(line);
                if (!(ss >> numLinearNodesX >> numLinearNodesY >> numLinearNodesZ)) {
                    std::cout << "Cannot parse grid file: " << nfilename << "\n";
                    return 0;
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
                if (!(ss >> this->Lx >> this->Ly >> this->Lz)) {
                    std::cout << "Cannot parse grid file: " << nfilename << "\n";
                    return 0;
                }
                this->hx = this->Lx / (numLinearNodesX - 1);
                this->hy = this->Ly / (numLinearNodesY - 1);
                this->hz = this->Lz / (numLinearNodesZ - 1);
            }
        } else {
            std::cout << "Cannot parse grid file: " << nfilename << "\n";
            return 0;
        }

        fin.close();

        //open particle file
        fin.open(pfilename);
        if (fin.is_open()) {
            std::string line;

            //parse header
            if (getline(fin, line)) {
                std::stringstream ss(line);
                if (!(ss >> numBodies)) {
                    std::cout << "Cannot parse particle file: " << pfilename << "\n";
                    std::cout << "Expected body count at file header." << "\n";
                    return 0;
                }
            }
            //create body objects
            double np = 0;
            for (size_t b=0;b<numBodies;b++){
                if (getline(fin, line)) {
                    std::stringstream ss(line);
                    if (!(ss >> np)) {
                        std::cout << "Cannot parse particle file: " << pfilename << "\n";
                        std::cout << "Expected particle counts after body count at file header." << "\n";
                        return 0;
                    }
                    this->bodies.push_back(Body(numNodes, np, numElements, b));
                    numParticles += np;
                }
            }

            this->num_particles = numParticles;
            this->num_bodies = numBodies;
            std::cout << "Bodies created (" << numBodies << ").\n";

            //assign particle to bodies for each particle in particle file
            while (getline(fin, line)) {
                std::stringstream ss(line);
                double b, pID, m, v, x, y, z, x_t, y_t, z_t;
                size_t idOut;

                if (!(ss >> b >> pID >> m >> v >> x >> y >> z >> x_t >> y_t >> z_t)) {
                    std::cout << "Cannot parse particle file: " << pfilename << "\n";
                    std::cout << "Line parsing failed: " << ss.str() << "\n";
                    return 0;
                }

                assert(b<this->num_bodies);
                assert(pID<this->bodies[b].p);

                this->bodies[b].particles.addParticle(m, v, x, y, z, x_t, y_t, z_t, (size_t)pID); //0-index particles
            }
        } else {
            std::cout << "Cannot parse particle file: " << pfilename << "\n";
            std::cout << "Unable to open file." << "\n";
            return 0;
        }
        fin.close();

        std::cout << "Particles created (" << numParticles << ").\n";

        //assign nodes to bodies
        for (size_t i = 0; i < numBodies; i++) {
            for (size_t nn = 0; nn < numNodes; nn++) {
                double x, y, z;
                job_t::node_number_to_coords(&x, &y, &z, nn, numLinearNodesX, numLinearNodesY, numLinearNodesZ, this->hx, this->hy, this->hz);
                this->bodies[i].nodes.addNode(x, y, z, nn);
            }

        }

        std::cout << "Nodes created (" << numNodes << ").\n";

        //create element object (using cube configuration for now)
        this->elements.resizeElements(numElements,8);

        //assign elements to bodies and nodes to elements
        for (size_t ne = 0; ne < numElements; ne++) {
            Eigen::VectorXi nodeIDs(8);
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

            this->elements.addElement(nodeIDs, ne);
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
                return 0;
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
                return 0;
            }
            this->hx = this->Lx / (numLinearNodesX - 1);
            this->hy = this->Ly / (numLinearNodesY - 1);
            this->hz = this->Lz;
        }
    } else {
        std::cout << "Cannot parse grid file: " << nfilename << "\n";
        return 0;
    }

    fin.close();

    //open particle file
    fin.open(pfilename);
    if (fin.is_open()) {
        std::string line;

        //parse header
        if (getline(fin, line)) {
            std::stringstream ss(line);
            if (!(ss >> numBodies)) {
                std::cout << "Cannot parse particle file: " << pfilename << "\n";
                std::cout << "Expected body count at file header." << "\n";
                return 0;
            }
        }
        //create body objects
        double np = 0;
        for (size_t b=0;b<numBodies;b++){
            if (getline(fin, line)) {
                std::stringstream ss(line);
                if (!(ss >> np)) {
                    std::cout << "Cannot parse particle file: " << pfilename << "\n";
                    std::cout << "Expected particle counts after body count at file header." << "\n";
                    return 0;
                }
                this->bodies.push_back(Body(numNodes, np, numElements, b));
                numParticles += np;
            }
        }

        this->num_particles = numParticles;
        this->num_bodies = numBodies;
        std::cout << "Bodies created (" << numBodies << ").\n";

        //assign particle to bodies for each particle in particle file
        while (getline(fin, line)) {
            std::stringstream ss(line);
            double b, pID, m, v, x, y, z, x_t, y_t, z_t;
            size_t idOut;

            if (!(ss >> b >> pID >> m >> v >> x >> y >> z >> x_t >> y_t >> z_t)) {
                std::cout << "Cannot parse particle file: " << pfilename << "\n";
                std::cout << "Line parsing failed: " << ss.str() << "\n";
                return 0;
            }

            assert(b<this->num_bodies);
            assert(pID<this->bodies[b].p);

            this->bodies[b].particles.addParticle(m, v, x, y, z, x_t, y_t, z_t, (size_t)pID); //0-index particles
        }
    } else {
        std::cout << "Cannot parse particle file: " << pfilename << "\n";
        std::cout << "Unable to open file.\n";
        return 0;
    }
    fin.close();

    std::cout << "Particles created (" << numParticles << ").\n";

    //assign nodes to bodies
    for (size_t i=0; i<numBodies; i++){
        for (size_t nn=0; nn<numNodes; nn++){
            double x, y, z;
            job_t::node_number_to_coords(&x,&y,&z,nn,numLinearNodesX,numLinearNodesY,numLinearNodesZ,this->hx,this->hy,this->hz); //ok for 2d as nodes count in x,y first
            this->bodies[i].nodes.addNode(x,y,z,nn);
        }

    }

    std::cout << "Nodes created (" << numNodes << ").\n";

    //create element object (using cube configuration for now)
    this->elements.resizeElements(numElements,8);

    //assign elements to bodies and nodes to elements
    for(size_t ne=0; ne<numElements; ne++){
        Eigen::VectorXi nodeIDs(8);
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

        this->elements.addElement(nodeIDs,ne);
    }

    std::cout << "Elements created (" << numElements << ").\n";
    return 1;
}

int job_t::assignDefaultMaterials() {
    std::string filename = "isolin.so";
    std::vector<double> matprops = {1e9,0.3};
    for (size_t i=0;i<this->num_bodies;i++){
        this->bodies[i].defineMaterial(filename,"./",matprops,std::vector<int>());
    }
    std::cout << "Materials assigned (" << this->num_bodies << ").\n";
    return 1;
}

int job_t::assignMaterial(std::string filename, std::string filepath, size_t id, std::vector<double> fp64props, std::vector<int> intprops) {
    this->bodies[id].defineMaterial(filename,filepath,fp64props,intprops);
    std::cout << "Material reassigned [" << id << "].\n";
    return 1;
}

int job_t::assignDefaultBoundaryConditions(){
    //default value
    std::string filename = "boxBC.so";
    this->boundary.setBoundary(filename,"./",std::vector<double>(),std::vector<int>());
    this->boundary.bc_init(this);

    std::cout << "Boundary Conditions assigned (1).\n";

    return 1;
}

int job_t::assignBoundaryConditions(std::string bcFile, std::string bcPath, std::vector<double> bcfp64, std::vector<int> bcint){

    this->boundary.setBoundary(bcFile,bcPath,bcfp64,bcint);
    this->boundary.bc_init(this);

    std::cout << "Boundary Conditions assigned (1).\n";

    return 1;
}

int job_t::assignDefaultContacts() {
    std::string filename;
    if (this->use_3d==1){
        filename = "huang.so";
    } else {
        filename = "huang2d.so";
    }
    this->num_contacts = this->num_bodies*(this->num_bodies-1)/2;
    this->contacts.resize(this->num_contacts);
    std::vector<int> bodyIDs = {0,1};
    std::vector<double> contactProps = {0.4};
    int b1 = 0;
    int b2 = 0;
    for(size_t i=0;i<this->num_contacts;i++){
        b2++;
        if (b2 > this->num_bodies){
            b1++;
            b2 = b1+1;
        }
        bodyIDs = {b1,b2};
        this->contacts[i].setContact(filename,"./",i,bodyIDs,contactProps,std::vector<int>());
        this->contacts[i].contact_init(this,i);
    }
    std::cout << "Contact Rules assigned (" << this->num_contacts << ")\n";
    return 1;
}

int job_t::assignContact(std::string filename, std::string filepath, std::vector<int> bodyIDs, std::vector<double> fp64props, std::vector<int> intprops) {
    size_t id = this->num_contacts;
    this->num_contacts += 1;
    this->contacts.push_back(Contact());
    this->contacts[id].setContact(filename,filepath,id,bodyIDs,fp64props,intprops);
    this->contacts[id].contact_init(this,id);

    //fix adjusted contacts
    for (size_t c=0; c<id; c++){
        this->contacts[c].fixFunctionPointers();
    }

    std::cout << "Contact Rule reassigned [" << id << "]\n";
    return 1;
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
            if((this->bodies[b].particles.updateActive(this,p))==1) {
                this->bodies[b].particles.updateCorners(this,p);
                //check that all corners are in domain
                int useExtent = 1;
                if (this->use_cpdi==0 || this->bodies[b].particles.corner_elements.row(p).minCoeff() == -1){
                    //set corners to particle position
                    this->bodies[b].particles.resetCorners(this,p);
                    useExtent = 0;
                }

                //use cpdi (or count center 8 times if corners reset)
                for (size_t c=0;c<8;c++) {
                    size_t e = this->bodies[b].particles.corner_elements(p,c);
                    this->elements.calculatePhic(&(this->bodies[b]), e, p, c, useExtent);
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
    this->boundary.generate_dirichlet_bcs(this);
    for (size_t b=0; b<this->num_bodies; b++) {
        //use Eigen Map to point to particle array
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;
        
        Eigen::MatrixXd p_m_x_t(numRowsP,numColsP);
        Eigen::MatrixXd p_m_y_t(numRowsP,numColsP);
        Eigen::MatrixXd p_m_z_t(numRowsP,numColsP);
        p_m_x_t = this->bodies[b].particles.m.array()*this->bodies[b].particles.x_t.array();
        p_m_y_t = this->bodies[b].particles.m.array()*this->bodies[b].particles.y_t.array();
        p_m_z_t = this->bodies[b].particles.m.array()*this->bodies[b].particles.z_t.array();

        Eigen::MatrixXd p_m_bx(numRowsP,numColsP);
        Eigen::MatrixXd p_m_by(numRowsP,numColsP);
        Eigen::MatrixXd p_m_bz(numRowsP,numColsP);
        p_m_bx = this->bodies[b].particles.m.array()*this->bodies[b].particles.bx.array();
        p_m_by = this->bodies[b].particles.m.array()*this->bodies[b].particles.by.array();
        p_m_bz = this->bodies[b].particles.m.array()*this->bodies[b].particles.bz.array();

        //use Eigen Map to point to node array
        size_t numRowsN = this->bodies[b].n;
        size_t numColsN = 1;
        
        //use to create dummy pvec and ones
        Eigen::VectorXd pvec(numRowsP);

        //use Phi to map particles to nodes
        this->bodies[b].nodes.m = this->bodies[b].Phi*this->bodies[b].particles.m;
        this->bodies[b].nodes.mx_t = this->bodies[b].Phi*p_m_x_t;
        this->bodies[b].nodes.my_t = this->bodies[b].Phi*p_m_y_t;
        this->bodies[b].nodes.mz_t = this->bodies[b].Phi*p_m_z_t;
        this->bodies[b].nodes.fx = this->bodies[b].Phi*p_m_bx; //need to add stress
        this->bodies[b].nodes.fy = this->bodies[b].Phi*p_m_by; //need to add stress
        this->bodies[b].nodes.fz = this->bodies[b].Phi*p_m_bz; //need to add stress

        //use gradPhi to map particles to nodes
        this->bodies[b].nodes.contact_normal_x = this->bodies[b].gradPhiX*Eigen::MatrixXd::Constant(numRowsP,numColsP,1);//this->bodies[b].particle_m;
        this->bodies[b].nodes.contact_normal_y = this->bodies[b].gradPhiY*Eigen::MatrixXd::Constant(numRowsP,numColsP,1);//this->bodies[b].particle_m;
        this->bodies[b].nodes.contact_normal_z = this->bodies[b].gradPhiZ*Eigen::MatrixXd::Constant(numRowsP,numColsP,1);//this->bodies[b].particle_m;

        /*Eigen::VectorXd normMag(numRowsN);
        normMag = this->bodies[b].node_contact_normal_x.array().square()
                + this->bodies[b].node_contact_normal_y.array().square()
                + this->bodies[b].node_contact_normal_z.array().square();

        this->bodies[b].node_contact_normal_x = this->bodies[b].node_contact_normal_x.array()/normMag.array().sqrt();
        this->bodies[b].node_contact_normal_y = this->bodies[b].node_contact_normal_y.array()/normMag.array().sqrt();
        this->bodies[b].node_contact_normal_z = this->bodies[b].node_contact_normal_z.array()/normMag.array().sqrt();
         */

         //a poor smoothing technique which better conserves particle histories
        /*Eigen::VectorXd trTOLD = this->bodies[b].particles.T.col(XX) + this->bodies[b].particles.T.col(YY) + this->bodies[b].particles.T.col(ZZ);
        Eigen::VectorXd trT = trTOLD;
        if (this->use_smoothing == 1){
            //smoothing [see Mast et. al. 2012 for kinematic locking solution]
            Eigen::VectorXd alpha(this->num_nodes);
            pvec = trTOLD.array() * this->bodies[b].particles.m.array();
            alpha = (this->bodies[b].Phi * pvec).array();

            alpha = alpha.array()/this->bodies[b].nodes.m.array();

            for (size_t i=0;i<this->num_nodes;i++){
                if (this->bodies[b].nodes.m[i] == 0){
                    alpha[i] = 0;
                }
            }
            trT = this->bodies[b].Phi.transpose() * alpha;
        }*/

        //pvec = this->bodies[b].particles.v.array() * (this->bodies[b].particles.T.col(XX) + (trT - trTOLD)/3.0).array();
        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.T.col(XX).array();
        this->bodies[b].nodes.fx -= this->bodies[b].gradPhiX*pvec;

        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.T.col(XY).array();
        this->bodies[b].nodes.fx -= this->bodies[b].gradPhiY*pvec;
        this->bodies[b].nodes.fy -= this->bodies[b].gradPhiX*pvec;

        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.T.col(XZ).array();
        this->bodies[b].nodes.fx -= this->bodies[b].gradPhiZ*pvec;
        this->bodies[b].nodes.fz -= this->bodies[b].gradPhiX*pvec;

        //pvec = this->bodies[b].particles.v.array() * (this->bodies[b].particles.T.col(YY) + (trT - trTOLD)/3.0).array();
        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.T.col(YY).array();
        this->bodies[b].nodes.fy -= this->bodies[b].gradPhiY*pvec;

        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.T.col(YZ).array();
        this->bodies[b].nodes.fy -= this->bodies[b].gradPhiZ*pvec;
        this->bodies[b].nodes.fz -= this->bodies[b].gradPhiY*pvec;

        //pvec = this->bodies[b].particles.v.array() * (this->bodies[b].particles.T.col(ZZ) + (trT - trTOLD)/3.0).array();
        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.T.col(ZZ).array();
        this->bodies[b].nodes.fz -= this->bodies[b].gradPhiZ*pvec;

    }
    return;
}

void job_t::addContactForces(){
    //if (this->use_3d==1) {
    for (size_t b = 0; b < this->num_bodies; b++) {
        this->bodies[b].nodes.contact_mx_t = this->bodies[b].nodes.mx_t;
        this->bodies[b].nodes.contact_my_t = this->bodies[b].nodes.my_t;
        this->bodies[b].nodes.contact_mz_t = this->bodies[b].nodes.mz_t;

        //the following appear unused
        this->bodies[b].nodes.contact_x_t = this->bodies[b].nodes.mx_t.array()/this->bodies[b].nodes.m.array();
        this->bodies[b].nodes.contact_y_t = this->bodies[b].nodes.my_t.array()/this->bodies[b].nodes.m.array();
        this->bodies[b].nodes.contact_z_t = this->bodies[b].nodes.mz_t.array()/this->bodies[b].nodes.m.array();

        for (size_t i=0;i<this->num_nodes;i++){
            if (this->bodies[b].nodes.m[i] == 0){
                this->bodies[b].nodes.contact_x_t[i] = 0;
                this->bodies[b].nodes.contact_y_t[i] = 0;
                this->bodies[b].nodes.contact_z_t[i] = 0;
            }
        }

        this->bodies[b].nodes.contact_fx = this->bodies[b].nodes.fx;
        this->bodies[b].nodes.contact_fy = this->bodies[b].nodes.fy;
        this->bodies[b].nodes.contact_fz = this->bodies[b].nodes.fz;
    }
    for (size_t c=0;c<this->num_contacts;c++){
        this->contacts[c].resolve_contact(this,c);
    }
    /*
    if (this->num_bodies > 1) {
        //look for contacts if there are two bodies
        for (size_t i = 0; i < this->num_nodes; i++) {
            //test every node for contact
            if (this->bodies[0].nodes.m[i] > TOL && this->bodies[1].nodes.m[i] > TOL) {
                //use normal from body 1
                Eigen::Vector3d n1i;
                n1i << this->bodies[0].nodes.contact_normal_x[i],
                        this->bodies[0].nodes.contact_normal_y[i],
                        this->bodies[0].nodes.contact_normal_z[i];
                //enforce unit length
                n1i /= sqrt(n1i.dot(n1i));

                //determine 'center of mass' velocity
                double m1 = this->bodies[0].nodes.m[i];
                double m2 = this->bodies[1].nodes.m[i];
                Eigen::Vector3d mv1i;
                Eigen::Vector3d mv2i;
                Eigen::Vector3d vCMi;
                mv1i << this->bodies[0].nodes.contact_mx_t[i],
                        this->bodies[0].nodes.contact_my_t[i],
                        this->bodies[0].nodes.contact_mz_t[i];
                mv2i << this->bodies[1].nodes.contact_mx_t[i],
                        this->bodies[1].nodes.contact_my_t[i],
                        this->bodies[1].nodes.contact_mz_t[i];
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
                this->bodies[0].nodes.contact_fx[i] = fcti[0];
                this->bodies[0].nodes.contact_fy[i] = fcti[1];
                this->bodies[0].nodes.contact_fz[i] = fcti[2];

                this->bodies[1].nodes.contact_fx[i] = -fcti[0];
                this->bodies[1].nodes.contact_fy[i] = -fcti[1];
                this->bodies[1].nodes.contact_fz[i] = -fcti[2];

                if (this->use_implicit==0) {
                    //adjust nodal velocities for non-penetration
                    mv1i = mv1i - n1i.dot(mv1i - m1 * vCMi) * n1i;
                    mv2i = mv2i - n1i.dot(mv2i - m2 * vCMi) * n1i;

                    this->bodies[0].nodes.contact_mx_t[i] = mv1i[0];
                    this->bodies[0].nodes.contact_my_t[i] = mv1i[1];
                    this->bodies[0].nodes.contact_mz_t[i] = mv1i[2];

                    this->bodies[0].nodes.contact_x_t[i] = mv1i[0] / m1;
                    this->bodies[0].nodes.contact_y_t[i] = mv1i[1] / m1;
                    this->bodies[0].nodes.contact_z_t[i] = mv1i[2] / m1;

                    this->bodies[1].nodes.contact_mx_t[i] = mv2i[0];
                    this->bodies[1].nodes.contact_my_t[i] = mv2i[1];
                    this->bodies[1].nodes.contact_mz_t[i] = mv2i[2];

                    this->bodies[1].nodes.contact_x_t[i] = mv2i[0] / m2;
                    this->bodies[1].nodes.contact_y_t[i] = mv2i[1] / m2;
                    this->bodies[1].nodes.contact_z_t[i] = mv2i[2] / m2;
                }
            }
        }
    }*/

    return;
    //} else {
    //    return this->addContactForces2D();
    //}
}

void job_t::addContactForces2D(){
    //resolve conflicts between grid velocites
    //implement later 10/21/16
    for (size_t b=0;b<this->num_bodies;b++) {
        this->bodies[b].nodes.contact_mx_t = this->bodies[b].nodes.mx_t;
        this->bodies[b].nodes.contact_my_t = this->bodies[b].nodes.my_t;
        //this->bodies[b].nodes.contact_mz_t = this->bodies[b].nodes.mz_t;

        //the following appear unused
        this->bodies[b].nodes.contact_x_t = this->bodies[b].nodes.mx_t.array()/this->bodies[b].nodes.m.array();
        this->bodies[b].nodes.contact_y_t = this->bodies[b].nodes.my_t.array()/this->bodies[b].nodes.m.array();
        //this->bodies[b].nodes.contact_z_t = this->bodies[b].nodes.mz_t.array()/this->bodies[b].nodes.m.array();
        for (size_t i=0;i<this->num_nodes;i++){
            if (this->bodies[b].nodes.m[i] == 0){
                this->bodies[b].nodes.contact_x_t[i] = 0;
                this->bodies[b].nodes.contact_y_t[i] = 0;
            }
        }

        this->bodies[b].nodes.contact_fx = this->bodies[b].nodes.fx;
        this->bodies[b].nodes.contact_fy = this->bodies[b].nodes.fy;
        //this->bodies[b].nodes.contact_fz = this->bodies[b].nodes.fz;
    }

    for (size_t c=0;c<this->num_contacts;c++){
        this->contacts[c].resolve_contact(this,c);
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
            double m = this->bodies[b].nodes.m[i];
            if (m > TOL) {
                this->bodies[b].nodes.contact_mx_t[i] += this->dt * this->bodies[b].nodes.contact_fx[i];
                this->bodies[b].nodes.contact_my_t[i] += this->dt * this->bodies[b].nodes.contact_fy[i];
                this->bodies[b].nodes.contact_mz_t[i] += this->dt * this->bodies[b].nodes.contact_fz[i];

                this->bodies[b].nodes.contact_x_t[i] = this->bodies[b].nodes.contact_mx_t[i] / m;
                this->bodies[b].nodes.contact_y_t[i] = this->bodies[b].nodes.contact_my_t[i] / m;
                this->bodies[b].nodes.contact_z_t[i] = this->bodies[b].nodes.contact_mz_t[i] / m;
            } else {
                this->bodies[b].nodes.contact_mx_t[i] = 0;
                this->bodies[b].nodes.contact_my_t[i] = 0;
                this->bodies[b].nodes.contact_mz_t[i] = 0;

                this->bodies[b].nodes.contact_x_t[i] = 0;
                this->bodies[b].nodes.contact_y_t[i] = 0;
                this->bodies[b].nodes.contact_z_t[i] = 0;
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
            double m = this->bodies[b].nodes.m[i];
            if (m > TOL) {
                this->bodies[b].nodes.contact_mx_t[i] += this->dt * this->bodies[b].nodes.contact_fx[i];
                this->bodies[b].nodes.contact_my_t[i] += this->dt * this->bodies[b].nodes.contact_fy[i];
                this->bodies[b].nodes.contact_mz_t[i] = 0;//this->dt * this->bodies[b].nodes.contact_fz[i];

                this->bodies[b].nodes.contact_x_t[i] = this->bodies[b].nodes.contact_mx_t[i] / m;
                this->bodies[b].nodes.contact_y_t[i] = this->bodies[b].nodes.contact_my_t[i] / m;
                this->bodies[b].nodes.contact_z_t[i] = 0;//this->bodies[b].nodes.contact_mz_t[i] / m;
            } else {
                this->bodies[b].nodes.contact_mx_t[i] = 0;
                this->bodies[b].nodes.contact_my_t[i] = 0;
                this->bodies[b].nodes.contact_mz_t[i] = 0;

                this->bodies[b].nodes.contact_x_t[i] = 0;
                this->bodies[b].nodes.contact_y_t[i] = 0;
                this->bodies[b].nodes.contact_z_t[i] = 0;
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
            double m = this->bodies[b].nodes.m[i];
            if (m != 0) {
                this->bodies[b].nodes.ux[i] = this->dt * this->bodies[b].nodes.contact_mx_t[i] / m;
                this->bodies[b].nodes.uy[i] = this->dt * this->bodies[b].nodes.contact_my_t[i] / m;
                this->bodies[b].nodes.uz[i] = this->dt * this->bodies[b].nodes.contact_mz_t[i] / m;

                this->bodies[b].nodes.diff_x_t[i] = this->dt * this->bodies[b].nodes.contact_fx[i] / m;
                this->bodies[b].nodes.diff_y_t[i] = this->dt * this->bodies[b].nodes.contact_fy[i] / m;
                this->bodies[b].nodes.diff_z_t[i] = this->dt * this->bodies[b].nodes.contact_fz[i] / m;
            } else {
                this->bodies[b].nodes.ux[i] = 0;
                this->bodies[b].nodes.uy[i] = 0;
                this->bodies[b].nodes.uz[i] = 0;

                this->bodies[b].nodes.diff_x_t[i] = 0;
                this->bodies[b].nodes.diff_y_t[i] = 0;
                this->bodies[b].nodes.diff_z_t[i] = 0;
            }
        }

        //map back to particles using S transpose
        this->bodies[b].particles.x += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.ux;
        this->bodies[b].particles.y += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.uy;
        this->bodies[b].particles.z += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.uz;

        this->bodies[b].particles.ux += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.ux;
        this->bodies[b].particles.uy += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.uy;
        this->bodies[b].particles.uz += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.uz;

        this->bodies[b].particles.x_t += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.diff_x_t;
        this->bodies[b].particles.y_t += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.diff_y_t;
        this->bodies[b].particles.z_t += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.diff_z_t;

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
            double m = this->bodies[b].nodes.m[i];
            if (m != 0) {
                this->bodies[b].nodes.ux[i] = this->dt * (this->bodies[b].nodes.mx_t_k[i] + this->bodies[b].nodes.contact_mx_t[i]) / (2*m);
                this->bodies[b].nodes.uy[i] = this->dt * (this->bodies[b].nodes.my_t_k[i] + this->bodies[b].nodes.contact_my_t[i]) / (2*m);
                this->bodies[b].nodes.uz[i] = this->dt * (this->bodies[b].nodes.mz_t_k[i] + this->bodies[b].nodes.contact_mz_t[i]) / (2*m);

                /*this->bodies[b].nodes.ux[i] = this->dt * this->bodies[b].nodes.contact_mx_t[i] / m;
                this->bodies[b].nodes.uy[i] = this->dt * this->bodies[b].nodes.contact_my_t[i] / m;
                this->bodies[b].nodes.uz[i] = this->dt * this->bodies[b].nodes.contact_mz_t[i] / m;*/

                this->bodies[b].nodes.diff_x_t[i] = this->dt * this->bodies[b].nodes.contact_fx[i] / m;
                this->bodies[b].nodes.diff_y_t[i] = this->dt * this->bodies[b].nodes.contact_fy[i] / m;
                this->bodies[b].nodes.diff_z_t[i] = this->dt * this->bodies[b].nodes.contact_fz[i] / m;
            } else {
                this->bodies[b].nodes.ux[i] = 0;
                this->bodies[b].nodes.uy[i] = 0;
                this->bodies[b].nodes.uz[i] = 0;

                this->bodies[b].nodes.diff_x_t[i] = 0;
                this->bodies[b].nodes.diff_y_t[i] = 0;
                this->bodies[b].nodes.diff_z_t[i] = 0;
            }
        }

        //map back to particles using S transpose
        this->bodies[b].particles.x += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.ux;
        this->bodies[b].particles.y += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.uy;
        this->bodies[b].particles.z += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.uz;

        this->bodies[b].particles.ux += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.ux;
        this->bodies[b].particles.uy += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.uy;
        this->bodies[b].particles.uz += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.uz;

        this->bodies[b].particles.x_t += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.diff_x_t;
        this->bodies[b].particles.y_t += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.diff_y_t;
        this->bodies[b].particles.z_t += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.diff_z_t;

    }
    return;
    //} else {
    //    return this->moveParticlesExplicit2D();
    //}
}

void job_t::moveParticlesImplicit2D(){
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
            double m = this->bodies[b].nodes.m[i];
            if (m != 0) {
                this->bodies[b].nodes.ux[i] = this->dt * (this->bodies[b].nodes.mx_t_k[i] + this->bodies[b].nodes.contact_mx_t[i]) / (2*m);
                this->bodies[b].nodes.uy[i] = this->dt * (this->bodies[b].nodes.my_t_k[i] + this->bodies[b].nodes.contact_my_t[i]) / (2*m);
                this->bodies[b].nodes.uz[i] = 0;
                /*this->bodies[b].nodes.ux[i] = this->dt * this->bodies[b].nodes.contact_mx_t[i] / m;
                this->bodies[b].nodes.uy[i] = this->dt * this->bodies[b].nodes.contact_my_t[i] / m;
                this->bodies[b].nodes.uz[i] = 0;*/

                this->bodies[b].nodes.diff_x_t[i] = this->dt * this->bodies[b].nodes.contact_fx[i] / m;
                this->bodies[b].nodes.diff_y_t[i] = this->dt * this->bodies[b].nodes.contact_fy[i] / m;
                this->bodies[b].nodes.diff_z_t[i] = 0;
            } else {
                this->bodies[b].nodes.ux[i] = 0;
                this->bodies[b].nodes.uy[i] = 0;
                this->bodies[b].nodes.uz[i] = 0;

                this->bodies[b].nodes.diff_x_t[i] = 0;
                this->bodies[b].nodes.diff_y_t[i] = 0;
                this->bodies[b].nodes.diff_z_t[i] = 0;
            }
        }

        //map back to particles using S transpose
        this->bodies[b].particles.x += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.ux;
        this->bodies[b].particles.y += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.uy;

        this->bodies[b].particles.ux += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.ux;
        this->bodies[b].particles.uy += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.uy;

        this->bodies[b].particles.x_t += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.diff_x_t;
        this->bodies[b].particles.y_t += this->bodies[b].Phi.transpose() * this->bodies[b].nodes.diff_y_t;

    }
    return;
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
            double m = this->bodies[b].nodes.m[i];
            if (m!=0){
                this->bodies[b].nodes.ux[i] = this->dt * this->bodies[b].nodes.contact_mx_t[i] / m;
                this->bodies[b].nodes.uy[i] = this->dt * this->bodies[b].nodes.contact_my_t[i] / m;
                //this->bodies[b].nodes.uz[i] = this->dt * this->bodies[b].nodes.contact_mz_t[i] / m;

                this->bodies[b].nodes.diff_x_t[i] = this->dt * this->bodies[b].nodes.contact_fx[i] / m;
                this->bodies[b].nodes.diff_y_t[i] = this->dt * this->bodies[b].nodes.contact_fy[i] / m;
                //this->bodies[b].nodes.diff_z_t[i] = this->dt * this->bodies[b].nodes.contact_fz[i] / m;
            } else {
                this->bodies[b].nodes.ux[i] = 0;
                this->bodies[b].nodes.uy[i] = 0;
                this->bodies[b].nodes.uz[i] = 0;

                this->bodies[b].nodes.diff_x_t[i] = 0;
                this->bodies[b].nodes.diff_y_t[i] = 0;
                this->bodies[b].nodes.diff_z_t[i] = 0;
            }
        }

        //map back to particles using S transpose
        this->bodies[b].particles.x += this->bodies[b].Phi.transpose()*this->bodies[b].nodes.ux;
        this->bodies[b].particles.y += this->bodies[b].Phi.transpose()*this->bodies[b].nodes.uy;
        //this->bodies[b].particles.z += this->bodies[b].Phi.transpose()*this->bodies[b].nodes.uz;

        this->bodies[b].particles.ux += this->bodies[b].Phi.transpose()*this->bodies[b].nodes.ux;
        this->bodies[b].particles.uy += this->bodies[b].Phi.transpose()*this->bodies[b].nodes.uy;
        //this->bodies[b].particles.uz += this->bodies[b].Phi.transpose()*this->bodies[b].nodes.uz;

        this->bodies[b].particles.x_t += this->bodies[b].Phi.transpose()*this->bodies[b].nodes.diff_x_t;
        this->bodies[b].particles.y_t += this->bodies[b].Phi.transpose()*this->bodies[b].nodes.diff_y_t;
        //this->bodies[b].particles.z_t += this->bodies[b].Phi.transpose()*this->bodies[b].nodes.diff_z_t;

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
        //unsafe with current periodic BC
        /*this->bodies[b].nodes.contact_x_t =
                this->bodies[b].nodes.contact_mx_t.array() / this->bodies[b].nodes.m.array();
        this->bodies[b].nodes.contact_y_t =
                this->bodies[b].nodes.contact_my_t.array() / this->bodies[b].nodes.m.array();
        this->bodies[b].nodes.contact_z_t =
                this->bodies[b].nodes.contact_mz_t.array() / this->bodies[b].nodes.m.array();*/
        this->bodies[b].nodes.contact_x_t.setZero();
        this->bodies[b].nodes.contact_y_t.setZero();
        this->bodies[b].nodes.contact_z_t.setZero();
        for (size_t i=0;i<this->num_nodes;i++){
            if (this->bodies[b].nodes.m[i] != 0){
                this->bodies[b].nodes.contact_x_t[i] =
                        this->bodies[b].nodes.contact_mx_t[i] / this->bodies[b].nodes.m[i];
                this->bodies[b].nodes.contact_y_t[i] =
                        this->bodies[b].nodes.contact_my_t[i] / this->bodies[b].nodes.m[i];
                this->bodies[b].nodes.contact_z_t[i] =
                        this->bodies[b].nodes.contact_mz_t[i] / this->bodies[b].nodes.m[i];
            }
        }

        //use to create dummy pvec
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;
        Eigen::VectorXd pvec(numRowsP);

        //calculate particle[i].L[9]
        pvec = this->bodies[b].gradPhiX.transpose() * this->bodies[b].nodes.contact_x_t;
        this->bodies[b].particles.L.col(XX) << pvec;

        pvec = this->bodies[b].gradPhiY.transpose() * this->bodies[b].nodes.contact_x_t;
        this->bodies[b].particles.L.col(XY) << pvec;

        pvec = this->bodies[b].gradPhiZ.transpose() * this->bodies[b].nodes.contact_x_t;
        this->bodies[b].particles.L.col(XZ) << pvec;

        pvec = this->bodies[b].gradPhiX.transpose() * this->bodies[b].nodes.contact_y_t;
        this->bodies[b].particles.L.col(YX) << pvec;

        pvec = this->bodies[b].gradPhiY.transpose() * this->bodies[b].nodes.contact_y_t;
        this->bodies[b].particles.L.col(YY) << pvec;

        pvec = this->bodies[b].gradPhiZ.transpose() * this->bodies[b].nodes.contact_y_t;
        this->bodies[b].particles.L.col(YZ) << pvec;

        pvec = this->bodies[b].gradPhiX.transpose() * this->bodies[b].nodes.contact_z_t;
        this->bodies[b].particles.L.col(ZX) << pvec;

        pvec = this->bodies[b].gradPhiY.transpose() * this->bodies[b].nodes.contact_z_t;
        this->bodies[b].particles.L.col(ZY) << pvec;

        pvec = this->bodies[b].gradPhiZ.transpose() * this->bodies[b].nodes.contact_z_t;
        this->bodies[b].particles.L.col(ZZ) << pvec;
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
        //this->bodies[b].nodes.contact_x_t = this->bodies[b].nodes.contact_mx_t.array()/this->bodies[b].nodes.m.array();
        //this->bodies[b].nodes.contact_y_t = this->bodies[b].nodes.contact_my_t.array()/this->bodies[b].nodes.m.array();
        //this->bodies[b].nodes.contact_z_t = this->bodies[b].nodes.contact_mz_t.array()/this->bodies[b].nodes.m.array();
        this->bodies[b].nodes.contact_x_t.setZero();
        this->bodies[b].nodes.contact_y_t.setZero();
        //this->bodies[b].nodes.contact_z_t.setZero();
        for (size_t i=0;i<this->num_nodes;i++){
            if (this->bodies[b].nodes.m[i] != 0){
                this->bodies[b].nodes.contact_x_t[i] =
                        this->bodies[b].nodes.contact_mx_t[i] / this->bodies[b].nodes.m[i];
                this->bodies[b].nodes.contact_y_t[i] =
                        this->bodies[b].nodes.contact_my_t[i] / this->bodies[b].nodes.m[i];
                //this->bodies[b].nodes.contact_z_t[i] =
                //        this->bodies[b].nodes.contact_mz_t[i] / this->bodies[b].nodes.m[i];
            }
        }

        //use to create dummy pvec
        size_t numRowsP = this->bodies[b].p;
        size_t numColsP = 1;
        //Eigen::MatrixXd pvec(numRowsP,numColsP);
        Eigen::VectorXd pvec(numRowsP);

        //calculate particle[i].L[9]
        pvec = this->bodies[b].gradPhiX.transpose() * this->bodies[b].nodes.contact_x_t;
        this->bodies[b].particles.L.col(XX) << pvec;

        pvec = this->bodies[b].gradPhiY.transpose() * this->bodies[b].nodes.contact_x_t;
        this->bodies[b].particles.L.col(XY) << pvec;

        //pvec = this->bodies[b].gradPhiZ.transpose() * this->bodies[b].nodes.contact_x_t;
        this->bodies[b].particles.L.col(XZ).setZero();

        pvec = this->bodies[b].gradPhiX.transpose() * this->bodies[b].nodes.contact_y_t;
        this->bodies[b].particles.L.col(YX) << pvec;

        pvec = this->bodies[b].gradPhiY.transpose() * this->bodies[b].nodes.contact_y_t;
        this->bodies[b].particles.L.col(YY) << pvec;

        //pvec = this->bodies[b].gradPhiZ.transpose() * this->bodies[b].nodes.contact_y_t;
        this->bodies[b].particles.L.col(YZ).setZero();

        //pvec = this->bodies[b].gradPhiX.transpose() * this->bodies[b].nodes.contact_z_t;
        this->bodies[b].particles.L.col(ZX).setZero();

        //pvec = this->bodies[b].gradPhiY.transpose() * this->bodies[b].nodes.contact_z_t;
        this->bodies[b].particles.L.col(ZY).setZero();

        //pvec = this->bodies[b].gradPhiZ.transpose() * this->bodies[b].nodes.contact_z_t;
        this->bodies[b].particles.L.col(ZZ).setZero();

        //a poor smoothing technique which smooths strain rate instead of strain
        /*if (this->use_smoothing == 1){
            //smoothing [see Mast et. al. 2012 for kinematic locking solution]
            Eigen::VectorXd alpha(this->num_nodes);
            Eigen::VectorXd beta(this->num_nodes);
            for (size_t b = 0; b < this->num_bodies; b++) {
                Eigen::VectorXd pvec(this->bodies[b].p);
                pvec = (this->bodies[b].particles.L.col(XX) + this->bodies[b].particles.L.col(YY)).array();
                alpha = (this->bodies[b].Phi * pvec).array();

                pvec = Eigen::VectorXd::Ones(this->bodies[b].p);
                beta = this->bodies[b].Phi * pvec;

                alpha = alpha.array()/beta.array();

                for (size_t i=0;i<this->num_nodes;i++){
                    if (this->bodies[b].nodes.m[i] == 0){
                        alpha[i] = 0;
                    }
                }

                //trE and trT
                Eigen::VectorXd trL = this->bodies[b].Phi.transpose() * alpha;
                pvec = (trL - this->bodies[b].particles.L.col(XX) - this->bodies[b].particles.L.col(YY))/2.0;

                this->bodies[b].particles.L.col(XX) += pvec;
                this->bodies[b].particles.L.col(YY) += pvec;
            }
        }*/
    }

    return;
}

void job_t::updateDensity(){
    //update density of particles per sachiths code
    for (size_t b=0;b<this->num_bodies;b++){
        for (size_t i=0;i<this->bodies[b].p;i++){
            double trL = this->bodies[b].particles.L(i,XX)+this->bodies[b].particles.L(i,YY)+this->bodies[b].particles.L(i,ZZ);
            this->bodies[b].particles.v[i] *= std::exp(this->dt * trL);
        }
    }
    return;
}

void job_t::updateTrialDensity(){
    //update density of particles per sachiths code
    for (size_t b=0;b<this->num_bodies;b++){
        for (size_t i=0;i<this->bodies[b].p;i++){
            double trL = this->bodies[b].particles.L(i,XX)+this->bodies[b].particles.L(i,YY)+this->bodies[b].particles.L(i,ZZ);
            this->bodies[b].particles.v_trial[i] = this->bodies[b].particles.v[i] * std::exp(this->dt * trL);
        }
    }
    return;
}

void job_t::updateStress(){
    //calculate stress
    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].material.calculate_stress(&(this->bodies[b]),this->dt);
    }

    //return;

    if (this->use_smoothing == 1) {
        //smoothing [see Mast et. al. 2012 for kinematic locking solution]
        Eigen::VectorXd alpha(this->num_nodes);
        Eigen::VectorXd beta(this->num_nodes);
        for (size_t b = 0; b < this->num_bodies; b++) {
            Eigen::VectorXd pvec(this->bodies[b].p);
            //pvec = (this->bodies[b].particles.v - this->bodies[b].particles.v0).array() / this->bodies[b].particles.v0.array() *
            pvec = (this->bodies[b].particles.v.array() / this->bodies[b].particles.v0.array()).array().log() * this->bodies[b].particles.m.array();
            alpha = (this->bodies[b].Phi * pvec).array() / this->bodies[b].nodes.m.array();

            pvec = (this->bodies[b].particles.T.col(XX) + this->bodies[b].particles.T.col(YY) +
                    this->bodies[b].particles.T.col(ZZ)).array() * this->bodies[b].particles.m.array();
            beta = (this->bodies[b].Phi * pvec).array() / this->bodies[b].nodes.m.array();

            for (size_t i=0;i<this->num_nodes;i++){
                if (this->bodies[b].nodes.m[i] == 0){
                    alpha[i] = 0;
                    beta[i] = 0;
                }
            }

            //trE and trT
            /*double mu = 1;
            pvec = (this->bodies[b].Phi.transpose() * alpha);
            Eigen::VectorXd trE = pvec.array() * mu + (this->bodies[b].particles.v.array() / this->bodies[b].particles.v0.array()).array().log() * (1-mu);
            pvec = (this->bodies[b].Phi.transpose() * beta);
            Eigen::VectorXd trT = pvec.array() * mu + (this->bodies[b].particles.T.col(XX) + this->bodies[b].particles.T.col(YY) + this->bodies[b].particles.T.col(ZZ)).array() * (1-mu);
            this->bodies[b].material.volumetric_smoothing(&(this->bodies[b]), trE, trT);
            */
            Eigen::VectorXd trE = (this->bodies[b].Phi.transpose() * alpha);
            Eigen::VectorXd trT = (this->bodies[b].Phi.transpose() * beta);
            this->bodies[b].material.volumetric_smoothing(&(this->bodies[b]), trE, trT);
            //adjust particle volume and density
            //this->bodies[b].particles.v = this->bodies[b].particles.v0.array() * trE.array().exp();
            //        this->bodies[b].particles.v0.array() + this->bodies[b].particles.v0.array() * trE.array();
        }
        /*Eigen::VectorXd alpha(this->num_elements);
        Eigen::VectorXd beta(this->num_elements);
        Eigen::VectorXd m(this->num_elements);
        for (size_t b=0;b<this->num_bodies;b++) {
            alpha.setZero();
            beta.setZero();
            m.setZero();
            this->bodies[b].particles.updateElementIDs(this);
            Eigen::VectorXd pvec(this->bodies[b].p);

            //trE
            pvec = (this->bodies[b].particles.v - this->bodies[b].particles.v0).array() / this->bodies[b].particles.v0.array() *
                   this->bodies[b].particles.m.array();
            for (size_t i = 0; i < this->bodies[b].p; i++) {
                size_t e = this->bodies[b].particles.elementIDs[i];
                if (this->bodies[b].particles.active[i]) {
                    alpha[e] += pvec[i];
                }
            }

            //trA
            pvec = (this->bodies[b].particles.T.col(XX) + this->bodies[b].particles.T.col(YY) +
                    this->bodies[b].particles.T.col(ZZ)).array() * this->bodies[b].particles.m.array();
            for (size_t i = 0; i < this->bodies[b].p; i++) {
                size_t e = this->bodies[b].particles.elementIDs[i];
                if (this->bodies[b].particles.active[i]) {
                    beta[e] += pvec[i];
                }
            }

            //m
            pvec = this->bodies[b].particles.m.array();
            for (size_t i = 0; i < this->bodies[b].p; i++) {
                size_t e = this->bodies[b].particles.elementIDs[i];
                if (this->bodies[b].particles.active[i]) {
                    m[e] += pvec[i];
                }
            }

            //trE and trT
            Eigen::VectorXd trE(this->bodies[b].p);
            Eigen::VectorXd trT(this->bodies[b].p);
            trE.setZero();
            trT.setZero();
            for (size_t i = 0; i < this->bodies[b].p; i++) {
                size_t e = this->bodies[b].particles.elementIDs[i];
                if (this->bodies[b].particles.active[i]) {
                    trE[i] = alpha[e] / m[e];
                    trT[i] = beta[e] / m[e];
                }
            }

            this->bodies[b].material.volumetric_smoothing(&(this->bodies[b]), trE, trT);

            //adjust particle volume and density
            this->bodies[b].particles.v =
                    this->bodies[b].particles.v0.array() + this->bodies[b].particles.v0.array() * trE.array();

        }*/
    }

    return;
}

void job_t::updateTrialStress(){
    //calculate stress
    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].material.calculate_stress_implicit(&(this->bodies[b]),this->dt);
    }

    if (this->use_smoothing == 1) {
        //smoothing [see Mast et. al. 2012 for kinematic locking solution]
        Eigen::VectorXd alpha(this->num_nodes);
        Eigen::VectorXd beta(this->num_nodes);
        for (size_t b = 0; b < this->num_bodies; b++) {
            Eigen::VectorXd pvec(this->bodies[b].p);
            pvec = (this->bodies[b].particles.v_trial.array() / this->bodies[b].particles.v0.array()).array().log() * this->bodies[b].particles.m.array();
            alpha = (this->bodies[b].Phi * pvec).array() / this->bodies[b].nodes.m.array();

            pvec = (this->bodies[b].particles.Ttrial.col(XX) + this->bodies[b].particles.Ttrial.col(YY) +
                    this->bodies[b].particles.Ttrial.col(ZZ)).array() * this->bodies[b].particles.m.array();
            beta = (this->bodies[b].Phi * pvec).array() / this->bodies[b].nodes.m.array();


            for (size_t i=0;i<this->num_nodes;i++){
                if (this->bodies[b].nodes.m[i] == 0){
                    alpha[i] = 0;
                    beta[i] = 0;
                }
            }

            //trE and trT
            Eigen::VectorXd trE = this->bodies[b].Phi.transpose() * alpha;
            Eigen::VectorXd trT = this->bodies[b].Phi.transpose() * beta;
            this->bodies[b].material.volumetric_smoothing_implicit(&(this->bodies[b]), trE, trT);
        }
        /*Eigen::VectorXd alpha(this->num_elements);
        Eigen::VectorXd beta(this->num_elements);
        Eigen::VectorXd m(this->num_elements);
        for (size_t b=0;b<this->num_bodies;b++) {
            alpha.setZero();
            beta.setZero();
            m.setZero();
            this->bodies[b].particles.updateElementIDs(this);
            Eigen::VectorXd pvec(this->bodies[b].p);

            //trE
            pvec = (this->bodies[b].particles.v_trial - this->bodies[b].particles.v0).array() / this->bodies[b].particles.v0.array() *
                   this->bodies[b].particles.m.array();
            for (size_t i = 0; i < this->bodies[b].p; i++) {
                size_t e = this->bodies[b].particles.elementIDs[i];
                if (this->bodies[b].particles.active[i]) {
                    alpha[e] += pvec[i];
                }
            }

            //trA
            pvec = (this->bodies[b].particles.Ttrial.col(XX) + this->bodies[b].particles.Ttrial.col(YY) +
                    this->bodies[b].particles.Ttrial.col(ZZ)).array() * this->bodies[b].particles.m.array();
            for (size_t i = 0; i < this->bodies[b].p; i++) {
                size_t e = this->bodies[b].particles.elementIDs[i];
                if (this->bodies[b].particles.active[i]) {
                    beta[e] += pvec[i];
                }
            }

            //m
            pvec = this->bodies[b].particles.m.array();
            for (size_t i = 0; i < this->bodies[b].p; i++) {
                size_t e = this->bodies[b].particles.elementIDs[i];
                if (this->bodies[b].particles.active[i]) {
                    m[e] += pvec[i];
                }
            }

            //trE and trT
            Eigen::VectorXd trE(this->bodies[b].p);
            Eigen::VectorXd trT(this->bodies[b].p);
            trE.setZero();
            trT.setZero();
            for (size_t i = 0; i < this->bodies[b].p; i++) {
                size_t e = this->bodies[b].particles.elementIDs[i];
                if (this->bodies[b].particles.active[i]) {
                    trE[i] = alpha[e] / m[e];
                    trT[i] = beta[e] / m[e];
                }
            }

            this->bodies[b].material.volumetric_smoothing_implicit(&(this->bodies[b]), trE, trT);

            //adjust particle volume and density
            this->bodies[b].particles.v_trial =
                    this->bodies[b].particles.v0.array() + this->bodies[b].particles.v0.array() * trE.array();

        }*/
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
        p_m_bx = this->bodies[b].particles.m.array() * this->bodies[b].particles.bx.array();
        p_m_by = this->bodies[b].particles.m.array() * this->bodies[b].particles.by.array();
        p_m_bz = this->bodies[b].particles.m.array() * this->bodies[b].particles.bz.array();

        this->bodies[b].nodes.fx = this->bodies[b].Phi * p_m_bx; //need to add stress
        this->bodies[b].nodes.fy = this->bodies[b].Phi * p_m_by; //need to add stress
        this->bodies[b].nodes.fz = this->bodies[b].Phi * p_m_bz; //need to add stress

        //use to create dummy pvec and ones
        Eigen::VectorXd pvec(numRowsP);

        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.Ttrial.col(XX).array();
        this->bodies[b].nodes.fx -= this->bodies[b].gradPhiX*pvec;

        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.Ttrial.col(XY).array();
        this->bodies[b].nodes.fx -= this->bodies[b].gradPhiY*pvec;
        this->bodies[b].nodes.fy -= this->bodies[b].gradPhiX*pvec;

        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.Ttrial.col(XZ).array();
        this->bodies[b].nodes.fx -= this->bodies[b].gradPhiZ*pvec;
        this->bodies[b].nodes.fz -= this->bodies[b].gradPhiX*pvec;

        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.Ttrial.col(YY).array();
        this->bodies[b].nodes.fy -= this->bodies[b].gradPhiY*pvec;

        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.Ttrial.col(YZ).array();
        this->bodies[b].nodes.fy -= this->bodies[b].gradPhiZ*pvec;
        this->bodies[b].nodes.fz -= this->bodies[b].gradPhiY*pvec;

        pvec = this->bodies[b].particles.v.array() * this->bodies[b].particles.Ttrial.col(ZZ).array();
        this->bodies[b].nodes.fz -= this->bodies[b].gradPhiZ*pvec;
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

        this->bodies[b].nodes.Rx = this->bodies[b].nodes.x_t_trial.array()*this->bodies[b].nodes.m.array();
        this->bodies[b].nodes.Rx -= this->dt*this->bodies[b].nodes.fx_L;
        this->bodies[b].nodes.Rx -= this->bodies[b].nodes.mx_t_k;

        this->bodies[b].nodes.Ry = this->bodies[b].nodes.y_t_trial.array()*this->bodies[b].nodes.m.array();
        this->bodies[b].nodes.Ry -= this->dt*this->bodies[b].nodes.fy_L;
        this->bodies[b].nodes.Ry -= this->bodies[b].nodes.my_t_k;

        this->bodies[b].nodes.Rz = this->bodies[b].nodes.z_t_trial.array()*this->bodies[b].nodes.m.array();
        this->bodies[b].nodes.Rz -= this->dt*this->bodies[b].nodes.fz_L;
        this->bodies[b].nodes.Rz -= this->bodies[b].nodes.mz_t_k;

        for (size_t i=0;i<this->num_nodes;i++){
            if (this->bodies[b].nodes.m[i] == 0){//<= TOL){
                //nodal boundary
                this->bodies[b].nodes.Rx[i] = 0;
                this->bodies[b].nodes.Ry[i] = 0;
                this->bodies[b].nodes.Rz[i] = 0;
            /*} else if (this->bodies[0].nodes.m[i] > TOL && this->bodies[1].nodes.m[i] > TOL){
                //nodal contact
                this->bodies[b].nodes.Rx[i] = 0;
                this->bodies[b].nodes.Ry[i] = 0;
                this->bodies[b].nodes.Rz[i] = 0;*/
            } else {
                if (this->u_dirichlet_mask[NODAL_DOF*i + XDOF_IDX] != 0){
                    this->bodies[b].nodes.Rx[i] = 0;
                } else {
                    this->bodies[b].nodes.Rx[i] /= this->bodies[b].nodes.m[i];
                }

                if (this->u_dirichlet_mask[NODAL_DOF*i + YDOF_IDX] != 0){
                    this->bodies[b].nodes.Ry[i] = 0;
                } else {
                    this->bodies[b].nodes.Ry[i] /= this->bodies[b].nodes.m[i];
                }

                if (this->u_dirichlet_mask[NODAL_DOF*i + ZDOF_IDX] != 0){
                    this->bodies[b].nodes.Rz[i] = 0;
                } else {
                    this->bodies[b].nodes.Rz[i] /= this->bodies[b].nodes.m[i];
                }
                //this->bodies[b].Ry[i] /= this->bodies[b].nodes.m[i];
                //this->bodies[b].Rz[i] /= this->bodies[b].nodes.m[i];
            }

            /*if (!std::isfinite(this->bodies[b].nodes.Rx[i])){
                std::cout << i << ", " << this->bodies[b].nodes.m[i] << ", " << this->bodies[b].nodes.fx_L[i] << ", " << this->bodies[b].nodes.mx_t_k[i] << ", " << this->bodies[b].nodes.x_t_trial[i] << std::endl;
            }*/

        }
    }
}

void job_t::moveGridImplicitCG() {
    //full conjugate gradient mathod for solving Js = -F(v)

    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].nodes.mx_t_k = this->bodies[b].nodes.contact_mx_t;
        this->bodies[b].nodes.my_t_k = this->bodies[b].nodes.contact_my_t;
        this->bodies[b].nodes.mz_t_k = this->bodies[b].nodes.contact_mz_t;
    }

    //move grid
    //this->moveGridExplicit();

    //save forces and initial trial velocity
    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].nodes.x_t_n.setZero();
        this->bodies[b].nodes.y_t_n.setZero();
        this->bodies[b].nodes.z_t_n.setZero();

        for (size_t i=0;i<this->num_nodes;i++){
            if (this->bodies[b].nodes.m[i] > TOL) {
                this->bodies[b].nodes.x_t_n[i] = this->bodies[b].nodes.contact_mx_t[i] / this->bodies[b].nodes.m[i];
                this->bodies[b].nodes.y_t_n[i] = this->bodies[b].nodes.contact_my_t[i] / this->bodies[b].nodes.m[i];
                this->bodies[b].nodes.z_t_n[i] = this->bodies[b].nodes.contact_mz_t[i] / this->bodies[b].nodes.m[i];
            }
        }

        this->bodies[b].nodes.x_t_trial = this->bodies[b].nodes.x_t_n;
        this->bodies[b].nodes.y_t_trial = this->bodies[b].nodes.y_t_n;
        this->bodies[b].nodes.z_t_trial = this->bodies[b].nodes.z_t_n;
    }

    //calculate L on particles from initial contact velocity
    this->calculateStrainRate();

    //update particle densities
    this->updateTrialDensity();

    //add body forces
    time_varying_loads(this);

    //material stress update
    this->updateTrialStress();

    //map particle stress back to nodes
    this->mapTrialStress2Grid();

    //add contact forces
    this->addContactForces();

    //enforce boundary conditions
    this->addBoundaryConditions();

    //save explicit solution
    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].nodes.fx_L = this->bodies[b].nodes.contact_fx;
        this->bodies[b].nodes.fy_L = this->bodies[b].nodes.contact_fy;
        this->bodies[b].nodes.fz_L = this->bodies[b].nodes.contact_fz;
    }

    //calculate residual
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
            this->bodies[b].nodes.Rvx = this->bodies[b].nodes.Rx;
            this->bodies[b].nodes.Rvy = this->bodies[b].nodes.Ry;
            this->bodies[b].nodes.Rvz = this->bodies[b].nodes.Rz;

            this->bodies[b].nodes.rk.setZero();

            this->bodies[b].nodes.rk << -this->bodies[b].nodes.Rvx, -this->bodies[b].nodes.Rvy, -this->bodies[b].nodes.Rvz;

            this->bodies[b].nodes.rkMin = this->bodies[b].nodes.rk;
            this->bodies[b].nodes.skMin = this->bodies[b].nodes.sk;

            this->bodies[b].nodes.rhok = this->bodies[b].nodes.rk.squaredNorm();

            rhoSum += this->bodies[b].nodes.rhok;

            this->bodies[b].nodes.pk = this->bodies[b].nodes.rk / std::sqrt(this->bodies[b].nodes.rhok);

            this->bodies[b].nodes.sk.setZero();
        }

        std::cout << "\nRHO: " << rhoSum << " ?< " << rhoTOL << std::endl;

        //solve for 's' to iterate 'v' [Sulsky 2003]
        size_t k = 0;
        while (k < 100 && rhoSum > rhoTOL && rhoSum < R_MAX) { //(this->num_bodies * this->num_nodes * R_TOL)) {
            double h = this->linearStepSize;//-6; //per paper suggestion

            for (size_t b = 0; b < this->num_bodies; b++) {
                double vNorm = std::sqrt(
                        this->bodies[b].nodes.x_t_n.squaredNorm() + this->bodies[b].nodes.y_t_n.squaredNorm()
                        + this->bodies[b].nodes.z_t_n.squaredNorm());
                if (vNorm <= 0) {
                    vNorm = 1.0;
                }
                double sNorm = this->bodies[b].nodes.pk.norm();
                if (sNorm>0) {
                    this->bodies[b].nodes.x_t_trial =
                            this->bodies[b].nodes.x_t_n + h * vNorm * this->bodies[b].nodes.pk.segment(0,this->bodies[b].n) / sNorm;
                    this->bodies[b].nodes.y_t_trial =
                            this->bodies[b].nodes.y_t_n + h * vNorm * this->bodies[b].nodes.pk.segment(this->bodies[b].n,this->bodies[b].n)/ sNorm;
                    this->bodies[b].nodes.z_t_trial =
                            this->bodies[b].nodes.z_t_n + h * vNorm * this->bodies[b].nodes.pk.segment(2*this->bodies[b].n,this->bodies[b].n) / sNorm;
                } else {
                    this->bodies[b].nodes.x_t_trial = this->bodies[b].nodes.x_t_n;
                    this->bodies[b].nodes.y_t_trial = this->bodies[b].nodes.y_t_n;
                    this->bodies[b].nodes.z_t_trial = this->bodies[b].nodes.z_t_n;
                }

                this->bodies[b].nodes.mx_t = this->bodies[b].nodes.x_t_trial.array() * this->bodies[b].nodes.m.array();
                this->bodies[b].nodes.my_t = this->bodies[b].nodes.y_t_trial.array() * this->bodies[b].nodes.m.array();
                this->bodies[b].nodes.mz_t = this->bodies[b].nodes.z_t_trial.array() * this->bodies[b].nodes.m.array();

                this->bodies[b].nodes.x_t = this->bodies[b].nodes.x_t_trial;
                this->bodies[b].nodes.y_t = this->bodies[b].nodes.y_t_trial;
                this->bodies[b].nodes.z_t = this->bodies[b].nodes.z_t_trial;
            }

            //add contact forces to trial velocity
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
                this->bodies[b].nodes.fx_L = this->bodies[b].nodes.contact_fx;
                this->bodies[b].nodes.fy_L = this->bodies[b].nodes.contact_fy;
                this->bodies[b].nodes.fz_L = this->bodies[b].nodes.contact_fz;
            }

            this->calculateImplicitResidual();

            for (size_t b = 0; b < this->num_bodies; b++) {
                double vNorm = std::sqrt(
                        this->bodies[b].nodes.x_t_n.squaredNorm() + this->bodies[b].nodes.y_t_n.squaredNorm() +
                        this->bodies[b].nodes.z_t_n.squaredNorm());
                if (vNorm <= 0) {
                    vNorm = 1.0;
                }
                double sNorm = this->bodies[b].nodes.pk.norm();
                this->bodies[b].nodes.DhRx = sNorm / (h * vNorm) * (this->bodies[b].nodes.Rx - this->bodies[b].nodes.Rvx);
                this->bodies[b].nodes.DhRy = sNorm / (h * vNorm) * (this->bodies[b].nodes.Ry - this->bodies[b].nodes.Rvy);
                this->bodies[b].nodes.DhRz = sNorm / (h * vNorm) * (this->bodies[b].nodes.Rz - this->bodies[b].nodes.Rvz);

                this->bodies[b].nodes.wk << this->bodies[b].nodes.DhRx, this->bodies[b].nodes.DhRy, this->bodies[b].nodes.DhRz;
                this->bodies[b].nodes.ak = this->bodies[b].nodes.rhok / (this->bodies[b].nodes.pk.transpose() * this->bodies[b].nodes.wk);

                this->bodies[b].nodes.sk += this->bodies[b].nodes.ak * this->bodies[b].nodes.pk;

                this->bodies[b].nodes.rk -= this->bodies[b].nodes.ak * this->bodies[b].nodes.wk;

                this->bodies[b].nodes.bk = this->bodies[b].nodes.rk.squaredNorm() / this->bodies[b].nodes.rhok;

                this->bodies[b].nodes.rhok = this->bodies[b].nodes.rk.squaredNorm();

                this->bodies[b].nodes.pk = this->bodies[b].nodes.rk + this->bodies[b].nodes.bk * this->bodies[b].nodes.pk;

                //store if minimum
                if (this->bodies[b].nodes.rk.norm() < this->bodies[b].nodes.rkMin.norm() || k == 0){
                    this->bodies[b].nodes.rkMin = this->bodies[b].nodes.rk;
                    this->bodies[b].nodes.skMin = this->bodies[b].nodes.sk;
                }
            }
            k += 1;
            rhoSum = 0;
            for (size_t b = 0; b < this->num_bodies; b++) {
                rhoSum += this->bodies[b].nodes.rhok;
            }
            std::cout << "\rn: " << nIter << " k: " << k <<  " r: " << rhoSum << "      \r" << std::flush;
        }

        for (size_t b = 0; b < this->num_bodies; b++) {
            /*this->bodies[b].nodes.x_t_n += this->bodies[b].nodes.sk.segment(0,this->bodies[b].n);
            this->bodies[b].nodes.y_t_n += this->bodies[b].nodes.sk.segment(this->bodies[b].n,this->bodies[b].n);
            this->bodies[b].nodes.z_t_n += this->bodies[b].nodes.sk.segment(2*this->bodies[b].n,this->bodies[b].n);*/
            this->bodies[b].nodes.x_t_n += this->bodies[b].nodes.skMin.segment(0,this->bodies[b].n);
            this->bodies[b].nodes.y_t_n += this->bodies[b].nodes.skMin.segment(this->bodies[b].n,this->bodies[b].n);
            this->bodies[b].nodes.z_t_n += this->bodies[b].nodes.skMin.segment(2*this->bodies[b].n,this->bodies[b].n);

            this->bodies[b].nodes.x_t_trial = this->bodies[b].nodes.x_t_n;
            this->bodies[b].nodes.y_t_trial = this->bodies[b].nodes.y_t_n;
            this->bodies[b].nodes.z_t_trial = this->bodies[b].nodes.z_t_n;

            this->bodies[b].nodes.mx_t = this->bodies[b].nodes.x_t_trial.array() * this->bodies[b].nodes.m.array();
            this->bodies[b].nodes.my_t = this->bodies[b].nodes.y_t_trial.array() * this->bodies[b].nodes.m.array();
            this->bodies[b].nodes.mz_t = this->bodies[b].nodes.z_t_trial.array() * this->bodies[b].nodes.m.array();

            this->bodies[b].nodes.x_t = this->bodies[b].nodes.x_t_trial;
            this->bodies[b].nodes.y_t = this->bodies[b].nodes.y_t_trial;
            this->bodies[b].nodes.z_t = this->bodies[b].nodes.z_t_trial;
        }

        //std::cout << "n: " << nIter << " k: " << k <<  " r: " << rhoSum << std::endl;
        nIter += 1;

        //add contact forces to trial velocity
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
            this->bodies[b].nodes.fx_L = this->bodies[b].nodes.contact_fx;
            this->bodies[b].nodes.fy_L = this->bodies[b].nodes.contact_fy;
            this->bodies[b].nodes.fz_L = this->bodies[b].nodes.contact_fz;
        }

        this->calculateImplicitResidual();

        //setup for next iteration
        rhoSum = 0;
        for (size_t b = 0; b < this->num_bodies; b++) {
            this->bodies[b].nodes.rk << -this->bodies[b].nodes.Rx, -this->bodies[b].nodes.Ry, -this->bodies[b].nodes.Rz;

            this->bodies[b].nodes.rhok = this->bodies[b].nodes.rk.squaredNorm();

            rhoSum += this->bodies[b].nodes.rhok;
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
        this->bodies[b].nodes.mx_t_k = this->bodies[b].nodes.contact_mx_t;
        this->bodies[b].nodes.my_t_k = this->bodies[b].nodes.contact_my_t;
        this->bodies[b].nodes.mz_t_k = this->bodies[b].nodes.contact_mz_t;
    }

    //move grid
    //this->moveGridExplicit();

    //save forces and initial trial velocity
    for (size_t b=0;b<this->num_bodies;b++){
        this->bodies[b].nodes.x_t_n.setZero();
        this->bodies[b].nodes.y_t_n.setZero();
        this->bodies[b].nodes.z_t_n.setZero();

        double activeNodes = 0;
        for (size_t i=0;i<this->num_nodes;i++){
            if (this->bodies[b].nodes.m[i] > TOL) {
                activeNodes += 1;
                this->bodies[b].nodes.x_t_n[i] = this->bodies[b].nodes.contact_mx_t[i] / this->bodies[b].nodes.m[i];
                this->bodies[b].nodes.y_t_n[i] = this->bodies[b].nodes.contact_my_t[i] / this->bodies[b].nodes.m[i];
                this->bodies[b].nodes.z_t_n[i] = this->bodies[b].nodes.contact_mz_t[i] / this->bodies[b].nodes.m[i];
            }
        }
        //std::cout << "\nactive nodes: " << activeNodes << std::endl;

        this->bodies[b].nodes.x_t_trial = this->bodies[b].nodes.x_t_n;
        this->bodies[b].nodes.y_t_trial = this->bodies[b].nodes.y_t_n;
        this->bodies[b].nodes.z_t_trial = this->bodies[b].nodes.z_t_n;
    }

    //calculate L on particles from initial contact velocities
    this->calculateStrainRate();

    //update particle densities
    this->updateTrialDensity();

    //add body forces
    time_varying_loads(this);

    //material stress update
    this->updateTrialStress();

    //map particle stress back to nodes
    this->mapTrialStress2Grid();

    //add contact forces
    this->addContactForces();

    //enforce boundary conditions
    this->addBoundaryConditions();

    //save final forces on nodes
    for (size_t b=0;b<this->num_bodies;b++){

        this->bodies[b].nodes.fx_L = this->bodies[b].nodes.contact_fx;
        this->bodies[b].nodes.fy_L = this->bodies[b].nodes.contact_fy;
        this->bodies[b].nodes.fz_L = this->bodies[b].nodes.contact_fz;
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
            this->bodies[b].nodes.Rvx = this->bodies[b].nodes.Rx;
            this->bodies[b].nodes.Rvy = this->bodies[b].nodes.Ry;
            this->bodies[b].nodes.Rvz = this->bodies[b].nodes.Rz;

            //initialize bicgstab variables
            this->bodies[b].nodes.sk.setZero();
            this->bodies[b].nodes.wk.setZero();
            this->bodies[b].nodes.pk.setZero();
            this->bodies[b].nodes.ak = 1;
            this->bodies[b].nodes.rhok = 1;
            this->bodies[b].nodes.ok = 1;
            this->bodies[b].nodes.rk.setZero();

            this->bodies[b].nodes.rk << -this->bodies[b].nodes.Rvx, -this->bodies[b].nodes.Rvy, -this->bodies[b].nodes.Rvz;
            this->bodies[b].nodes.r0 = this->bodies[b].nodes.rk;

            this->bodies[b].nodes.rkMin = this->bodies[b].nodes.rk;
            this->bodies[b].nodes.skMin = this->bodies[b].nodes.sk;

            rhoSum += this->bodies[b].nodes.rk.squaredNorm();
        }

        std::cout << "RHO: " << rhoSum << " ?< " << rhoTOL << std::endl;/*<<  " vnorm: " << this->bodies[0].node_x_t_n.lpNorm<Eigen::Infinity>() <<
        "," << this->bodies[0].node_y_t_n.lpNorm<Eigen::Infinity>() <<
        "," << this->bodies[0].node_z_t_n.lpNorm<Eigen::Infinity>() << std::endl;*/

        //solve for 's' to iterate 'v' [Sulsky 2003]
        size_t k = 0;
        bool isConverging = true;
        while (k < 100 && rhoSum > rhoTOL && rhoSum < R_MAX && isConverging) {
            double h = this->linearStepSize;

            for (size_t b = 0; b < this->num_bodies; b++) {

                this->bodies[b].nodes.bk = (this->bodies[b].nodes.r0.dot(this->bodies[b].nodes.rk))/(this->bodies[b].nodes.rhok)*this->bodies[b].nodes.ak/this->bodies[b].nodes.ok;
                this->bodies[b].nodes.rhok = this->bodies[b].nodes.r0.dot(this->bodies[b].nodes.rk);
                this->bodies[b].nodes.pk = this->bodies[b].nodes.rk + this->bodies[b].nodes.bk*(this->bodies[b].nodes.pk - this->bodies[b].nodes.ok*this->bodies[b].nodes.wk);

                //wk = DhDF(vn,pk)
                double vNorm = std::sqrt(
                        this->bodies[b].nodes.x_t_n.squaredNorm() + this->bodies[b].nodes.y_t_n.squaredNorm()
                        + this->bodies[b].nodes.z_t_n.squaredNorm());
                if (vNorm <= 0) {
                    vNorm = 1.0;
                }
                double sNorm = this->bodies[b].nodes.pk.norm();
                if (sNorm>0) {
                    this->bodies[b].nodes.x_t_trial =
                            this->bodies[b].nodes.x_t_n + h * vNorm * this->bodies[b].nodes.pk.segment(0,this->bodies[b].n) / sNorm;
                    this->bodies[b].nodes.y_t_trial =
                            this->bodies[b].nodes.y_t_n + h * vNorm * this->bodies[b].nodes.pk.segment(this->bodies[b].n,this->bodies[b].n)/ sNorm;
                    this->bodies[b].nodes.z_t_trial =
                            this->bodies[b].nodes.z_t_n + h * vNorm * this->bodies[b].nodes.pk.segment(2*this->bodies[b].n,this->bodies[b].n) / sNorm;
                } else {
                    isConverging = false;
                    /*std::cout << k << " 1 " << sNorm << std::endl;
                    this->bodies[b].nodes.x_t_trial = this->bodies[b].nodes.x_t_n;
                    this->bodies[b].nodes.y_t_trial = this->bodies[b].nodes.y_t_n;
                    this->bodies[b].nodes.z_t_trial = this->bodies[b].nodes.z_t_n;*/
                }

                this->bodies[b].nodes.mx_t = this->bodies[b].nodes.x_t_trial.array() * this->bodies[b].nodes.m.array();
                this->bodies[b].nodes.my_t = this->bodies[b].nodes.y_t_trial.array() * this->bodies[b].nodes.m.array();
                this->bodies[b].nodes.mz_t = this->bodies[b].nodes.z_t_trial.array() * this->bodies[b].nodes.m.array();

                this->bodies[b].nodes.x_t = this->bodies[b].nodes.x_t_trial;
                this->bodies[b].nodes.y_t = this->bodies[b].nodes.y_t_trial;
                this->bodies[b].nodes.z_t = this->bodies[b].nodes.z_t_trial;
            }

            //if failing to converge, break
            if(!isConverging){break;}

            //add contact forces to trial velocity
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
                this->bodies[b].nodes.fx_L = this->bodies[b].nodes.contact_fx;
                this->bodies[b].nodes.fy_L = this->bodies[b].nodes.contact_fy;
                this->bodies[b].nodes.fz_L = this->bodies[b].nodes.contact_fz;
            }

            this->calculateImplicitResidual();

            for (size_t b = 0; b < this->num_bodies; b++) {
                double vNorm = std::sqrt(
                        this->bodies[b].nodes.x_t_n.squaredNorm() + this->bodies[b].nodes.y_t_n.squaredNorm() +
                        this->bodies[b].nodes.z_t_n.squaredNorm());
                if (vNorm <= 0) {
                    vNorm = 1.0;
                }
                double sNorm = this->bodies[b].nodes.pk.norm();
                if (sNorm > 0) {
                    //isConverging = true;
                    this->bodies[b].nodes.DhRx = sNorm / (h * vNorm) * (this->bodies[b].nodes.Rx - this->bodies[b].nodes.Rvx);
                    this->bodies[b].nodes.DhRy = sNorm / (h * vNorm) * (this->bodies[b].nodes.Ry - this->bodies[b].nodes.Rvy);
                    this->bodies[b].nodes.DhRz = sNorm / (h * vNorm) * (this->bodies[b].nodes.Rz - this->bodies[b].nodes.Rvz);

                    this->bodies[b].nodes.wk << this->bodies[b].nodes.DhRx, this->bodies[b].nodes.DhRy, this->bodies[b].nodes.DhRz;
                    this->bodies[b].nodes.ak = this->bodies[b].nodes.rhok / (this->bodies[b].nodes.r0.dot(this->bodies[b].nodes.wk));
                    this->bodies[b].nodes.hk = this->bodies[b].nodes.sk + this->bodies[b].nodes.ak * this->bodies[b].nodes.pk;

                    if (!std::isfinite(this->bodies[b].nodes.ak)) {
                        std::cout << "\nak is infinite\n";
                    }

                    //check hk for convergence?

                    this->bodies[b].nodes.qk = this->bodies[b].nodes.rk - this->bodies[b].nodes.ak * this->bodies[b].nodes.wk;
                } else {
                    isConverging = false;
                    //std::cout << "sNorm <= " << TOL << std::endl;
                    /*std::cout << k << " 2 " << sNorm << std::endl;
                    isConverging = false;
                    this->bodies[b].nodes.hk = this->bodies[b].nodes.sk;
                    this->bodies[b].nodes.qk = this->bodies[b].nodes.rk;*/
                    //let ak retain former value for bk calculation at beginning of step
                }

                //calculate t = DhDF(vn,q)
                sNorm = this->bodies[b].nodes.qk.norm();
                if (sNorm>0) {
                    this->bodies[b].nodes.x_t_trial =
                            this->bodies[b].nodes.x_t_n + h * vNorm * this->bodies[b].nodes.qk.segment(0,this->bodies[b].n) / sNorm;
                    this->bodies[b].nodes.y_t_trial =
                            this->bodies[b].nodes.y_t_n + h * vNorm * this->bodies[b].nodes.qk.segment(this->bodies[b].n,this->bodies[b].n)/ sNorm;
                    this->bodies[b].nodes.z_t_trial =
                            this->bodies[b].nodes.z_t_n + h * vNorm * this->bodies[b].nodes.qk.segment(2*this->bodies[b].n,this->bodies[b].n) / sNorm;
                } else {
                    isConverging = false;
                    /*std::cout << k << " 3 " << sNorm << std::endl;
                    this->bodies[b].nodes.x_t_trial = this->bodies[b].nodes.x_t_n;
                    this->bodies[b].nodes.y_t_trial = this->bodies[b].nodes.y_t_n;
                    this->bodies[b].nodes.z_t_trial = this->bodies[b].nodes.z_t_n;*/
                }

                this->bodies[b].nodes.mx_t = this->bodies[b].nodes.x_t_trial.array() * this->bodies[b].nodes.m.array();
                this->bodies[b].nodes.my_t = this->bodies[b].nodes.y_t_trial.array() * this->bodies[b].nodes.m.array();
                this->bodies[b].nodes.mz_t = this->bodies[b].nodes.z_t_trial.array() * this->bodies[b].nodes.m.array();

                this->bodies[b].nodes.x_t = this->bodies[b].nodes.x_t_trial;
                this->bodies[b].nodes.y_t = this->bodies[b].nodes.y_t_trial;
                this->bodies[b].nodes.z_t = this->bodies[b].nodes.z_t_trial;
            }

            //if failing to converge, break
            if(!isConverging){break;}

            //add contact forces to trial velocity
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
                this->bodies[b].nodes.fx_L = this->bodies[b].nodes.contact_fx;
                this->bodies[b].nodes.fy_L = this->bodies[b].nodes.contact_fy;
                this->bodies[b].nodes.fz_L = this->bodies[b].nodes.contact_fz;
            }

            this->calculateImplicitResidual();

            for (size_t b = 0; b < this->num_bodies; b++) {
                double vNorm = std::sqrt(
                        this->bodies[b].nodes.x_t_n.squaredNorm() + this->bodies[b].nodes.y_t_n.squaredNorm() +
                        this->bodies[b].nodes.z_t_n.squaredNorm());
                if (vNorm <= 0) {
                    vNorm = 1.0;
                }
                double sNorm = this->bodies[b].nodes.qk.norm();
                if (sNorm > 0) {
                    //isConverging = true;
                    this->bodies[b].nodes.DhRx = sNorm / (h * vNorm) * (this->bodies[b].nodes.Rx - this->bodies[b].nodes.Rvx);
                    this->bodies[b].nodes.DhRy = sNorm / (h * vNorm) * (this->bodies[b].nodes.Ry - this->bodies[b].nodes.Rvy);
                    this->bodies[b].nodes.DhRz = sNorm / (h * vNorm) * (this->bodies[b].nodes.Rz - this->bodies[b].nodes.Rvz);

                    this->bodies[b].nodes.tk << this->bodies[b].nodes.DhRx, this->bodies[b].nodes.DhRy, this->bodies[b].nodes.DhRz;
                    this->bodies[b].nodes.ok =
                            this->bodies[b].nodes.tk.dot(this->bodies[b].nodes.qk) / (this->bodies[b].nodes.tk.dot(this->bodies[b].nodes.tk));
                    this->bodies[b].nodes.sk = this->bodies[b].nodes.hk + this->bodies[b].nodes.ok * this->bodies[b].nodes.qk;

                    //check convergence of sk?

                    if (!std::isfinite(this->bodies[b].nodes.ok)) {
                        std::cout << "ok is infinite\n";
                    }

                    this->bodies[b].nodes.rk = this->bodies[b].nodes.qk - this->bodies[b].nodes.ok * this->bodies[b].nodes.tk;

                    //store if minimum
                    if (this->bodies[b].nodes.rk.norm() < this->bodies[b].nodes.rkMin.norm() || k==0){
                        this->bodies[b].nodes.rkMin = this->bodies[b].nodes.rk;
                        this->bodies[b].nodes.skMin = this->bodies[b].nodes.sk;
                    }
                } else {
                    isConverging = false;
                    /*std::cout << k << " 4 " << sNorm << std::endl;
                    //std::cout << "sNorm <= " << TOL << std::endl;
                    isConverging = false;
                    this->bodies[b].nodes.sk = this->bodies[b].nodes.hk;
                    this->bodies[b].nodes.rk = this->bodies[b].nodes.qk;
                    //let ok retain previous value for calculation of bk at beginning of step*/
                }
            }

            //if failing to converge, break
            if(!isConverging){break;}

            k += 1;
            rhoSum = 0;
            for (size_t b = 0; b < this->num_bodies; b++) {
                rhoSum += this->bodies[b].nodes.rk.squaredNorm();
            }
            std::cout << "\rn: " << nIter << " k: " << k <<  " r: " << rhoSum << "      \r" << std::flush;
        }

        for (size_t b = 0; b < this->num_bodies; b++) {
            /*this->bodies[b].nodes.x_t_n += this->bodies[b].nodes.sk.segment(0,this->bodies[b].n);
            this->bodies[b].nodes.y_t_n += this->bodies[b].nodes.sk.segment(this->bodies[b].n,this->bodies[b].n);
            this->bodies[b].nodes.z_t_n += this->bodies[b].nodes.sk.segment(2*this->bodies[b].n,this->bodies[b].n);*/

            this->bodies[b].nodes.x_t_n += this->bodies[b].nodes.skMin.segment(0,this->bodies[b].n);
            this->bodies[b].nodes.y_t_n += this->bodies[b].nodes.skMin.segment(this->bodies[b].n,this->bodies[b].n);
            this->bodies[b].nodes.z_t_n += this->bodies[b].nodes.skMin.segment(2*this->bodies[b].n,this->bodies[b].n);

            this->bodies[b].nodes.x_t_trial = this->bodies[b].nodes.x_t_n;
            this->bodies[b].nodes.y_t_trial = this->bodies[b].nodes.y_t_n;
            this->bodies[b].nodes.z_t_trial = this->bodies[b].nodes.z_t_n;

            this->bodies[b].nodes.mx_t = this->bodies[b].nodes.x_t_trial.array() * this->bodies[b].nodes.m.array();
            this->bodies[b].nodes.my_t = this->bodies[b].nodes.y_t_trial.array() * this->bodies[b].nodes.m.array();
            this->bodies[b].nodes.mz_t = this->bodies[b].nodes.z_t_trial.array() * this->bodies[b].nodes.m.array();

            this->bodies[b].nodes.x_t = this->bodies[b].nodes.x_t_trial;
            this->bodies[b].nodes.y_t = this->bodies[b].nodes.y_t_trial;
            this->bodies[b].nodes.z_t = this->bodies[b].nodes.z_t_trial;
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
            this->bodies[b].nodes.fx_L = this->bodies[b].nodes.contact_fx;
            this->bodies[b].nodes.fy_L = this->bodies[b].nodes.contact_fy;
            this->bodies[b].nodes.fz_L = this->bodies[b].nodes.contact_fz;
        }

        this->calculateImplicitResidual();

        //setup for next iteration
        rhoSum = 0;
        for (size_t b = 0; b < this->num_bodies; b++) {
            this->bodies[b].nodes.rk << -this->bodies[b].nodes.Rx, -this->bodies[b].nodes.Ry, -this->bodies[b].nodes.Rz;
            rhoSum += this->bodies[b].nodes.rk.squaredNorm();
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

    //add boddy forces
    time_varying_loads(this);

    //material stress update
    this->updateStress();

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

    //add boddy forces
    time_varying_loads(this);

    //material stress update
    this->updateStress();

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

    //add boddy forces
    time_varying_loads(this);

    //material stress update
    this->updateStress();

    return 1;
}

int job_t::mpmStepUSLImplicit2D() {
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
    this->moveParticlesImplicit2D();

    //calculate L on particles
    this->calculateStrainRate2D();

    //update particle densities
    this->updateDensity();

    //add boddy forces
    time_varying_loads(this);

    //material stress update
    this->updateStress();

    return 1;
}
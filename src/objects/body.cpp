//
// Created by aaron on 8/26/16.
// body.cpp
//
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <eigen3/Eigen/Sparse>

#include "body.hpp"
#include "particle.hpp"
#include "node.hpp"

/*#define WHICH_ELEMENT WHICH_ELEMENT9
//very ugly but it should work
#define WHICH_ELEMENT9(px,py,pz,Nx,Ny,Nz,Lx,Ly,Lz,hx,hy,hz) \
    ((int)(((px)<Lx && (px)>=0 && (py)<Ly && (py)>=0 && (pz)<Lz && (pz)>=0)?((floor((px)/(hx)) + floor((py)/(hy))*((Nx)-1) + floor((pz)/(hz))*((Nx*Ny)-1)):(-1))))
*/

Body::Body(size_t numNodes, size_t numParticles, size_t numElements, size_t bodyID):
        n(numNodes),
        p(numParticles),
        e(numElements),
        id(bodyID),

        //objects
        particles(p),
        nodes(n)
        //elements(e),

        //Phi and gradPhi
        /*Phi(n,p),
        gradPhiX(n,p),
        gradPhiY(n,p),
        gradPhiZ(n,p)*/
{
    material = Material();
    Phi.resize(n,p);
    gradPhiX.resize(n,p);
    gradPhiY.resize(n,p);
    gradPhiZ.resize(n,p);

    //Phi and gradPhi
    PhiTriplets.reserve(8*8*p);
    gradPhiXTriplets.reserve(8*8*p);
    gradPhiYTriplets.reserve(8*8*p);
    gradPhiZTriplets.reserve(8*8*p);

}

/*void Body::addParticle(double mIn,double vIn,double xIn,double yIn,double zIn,double x_tIn,double y_tIn,double z_tIn, size_t idIn){
    particles.addParticle(mIn,vIn,xIn,yIn,zIn,x_tIn,y_tIn,z_tIn,idIn);
}*/

/*void Body::addNode(double xIn, double yIn, double zIn, size_t idIn) {
    //Add node from job. Zero out unset terms
    nodes.addNode(xIn,yIn,xIn,idIn);
}*/

void Body::defineMaterial(double * fp64_props, size_t num_fp64_props , int * int_props, size_t num_int_props) {
    //this->material = Material();
    this->material.fp64_props = fp64_props;
    this->material.int_props = int_props;
    this->material.num_fp64_props = num_fp64_props;
    this->material.num_int_props = num_int_props;

    this->material.material_init(this);
    return;
}

//Body::~Body(){
//
//}
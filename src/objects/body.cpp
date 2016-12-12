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
        nodes(n),
        elements(e),

        //mass
        node_m(n),

        //position
        node_x(n),
        node_y(n),
        node_z(n),

        //displacement
        node_ux(n),
        node_uy(n),
        node_uz(n),

        //velocity
        node_x_t(n),
        node_y_t(n),
        node_z_t(n),

        //velocity difference
        node_diff_x_t(n),
        node_diff_y_t(n),
        node_diff_z_t(n),

        //momentum
        node_mx_t(n),
        node_my_t(n),
        node_mz_t(n),

        //force
        node_fx(n),
        node_fy(n),
        node_fz(n),

        //density
        node_rho(n),

        //body contact resolution
        node_contact_mx_t(n),
        node_contact_my_t(n),
        node_contact_mz_t(n),

        node_contact_x_t(n),
        node_contact_y_t(n),
        node_contact_z_t(n),

        node_contact_fx(n),
        node_contact_fy(n),
        node_contact_fz(n),

        //node_real_contact_fx(n),
        //node_real_contact_fy(n),
        //node_real_contact_fz(n),

        node_contact_normal_x(n),
        node_contact_normal_y(n),
        node_contact_normal_z(n),

        //implicit states
        node_mx_t_k(n),
        node_my_t_k(n),
        node_mz_t_k(n),

        node_x_t_trial(n),
        node_y_t_trial(n),
        node_z_t_trial(n),

        node_x_t_explicit(n),
        node_y_t_explicit(n),
        node_z_t_explicit(n),

        node_x_t_n(n),
        node_y_t_n(n),
        node_z_t_n(n),

        node_fx_k(n),
        node_fy_k(n),
        node_fz_k(n),

        node_fx_L(n),
        node_fy_L(n),
        node_fz_L(n),

        Rx(n),
        Ry(n),
        Rz(n),

        Rvx(n),
        Rvy(n),
        Rvz(n),

        DhRx(n),
        DhRy(n),
        DhRz(n),

        //ak(n),
        //rhok(n),
        //bk(n),
        wk(3*n),
        sk(3*n),
        pk(3*n),
        rk(3*n),
        r0(3*n),

        //position
        particle_x(p),
        particle_y(p),
        particle_z(p),

        //volume
        particle_v(p),
        particle_v_trial(p),
        particle_v0(p),
        particle_v_averaging(p),

        //half side length
        particle_a(p),

        //mass
        particle_m(p),

        //velocity
        particle_x_t(p),
        particle_y_t(p),
        particle_z_t(p),

        //body forces
        particle_bx(p),
        particle_by(p),
        particle_bz(p),

        //displacements
        particle_ux(p),
        particle_uy(p),
        particle_uz(p),

        //active
        particle_active(p)

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

void Body::addParticle(double mIn,double vIn,double xIn,double yIn,double zIn,double x_tIn,double y_tIn,double z_tIn, size_t idIn){
    //Add particle from file. Zero out unset terms.
    //position
    particle_x[idIn] = xIn;
    particle_y[idIn] = yIn;
    particle_z[idIn] = zIn;

    //volume
    particle_v[idIn] = vIn;
    particle_v_trial[idIn] = vIn;
    particle_v0[idIn] = vIn;
    particle_v_averaging[idIn] = 0.125*vIn; //from Sachith's code

    //half side length
    particle_a[idIn] = 0.5*cbrt(0.125*vIn); //from Sachith's code

    //mass
    particle_m[idIn] = mIn;

    //velocity
    particle_x_t[idIn] = x_tIn;
    particle_y_t[idIn] = y_tIn;
    particle_z_t[idIn] = z_tIn;

    //body forces
    particle_bx[idIn] = 0;
    particle_by[idIn] = 0;
    particle_bz[idIn] = 0;

    //displacements
    particle_ux[idIn] = 0;
    particle_uy[idIn] = 0;
    particle_uz[idIn] = 0;

    //active
    particle_active[idIn] = 0;

    //create Particle object
    this->particles[idIn] = Particle(this,idIn);

}

void Body::addNode(double xIn, double yIn, double zIn, size_t idIn) {
    //Add node from job. Zero out unset terms
    //mass
    node_m[idIn] = 0;

    //position
    node_x[idIn] = xIn;
    node_y[idIn] = yIn;
    node_z[idIn] = zIn;

    //displacement
    node_ux[idIn] = 0;
    node_uy[idIn] = 0;
    node_uz[idIn] = 0;

    //velocity
    node_x_t[idIn] = 0;
    node_y_t[idIn] = 0;
    node_z_t[idIn] = 0;

    //velocity difference
    node_diff_x_t[idIn] = 0;
    node_diff_y_t[idIn] = 0;
    node_diff_z_t[idIn] = 0;

    //momentum
    node_mx_t[idIn] = 0;
    node_my_t[idIn] = 0;
    node_mz_t[idIn] = 0;

    //force
    node_fx[idIn] = 0;
    node_fy[idIn] = 0;
    node_fz[idIn] = 0;

    //density
    node_rho[idIn] = 0;

    //body contact resolution
    node_contact_mx_t[idIn] = 0;
    node_contact_my_t[idIn] = 0;
    node_contact_mz_t[idIn] = 0;

    node_contact_x_t[idIn] = 0;
    node_contact_y_t[idIn] = 0;
    node_contact_z_t[idIn] = 0;

    node_contact_fx[idIn] = 0;
    node_contact_fy[idIn] = 0;
    node_contact_fz[idIn] = 0;

    //node_real_contact_fx[idIn] = 0;
    //node_real_contact_fy[idIn] = 0;
    //node_real_contact_fz[idIn] = 0;

    node_contact_normal_x[idIn] = 0;
    node_contact_normal_y[idIn] = 0;
    node_contact_normal_z[idIn] = 0;

    //implicit states
    node_x_t_trial[idIn] = 0;
    node_y_t_trial[idIn] = 0;
    node_z_t_trial[idIn] = 0;

    node_x_t_explicit[idIn] = 0;
    node_y_t_explicit[idIn] = 0;
    node_z_t_explicit[idIn] = 0;

    node_fx_k[idIn] = 0;
    node_fy_k[idIn] = 0;
    node_fz_k[idIn] = 0;

    node_fx_L[idIn] = 0;
    node_fy_L[idIn] = 0;
    node_fz_L[idIn] = 0;

    //nodal residuals
    Rx[idIn] = 0;
    Ry[idIn] = 0;
    Rz[idIn] = 0;

    DhRx[idIn] = 0;
    DhRy[idIn] = 0;
    DhRz[idIn] = 0;
    
    //implicit algorithm
    sk[idIn]=0;
    sk[idIn+n]=0;
    sk[idIn+2*n]=0;

    /*ak[idIn]=0;
    sk[idIn]=0;
    rk[idIn]=0;
    rhok[idIn]=0;
    bk[idIn]=0;
    pk[idIn]=0;*/

    //create Nodes object
    this->nodes[idIn] = Node(this,idIn);
}

void Body::addElement(size_t * nodeIDs, size_t idIn) {
    this->elements[idIn] = Element(8,nodeIDs,idIn);
    return;
}

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
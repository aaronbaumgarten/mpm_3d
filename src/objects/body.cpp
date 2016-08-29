//
// Created by aaron on 8/26/16.
// body.cpp
//
#include <iostream>
#include <stdlib.h>

#include "../../include/body.hpp"
#include "../../include/particle.hpp"
#include "../../include/node.hpp"

#define WHICH_ELEMENT WHICH_ELEMENT9
//very ugly but it should work
#define WHICH_ELEMENT9(px,py,pz,Nx,Ny,Nz,Lx,Ly,Lz,hx,hy,hz) \
    ((int)(((px)<Lx && (px)>=0 && (py)<Ly && (py)>=0 && (pz)<Lz && (pz)>=0)?((floor((px)/(hx)) + floor((py)/(hy))*((Nx)-1) + floor((pz)/(hz))*((Nx*Ny)-1)):(-1))))

Body::Body(size_t numNodes, size_t numParticles, size_t bodyID):
        n(numParticles),
        p(numParticles),
        id(bodyID),

        //objects
        particles(p),
        nodes(n),

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
        node_contact_x_t(n),
        node_contact_y_t(n),
        node_contact_z_t(n),

        node_contact_fx(n),
        node_contact_fy(n),
        node_contact_fz(n),

        node_real_contact_fx(n),
        node_real_contact_fy(n),
        node_real_contact_fz(n),

        node_contact_normal_x(n),
        node_contact_normal_y(n),
        node_contact_normal_z(n),

        //position
        particle_x(p),
        particle_y(p),
        particle_z(p),

        //volume
        particle_v(p),
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
{ }

void Body::addParticle(double mIn,double vIn,double xIn,double yIn,double zIn,double x_tIn,double y_tIn,double z_tIn, size_t idIn){
    //Add particle from file. Zero out unset terms.
    //position
    particle_x[idIn] = xIn;
    particle_y[idIn] = yIn;
    particle_z[idIn] = zIn;

    //volume
    particle_v[idIn] = vIn;
    particle_v0[idIn] = vIn;
    particle_v_averaging[idIn] = vIn;

    //half side length
    particle_a[idIn] = 0.25; //from Sachith's code

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

    //create references in Particle object
    //must be written this way to ensure that particle.hpp doesn't reference body.hpp
    this->particles[idIn] = Particle(this,idIn);
}

//Body::~Body(){
//
//}
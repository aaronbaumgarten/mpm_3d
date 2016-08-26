//
// Created by aaron on 8/26/16.
// body.cpp
//

#include <stdlib.h>

#include "../../include/body.hpp"
#include "../../include/particle.hpp"
#include "../../include/node.hpp"


Body::Body(size_t numNodes, size_t numParticles, size_t bodyID):
        n(numParticles),
        p(numParticles),
        id(bodyID),

        //objects
        particle(p),
        node(n),

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
        particle_uz(p)
{ }

//Body::~Body(){
//
//}
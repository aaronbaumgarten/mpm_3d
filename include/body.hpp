//
// Created by aaron on 8/26/16.
// body.hpp
//

#ifndef MPM_3D_BODY_HPP
#define MPM_3D_BODY_HPP

#include <vector>

#include "particle.hpp"
#include "node.hpp"

class Body{
public:
    //nodes
    size_t n;

    //particles
    size_t p;

    //body id
    size_t id;

    //***** objects *****
    std::vector<Particle> particles;
    std::vector<Node> nodes;

    //***** nodal vectors *****
    //mass
    std::vector<double> node_m;

    //position
    std::vector<double> node_x;
    std::vector<double> node_y;
    std::vector<double> node_z;

    //displacement
    std::vector<double> node_ux;
    std::vector<double> node_uy;
    std::vector<double> node_uz;

    //velocity
    std::vector<double> node_x_t;
    std::vector<double> node_y_t;
    std::vector<double> node_z_t;

    //velocity difference
    std::vector<double> node_diff_x_t;
    std::vector<double> node_diff_y_t;
    std::vector<double> node_diff_z_t;

    //momentum
    std::vector<double> node_mx_t;
    std::vector<double> node_my_t;
    std::vector<double> node_mz_t;

    //force
    std::vector<double> node_fx;
    std::vector<double> node_fy;
    std::vector<double> node_fz;

    //density
    std::vector<double> node_rho;

    //body contact resolution
    std::vector<double> node_contact_x_t;
    std::vector<double> node_contact_y_t;
    std::vector<double> node_contact_z_t;

    std::vector<double> node_contact_fx;
    std::vector<double> node_contact_fy;
    std::vector<double> node_contact_fz;

    std::vector<double> node_real_contact_fx;
    std::vector<double> node_real_contact_fy;
    std::vector<double> node_real_contact_fz;

    std::vector<double> node_contact_normal_x;
    std::vector<double> node_contact_normal_y;
    std::vector<double> node_contact_normal_z;

    //***** particle values *****
    //position
    std::vector<double> particle_x;
    std::vector<double> particle_y;
    std::vector<double> particle_z;

    //volume
    std::vector<double> particle_v;
    std::vector<double> particle_v0;
    std::vector<double> particle_v_averaging;

    //half side length
    std::vector<double> particle_a;

    //mass
    std::vector<double> particle_m;

    //velocity
    std::vector<double> particle_x_t;
    std::vector<double> particle_y_t;
    std::vector<double> particle_z_t;

    //body forces
    std::vector<double> particle_bx;
    std::vector<double> particle_by;
    std::vector<double> particle_bz;

    //displacements
    std::vector<double> particle_ux;
    std::vector<double> particle_uy;
    std::vector<double> particle_uz;

    //active
    std::vector<int> particle_active;

    //construcors
    Body(size_t,size_t,size_t);
    //destructors
    //~Body();

    //functions
    void addParticle(double,double,double,double,double,double,double,double,size_t);
    void addNode();
    void addElement();
    void mapParticles2Elements();
    void mapParticles2Nodes();
};

#endif //MPM_3D_BODY_HPP

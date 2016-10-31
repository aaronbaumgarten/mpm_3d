//
// Created by aaron on 8/26/16.
// body.hpp
//

#ifndef MPM_3D_BODY_HPP
#define MPM_3D_BODY_HPP

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

#include "particle.hpp"
#include "node.hpp"
#include "element.hpp"
#include "material.hpp"

class Body{
public:
    //material
    Material material;

    //nodes
    size_t n;

    //particles
    size_t p;

    //elements
    size_t e;

    //body id
    size_t id;

    //***** objects *****
    std::vector<Particle> particles;
    std::vector<Node> nodes;
    std::vector<Element> elements;

    //***** nodal vectors *****
    //mass
    Eigen::VectorXd node_m;

    //position
    Eigen::VectorXd node_x;
    Eigen::VectorXd node_y;
    Eigen::VectorXd node_z;

    //displacement
    Eigen::VectorXd node_ux;
    Eigen::VectorXd node_uy;
    Eigen::VectorXd node_uz;

    //velocity
    Eigen::VectorXd node_x_t;
    Eigen::VectorXd node_y_t;
    Eigen::VectorXd node_z_t;

    //velocity difference
    Eigen::VectorXd node_diff_x_t;
    Eigen::VectorXd node_diff_y_t;
    Eigen::VectorXd node_diff_z_t;

    //momentum
    Eigen::VectorXd node_mx_t;
    Eigen::VectorXd node_my_t;
    Eigen::VectorXd node_mz_t;

    //force
    Eigen::VectorXd node_fx;
    Eigen::VectorXd node_fy;
    Eigen::VectorXd node_fz;

    //density
    Eigen::VectorXd node_rho;

    //body contact resolution
    Eigen::VectorXd node_contact_mx_t;
    Eigen::VectorXd node_contact_my_t;
    Eigen::VectorXd node_contact_mz_t;

    Eigen::VectorXd node_contact_x_t;
    Eigen::VectorXd node_contact_y_t;
    Eigen::VectorXd node_contact_z_t;

    Eigen::VectorXd node_contact_fx;
    Eigen::VectorXd node_contact_fy;
    Eigen::VectorXd node_contact_fz;

    Eigen::VectorXd node_real_contact_fx;
    Eigen::VectorXd node_real_contact_fy;
    Eigen::VectorXd node_real_contact_fz;

    Eigen::VectorXd node_contact_normal_x;
    Eigen::VectorXd node_contact_normal_y;
    Eigen::VectorXd node_contact_normal_z;

    //***** particle values *****
    //position
    Eigen::VectorXd particle_x;
    Eigen::VectorXd particle_y;
    Eigen::VectorXd particle_z;

    //volume
    Eigen::VectorXd particle_v;
    Eigen::VectorXd particle_v0;
    Eigen::VectorXd particle_v_averaging;

    //half side length
    Eigen::VectorXd particle_a;

    //mass
    Eigen::VectorXd particle_m;

    //velocity
    Eigen::VectorXd particle_x_t;
    Eigen::VectorXd particle_y_t;
    Eigen::VectorXd particle_z_t;

    //body forces
    Eigen::VectorXd particle_bx;
    Eigen::VectorXd particle_by;
    Eigen::VectorXd particle_bz;

    //displacements
    Eigen::VectorXd particle_ux;
    Eigen::VectorXd particle_uy;
    Eigen::VectorXd particle_uz;

    //active
    Eigen::VectorXi particle_active;

    //S and gradS triplets
    std::vector<Eigen::Triplet<double>,Eigen::aligned_allocator<Eigen::Triplet<double>>> SipTriplets;
    std::vector<Eigen::Triplet<double>,Eigen::aligned_allocator<Eigen::Triplet<double>>> gradSipXTriplets;
    std::vector<Eigen::Triplet<double>,Eigen::aligned_allocator<Eigen::Triplet<double>>> gradSipYTriplets;
    std::vector<Eigen::Triplet<double>,Eigen::aligned_allocator<Eigen::Triplet<double>>> gradSipZTriplets;

    //Sip and gradSip
    Eigen::SparseMatrix<double> Sip;
    Eigen::SparseMatrix<double> gradSipX;
    Eigen::SparseMatrix<double> gradSipY;
    Eigen::SparseMatrix<double> gradSipZ;

    //construcors
    Body(size_t,size_t,size_t,size_t);
    //destructors
    //~Body();

    //functions
    void addParticle(double,double,double,double,double,double,double,double,size_t);
    void addNode(double,double,double,size_t);
    void addElement(size_t[8],size_t);
    void defineMaterial(double*,int*);
};

#endif //MPM_3D_BODY_HPP

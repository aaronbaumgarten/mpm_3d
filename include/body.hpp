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

    //Eigen::VectorXd node_real_contact_fx;
    //Eigen::VectorXd node_real_contact_fy;
    //Eigen::VectorXd node_real_contact_fz;

    Eigen::VectorXd node_contact_normal_x;
    Eigen::VectorXd node_contact_normal_y;
    Eigen::VectorXd node_contact_normal_z;

    //implicit states (not stored on nodes)
    Eigen::VectorXd node_mx_t_k;
    Eigen::VectorXd node_my_t_k;
    Eigen::VectorXd node_mz_t_k;

    Eigen::VectorXd node_x_t_trial;
    Eigen::VectorXd node_y_t_trial;
    Eigen::VectorXd node_z_t_trial;

    Eigen::VectorXd node_x_t_n;
    Eigen::VectorXd node_y_t_n;
    Eigen::VectorXd node_z_t_n;

    Eigen::VectorXd node_x_t_explicit;
    Eigen::VectorXd node_y_t_explicit;
    Eigen::VectorXd node_z_t_explicit;

    Eigen::VectorXd node_fx_k;
    Eigen::VectorXd node_fy_k;
    Eigen::VectorXd node_fz_k;

    Eigen::VectorXd node_fx_L;
    Eigen::VectorXd node_fy_L;
    Eigen::VectorXd node_fz_L;

    //nodal residuals (not stored on nodes)
    Eigen::VectorXd Rx;
    Eigen::VectorXd Ry;
    Eigen::VectorXd Rz;

    Eigen::VectorXd Rvx;
    Eigen::VectorXd Rvy;
    Eigen::VectorXd Rvz;

    Eigen::VectorXd DhRx;
    Eigen::VectorXd DhRy;
    Eigen::VectorXd DhRz;

    //implicit algorithm
    Eigen::VectorXd wk;
    //Eigen::VectorXd ak;//*
    double ak;
    Eigen::VectorXd sk;
    Eigen::VectorXd rk;
    Eigen::VectorXd r0;
    Eigen::VectorXd qk;
    Eigen::VectorXd tk;
    Eigen::VectorXd hk;
    double ok;
    //Eigen::VectorXd rhok;//*
    double rhok;
    //Eigen::VectorXd bk;//*
    double bk;
    Eigen::VectorXd pk;

    //***** particle values *****
    //position
    Eigen::VectorXd particle_x;
    Eigen::VectorXd particle_y;
    Eigen::VectorXd particle_z;

    //volume
    Eigen::VectorXd particle_v;
    Eigen::VectorXd particle_v_trial;
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
    std::vector<Eigen::Triplet<double>,Eigen::aligned_allocator<Eigen::Triplet<double>>> PhiTriplets;
    std::vector<Eigen::Triplet<double>,Eigen::aligned_allocator<Eigen::Triplet<double>>> gradPhiXTriplets;
    std::vector<Eigen::Triplet<double>,Eigen::aligned_allocator<Eigen::Triplet<double>>> gradPhiYTriplets;
    std::vector<Eigen::Triplet<double>,Eigen::aligned_allocator<Eigen::Triplet<double>>> gradPhiZTriplets;

    //Phi and gradPhi
    Eigen::SparseMatrix<double> Phi;
    Eigen::SparseMatrix<double> gradPhiX;
    Eigen::SparseMatrix<double> gradPhiY;
    Eigen::SparseMatrix<double> gradPhiZ;

    //construcors
    Body(size_t,size_t,size_t,size_t);
    //destructors
    //~Body();

    //functions
    void addParticle(double,double,double,double,double,double,double,double,size_t);
    void addNode(double,double,double,size_t);
    void addElement(size_t[8],size_t);
    void defineMaterial(double*,size_t,int*,size_t);
};

#endif //MPM_3D_BODY_HPP

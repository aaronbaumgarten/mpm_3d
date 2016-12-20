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
    Nodes nodes;

    //***** nodal vectors *****
    // stored on nodes now

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
    //void addNode(double,double,double,size_t);
    void defineMaterial(double*,size_t,int*,size_t);
};

#endif //MPM_3D_BODY_HPP

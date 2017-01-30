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
    Particles particles;
    Nodes nodes;

    //***** nodal vectors *****
    // stored on nodes now

    //***** particle values *****
    // stored on particles now

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
    //void addParticle(double,double,double,double,double,double,double,double,size_t);
    //void addNode(double,double,double,size_t);
    void defineMaterial(std::string,std::string,std::vector<double>,std::vector<int>);
};

#endif //MPM_3D_BODY_HPP

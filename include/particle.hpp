//
// Created by aaron on 8/26/16.
// particle.hpp (heavily borrows from mpm-2d-legacy)
//

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <Eigen/Core>

#ifndef MPM_3D_PARTICLE_HPP
#define MPM_3D_PARTICLE_HPP

//depth of state vector
#define DEPVAR 11

//shape function related indices in arrays
#define S_XIDX 0
#define S_YIDX 1
#define S_ZIDX 2

//degress of freedom
#define NODAL_DOF 3
#define XDOF_IDX 0
#define YDOF_IDX 1
#define ZDOF_IDX 2

//dimensions of the stress/strain tensors
#define NDIM 3

//indices for stress/strain tensors
#define XX (NDIM * S_XIDX + S_XIDX)
#define XY (NDIM * S_XIDX + S_YIDX)
#define XZ (NDIM * S_XIDX + S_ZIDX)
#define YX (NDIM * S_YIDX + S_XIDX)
#define YY (NDIM * S_YIDX + S_YIDX)
#define YZ (NDIM * S_YIDX + S_ZIDX)
#define ZX (NDIM * S_ZIDX + S_XIDX)
#define ZY (NDIM * S_ZIDX + S_YIDX)
#define ZZ (NDIM * S_ZIDX + S_ZIDX)

#define WHICH_ELEMENT WHICH_ELEMENT9
//very ugly but it should work
#define WHICH_ELEMENT9(px,py,pz,Nx,Ny,Nz,Lx,Ly,Lz,hx,hy,hz) \
    ((int)(((px)<Lx && (px)>=0 && (py)<Ly && (py)>=0 && (pz)<Lz && (pz)>=0)?(floor((px)/(hx)) + floor((py)/(hy))*((Nx-1)) + floor((pz)/(hz))*((Nx-1)*(Ny-1))):(-1)))

class Particles {
public:
    //number of particles
    size_t numParticles;

    //position
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd z;

    //volume
    Eigen::VectorXd v;
    Eigen::VectorXd v_trial;
    Eigen::VectorXd v0;
    Eigen::VectorXd v_averaging;

    //half side length
    Eigen::VectorXd a;

    //mass
    Eigen::VectorXd m;

    //velocity
    Eigen::VectorXd x_t;
    Eigen::VectorXd y_t;
    Eigen::VectorXd z_t;

    //body forces
    Eigen::VectorXd bx;
    Eigen::VectorXd by;
    Eigen::VectorXd bz;

    //full 3d stress tensor
    /* CAUTION: T is also used for templates */
    //double T[NDIM * NDIM];
    //double Ttrial[NDIM * NDIM];
    Eigen::MatrixXd T;
    Eigen::MatrixXd Ttrial;

    //velocity gradient
    //double L[NDIM*NDIM];
    Eigen::MatrixXd L;

    //full 3d deformation gradient
    //double F[NDIM * NDIM];
    Eigen::MatrixXd F;

    //full 3d plastic deformation gradient
    //double Fp[NDIM * NDIM];
    Eigen::MatrixXd Fp;

    //Left Cauchy Green Tensor
    Eigen::MatrixXd Be;

    //displacements
    Eigen::VectorXd ux;
    Eigen::VectorXd uy;
    Eigen::VectorXd uz;


    //state variables (used for constitutive laws)
    //double state[DEPVAR];
    Eigen::MatrixXd state;

    //flag for whether particle is active or not
    Eigen::VectorXi active;

    //matrix of corner positions corner[id][#][x, y or z]
    //double corner[8][3];
    Eigen::MatrixXd corner_x;
    Eigen::MatrixXd corner_y;
    Eigen::MatrixXd corner_z;

    //which element owns each corner
    //int corner_elements[8];
    Eigen::MatrixXi corner_elements;

    //construcors
    Particles(size_t p);
    Particles(){}

    //functions
    void addParticle(double,double,double,double,double,double,double,double,size_t);

    template<class jobT>
    void updateCorners(jobT* job, size_t id);

    template<class jobT>
    void updateAllCorners(jobT* job);

    template<class jobT>
    void resetCorners(jobT* job, size_t id);

    template<class jobT>
    void resetAllCorners(jobT* job);

    template<class jobT>
    int updateActive(jobT* job, size_t id);
};

#endif //MPM_3D_PARTICLE_HPP
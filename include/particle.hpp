//
// Created by aaron on 8/26/16.
// particle.hpp (heavily borrows from mpm-2d-legacy)
//

#ifndef MPM_3D_PARTICLE_HPP
#define MPM_3D_PARTICLE_HPP

//depth of state vector
#define DEPVAR 11

//shape function related indices in arrays
#define S_XIDX 0
#define S_YIDX 1
#define S_ZIDX 2

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

class Particle{
public:
    //position
    double x;
    double y;
    double z;

    //volume
    double v;
    double v0;
    double v_averaging;

    //half side length
    double a;

    //mass
    double m;

    //velocity
    double x_t;
    double y_t;
    double z_t;

    //body forces
    double bx;
    double by;
    double bz;

    //full 3d stress tensor
    /* CAUTION: T is also used for templates */
    double T[NDIM*NDIM];

    //full 3d deformation gradient
    double F[NDIM*NDIM];

    //full 3d plastic deformation gradient
    double Fp[NDIM*NDIM];

    //displacements
    double ux;
    double uy;
    double uz;

    //state variables (used for constitutive laws)
    double state[DEPVAR];

    //flag for whether particle is active or not
    int active;

    //body type (for contact models)?
    //size_t body;

    //matrix of corner positions corner[#][x, y or z]
    double corner[8][3];

    //which element owns each corner
    int corner_elements[8];

    //unique id for particle tracking
    size_t id;

    double *working;

    //material specific data structure?
    //void *material_data;

    size_t blocksize;

    //construcors
    Particle() {}
    //destructors
    ~Particle() {}
};

#endif //MPM_3D_PARTICLE_HPP
//
// Created by aaron on 8/26/16.
// particle.hpp (heavily borrows from mpm-2d-legacy)
//

#include <iostream>
#include <stdlib.h>

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

class Particle {
public:
    //dummy variables
    double rmDouble = 0.0;
    int rmInt = 0;

    //position

    double* x;
    double* y;
    double* z;

    //volume
    double* v;
    double* v0;
    double* v_averaging;

    //half side length
    double* a;

    //mass
    double* m;

    //velocity
    double* x_t;
    double* y_t;
    double* z_t;

    //body forces
    double* bx;
    double* by;
    double* bz;

    //full 3d stress tensor
    /* CAUTION: T is also used for templates */
    double T[NDIM * NDIM];

    //full 3d deformation gradient
    double F[NDIM * NDIM];

    //full 3d plastic deformation gradient
    double Fp[NDIM * NDIM];

    //displacements
    double* ux;
    double* uy;
    double* uz;


    //state variables (used for constitutive laws)
    double state[DEPVAR];

    //flag for whether particle is active or not
    int* active;

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
    template<class T>
    Particle(T* bd, size_t idIn):
            id(idIn),
            active(&(bd->particle_active[idIn])),

            x(&(bd->particle_x[idIn])),
            y(&(bd->particle_y[idIn])),
            z(&(bd->particle_z[idIn])),

            //volume
            v(&(bd->particle_v[idIn])),
            v0(&(bd->particle_v0[idIn])),
            v_averaging(&(bd->particle_v_averaging[idIn])),

            //half side length
            a(&(bd->particle_a[idIn])),

            //mass
            m(&(bd->particle_m[idIn])),

            //velocity
            x_t(&(bd->particle_x_t[idIn])),
            y_t(&(bd->particle_y_t[idIn])),
            z_t(&(bd->particle_z_t[idIn])),

            //body forces
            bx(&(bd->particle_bx[idIn])),
            by(&(bd->particle_by[idIn])),
            bz(&(bd->particle_bz[idIn])),

            //displacements
            ux(&(bd->particle_ux[idIn])),
            uy(&(bd->particle_uy[idIn])),
            uz(&(bd->particle_uz[idIn]))
    {};

    Particle() {}
};

#endif //MPM_3D_PARTICLE_HPP
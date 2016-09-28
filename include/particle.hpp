//
// Created by aaron on 8/26/16.
// particle.hpp (heavily borrows from mpm-2d-legacy)
//

#include <iostream>
#include <stdlib.h>
#include <math.h>

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

#define WHICH_ELEMENT WHICH_ELEMENT9
//very ugly but it should work
#define WHICH_ELEMENT9(px,py,pz,Nx,Ny,Nz,Lx,Ly,Lz,hx,hy,hz) \
    ((int)(((px)<Lx && (px)>=0 && (py)<Ly && (py)>=0 && (pz)<Lz && (pz)>=0)?(floor((px)/(hx)) + floor((py)/(hy))*((Nx)-1) + floor((pz)/(hz))*((Nx-1)*(Ny-1))):(-1)))

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

    //velocity gradient
    double L[NDIM*NDIM];

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
    template<class bodyT>
    Particle(bodyT* bd, size_t idIn):
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
    {
        for (int i=0; i<NDIM*NDIM; i++){
            T[i] = 0;
            L[i] = 0;
            F[i] = 0;
            Fp[i] = 0;
        }
        F[XX] = 1;
        F[YY] = 1;
        F[ZZ] = 1;

        Fp[XX] = 1;
        Fp[YY] = 1;
        Fp[ZZ] = 1;

        state[9] = 0;
        state[10] = 0;

        for (int i=0; i<8; i++){
            for (int j=0; j<3; j++){
                corner[i][j] = 0;
            }
        }
    };

    Particle() {}

    //functions
    template<class jobT>
    void updateCorners(jobT* job){
        int sx = 1;
        int sy = 1;
        int sz = 1;
        ///implement corner math really just to be clever
        for (size_t i=0;i<8;i++){
            if (i%4 == 0){
                sx *= -1;
            }
            if (i%2 == 0){
                sy *= -1;
            }
            sz *= -1;
            corner[i][0] = x[0]+sx*a[0];
            corner[i][1] = y[0]+sy*a[0];
            corner[i][2] = z[0]+sz*a[0];

            //assumes regular cartesian grid
            corner_elements[i] = WHICH_ELEMENT9(corner[i][0],corner[i][1],corner[i][2],
                                                job->Nx,job->Ny,job->Nz,
                                                job->Lx,job->Ly,job->Lz,
                                                job->hx,job->hy,job->hz);
        }
    }

    template<class jobT>
    void resetCorners(jobT* job){
        for (size_t i=0;i<8;i++){
            //resets corners to particle position
            corner[i][0] = x[0];
            corner[i][1] = y[0];
            corner[i][2] = z[0];
            corner_elements[i] = WHICH_ELEMENT9(corner[i][0],corner[i][1],corner[i][2],
                                                job->Nx,job->Ny,job->Nz,
                                                job->Lx,job->Ly,job->Lz,
                                                job->hx,job->hy,job->hz);
        }
    }

    template<class jobT>
    int updateActive(jobT* job){
        if (WHICH_ELEMENT9(x[0],y[0],z[0],
                       job->Nx,job->Ny,job->Nz,
                       job->Lx,job->Ly,job->Lz,
                       job->hx,job->hy,job->hz) == -1) {
            active[0] = 0;
        } else {
            active[0] = 1;
        }
        return active[0];
    }
};

#endif //MPM_3D_PARTICLE_HPP
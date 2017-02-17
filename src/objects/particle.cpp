//
// Created by aaron on 12/20/16.
// particle.cpp
//

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <eigen3/Eigen/Sparse>

#include "particle.hpp"
#include "body.hpp"
#include "process.hpp"

Particles::Particles(size_t p):
        //number of particles
        numParticles(p),

        //position
        x(p),
        y(p),
        z(p),

        //volume
        v(p),
        v_trial(p),
        v0(p),
        v_averaging(p),

        //half side length
        a(p),

        //mass
        m(p),

        //velocity
        x_t(p),
        y_t(p),
        z_t(p),

        //body forces
        bx(p),
        by(p),
        bz(p),

        //displacements
        ux(p),
        uy(p),
        uz(p),

        //active
        active(p),

        //stress
        T(p,9),
        Ttrial(p,9),

        //velocity gradient
        L(p,9),

        //deformation gradients
        F(p,9),
        Fp(p,9),
        Be(p,9),

        //state
        state(p,DEPVAR),

        //corner positions
        corner_x(p,8),
        corner_y(p,8),
        corner_z(p,8),

        //corner elements
        corner_elements(p,8),
        elementIDs(p)
{
    //position
    x.setZero();
    y.setZero();
    z.setZero();

    //volume
    v.setZero();
    v_trial.setZero();
    v0.setZero();
    v_averaging.setZero();

    //half side length
    a.setZero();

    //mass
    m.setZero();

    //velocity
    x_t.setZero();
    y_t.setZero();
    z_t.setZero();

    //body forces
    bx.setZero();
    by.setZero();
    bz.setZero();

    //displacements
    ux.setZero();
    uy.setZero();
    uz.setZero();

    //active
    active.setZero();

    //stress
    T.setZero();
    Ttrial.setZero();

    //velocity gradient
    L.setZero();

    //deformation gradients
    F.setZero();
    Fp.setZero();
    Be.setZero();

    //state
    state.setZero();

    //corner positions
    corner_x.setZero();
    corner_y.setZero();
    corner_z.setZero();

    //corner elements
    corner_elements.setZero();
    elementIDs.setZero();
};

//functions
void Particles::addParticle(double mIn,double vIn,double xIn,double yIn,double zIn,double x_tIn,double y_tIn,double z_tIn, size_t idIn){
    //Add particle from file. Zero out unset terms.
    //position
    this->x[idIn] = xIn;
    this->y[idIn] = yIn;
    this->z[idIn] = zIn;

    //volume
    this->v[idIn] = vIn;
    this->v_trial[idIn] = vIn;
    this->v0[idIn] = vIn;
    this->v_averaging[idIn] = 0.125*vIn; //from Sachith's code

    //half side length
    this->a[idIn] = 0.5*cbrt(0.125*vIn); //from Sachith's code

    //mass
    this->m[idIn] = mIn;

    //velocity
    this->x_t[idIn] = x_tIn;
    this->y_t[idIn] = y_tIn;
    this->z_t[idIn] = z_tIn;

    //body forces
    this->bx[idIn] = 0;
    this->by[idIn] = 0;
    this->bz[idIn] = 0;

    //displacements
    this->ux[idIn] = 0;
    this->uy[idIn] = 0;
    this->uz[idIn] = 0;

    //active
    this->active[idIn] = 0;
}

template<>
void Particles::updateCorners(job_t* job, size_t id){
    int S[3][8] = {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
    for (size_t i=0;i<8;i++) {
        this->corner_x(id, i) = this->x[id] + S[0][i] * this->a[id];
        this->corner_y(id, i) = this->y[id] + S[1][i] * this->a[id];
        this->corner_z(id, i) = this->z[id] + S[2][i] * this->a[id];
        this->corner_elements(id, i) = WHICH_ELEMENT9(this->corner_x(id, i), this->corner_y(id, i),
                                                      this->corner_z(id, i),
                                                      job->Nx, job->Ny, job->Nz,
                                                      job->Lx, job->Ly, job->Lz,
                                                      job->hx, job->hy, job->hz);
    }
    return;
}

template<>
void Particles::updateAllCorners(job_t* job){
    int S[3][8] = {{-1,1,-1,1,-1,1,-1,1},{-1,-1,1,1,-1,-1,1,1},{-1,-1,-1,-1,1,1,1,1}};
    for (size_t p=0;p<this->numParticles;p++){
        for (size_t i=0;i<8;i++) {
            this->corner_x(p, i) = this->x[p] + S[0][i] * this->a[p];
            this->corner_y(p, i) = this->y[p] + S[1][i] * this->a[p];
            this->corner_z(p, i) = this->z[p] + S[2][i] * this->a[p];
            this->corner_elements(p, i) = WHICH_ELEMENT9(this->corner_x(p, i), this->corner_y(p, i),
                                                         this->corner_z(p, i),
                                                         job->Nx, job->Ny, job->Nz,
                                                         job->Lx, job->Ly, job->Lz,
                                                         job->hx, job->hy, job->hz);
        }
    }
    return;
}

template<>
void Particles::resetCorners(job_t* job, size_t id){
    for (size_t i = 0; i < 8; i++) {
        //resets corners to particle position
        this->corner_x(id,i) = this->x[id];
        this->corner_y(id,i) = this->y[id];
        this->corner_z(id,i) = this->z[id];
        this->corner_elements(id,i) = WHICH_ELEMENT9(this->corner_x(id,i), this->corner_y(id,i), this->corner_z(id,i),
                                                    job->Nx, job->Ny, job->Nz,
                                                    job->Lx, job->Ly, job->Lz,
                                                    job->hx, job->hy, job->hz);
    }
    return;
}

template<>
void Particles::resetAllCorners(job_t* job){
    for (size_t p=0;p<this->numParticles;p++) {
        for (size_t i = 0; i < 8; i++) {
            //resets corners to particle position
            this->corner_x(p,i) = this->x[p];
            this->corner_y(p,i) = this->y[p];
            this->corner_z(p,i) = this->z[p];
            this->corner_elements(p,i) = WHICH_ELEMENT9(this->corner_x(p,i), this->corner_y(p,i), this->corner_z(p,i),
                                                        job->Nx, job->Ny, job->Nz,
                                                        job->Lx, job->Ly, job->Lz,
                                                        job->hx, job->hy, job->hz);
        }
    }
    return;
}

template<>
int Particles::updateActive(job_t* job, size_t id){
    if (WHICH_ELEMENT9(this->x[id],this->y[id],this->z[id],
                       job->Nx,job->Ny,job->Nz,
                       job->Lx,job->Ly,job->Lz,
                       job->hx,job->hy,job->hz) == -1) {
        this->active[id] = 0;
    } else {
        this->active[id] = 1;
    }
    return this->active[id];
}

template<>
void Particles::updateElementIDs(job_t *job) {
    for (size_t p=0;p<this->numParticles;p++){
        this->elementIDs[p] = WHICH_ELEMENT9(this->x[p], this->y[p], this->z[p],
                                             job->Nx, job->Ny, job->Nz,
                                             job->Lx, job->Ly, job->Lz,
                                             job->hx, job->hy, job->hz);
    }
    return;
}
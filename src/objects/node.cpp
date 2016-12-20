//
// Created by aaron on 12/20/16.
// node.cpp
//

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <Eigen/Core>
#include "node.hpp"

Nodes::Nodes(size_t n):
        num_nodes(n),

        //mass
        m(n),

        //position
        x(n),
        y(n),
        z(n),

        //displacement
        ux(n),
        uy(n),
        uz(n),

        //velocity
        x_t(n),
        y_t(n),
        z_t(n),

        //velocity difference
        diff_x_t(n),
        diff_y_t(n),
        diff_z_t(n),

        //momentum
        mx_t(n),
        my_t(n),
        mz_t(n),

        //force
        fx(n),
        fy(n),
        fz(n),

        //density
        rho(n),

        //body contact resolution
        contact_mx_t(n),
        contact_my_t(n),
        contact_mz_t(n),

        contact_x_t(n),
        contact_y_t(n),
        contact_z_t(n),

        contact_fx(n),
        contact_fy(n),
        contact_fz(n),

        contact_normal_x(n),
        contact_normal_y(n),
        contact_normal_z(n),

        //nodal initial momentum
        mx_t_k(n),
        my_t_k(n),
        mz_t_k(n),

        //nodal trial velocity
        x_t_trial(n),
        y_t_trial(n),
        z_t_trial(n),

        //nodal iterated velocity
        x_t_n(n),
        y_t_n(n),
        z_t_n(n),

        //nodal initial force
        //fx_k(n),
        //fy_k(n),
        //fz_k(n),

        //nodal final force
        fx_L(n),
        fy_L(n),
        fz_L(n),

        //nodal residuals
        Rx(n),
        Ry(n),
        Rz(n),

        //saved residuals
        Rvx(n),
        Rvy(n),
        Rvz(n),

        //directional residual derivative
        DhRx(n),
        DhRy(n),
        DhRz(n),

        //implicit algorithm
        wk(3*n),
        ak(1),//*
        sk(3*n),
        rk(3*n),
        r0(3*n),
        qk(3*n),
        tk(3*n),
        hk(3*n),
        ok(1),//*
        rhok(1),//*
        bk(1),//*
        pk(3*n)
{ }

void Nodes::addNode(double xIn, double yIn, double zIn, size_t idIn) {
    //Add node from job. Zero out unset terms
    //mass
    this->m[idIn] = 0;

    //position
    this->x[idIn] = xIn;
    this->y[idIn] = yIn;
    this->z[idIn] = zIn;

    //displacement
    this->ux[idIn] = 0;
    this->uy[idIn] = 0;
    this->uz[idIn] = 0;

    //velocity
    this->x_t[idIn] = 0;
    this->y_t[idIn] = 0;
    this->z_t[idIn] = 0;

    //velocity difference
    this->diff_x_t[idIn] = 0;
    this->diff_y_t[idIn] = 0;
    this->diff_z_t[idIn] = 0;

    //momentum
    this->mx_t[idIn] = 0;
    this->my_t[idIn] = 0;
    this->mz_t[idIn] = 0;

    //force
    this->fx[idIn] = 0;
    this->fy[idIn] = 0;
    this->fz[idIn] = 0;

    //density
    this->rho[idIn] = 0;

    //body contact resolution
    this->contact_mx_t[idIn] = 0;
    this->contact_my_t[idIn] = 0;
    this->contact_mz_t[idIn] = 0;

    this->contact_x_t[idIn] = 0;
    this->contact_y_t[idIn] = 0;
    this->contact_z_t[idIn] = 0;

    this->contact_fx[idIn] = 0;
    this->contact_fy[idIn] = 0;
    this->contact_fz[idIn] = 0;

    //this->real_contact_fx[idIn] = 0;
    //this->real_contact_fy[idIn] = 0;
    //this->real_contact_fz[idIn] = 0;

    this->contact_normal_x[idIn] = 0;
    this->contact_normal_y[idIn] = 0;
    this->contact_normal_z[idIn] = 0;

    //implicit states
    this->x_t_trial[idIn] = 0;
    this->y_t_trial[idIn] = 0;
    this->z_t_trial[idIn] = 0;

    //this->fx_k[idIn] = 0;
    //this->fy_k[idIn] = 0;
    //this->fz_k[idIn] = 0;

    this->fx_L[idIn] = 0;
    this->fy_L[idIn] = 0;
    this->fz_L[idIn] = 0;

    //nodal residuals
    this->Rx[idIn] = 0;
    this->Ry[idIn] = 0;
    this->Rz[idIn] = 0;

    this->DhRx[idIn] = 0;
    this->DhRy[idIn] = 0;
    this->DhRz[idIn] = 0;

    //implicit algorithm
    //this->sk[idIn]=0;
    //this->sk[idIn+n]=0;
    //this->sk[idIn+2*n]=0;

    /*ak[idIn]=0;
    sk[idIn]=0;
    rk[idIn]=0;
    rhok[idIn]=0;
    bk[idIn]=0;
    pk[idIn]=0;*/
}
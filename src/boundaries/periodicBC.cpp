//
// Created by aaron on 2/6/17.
// periodicBC.cpp
//

//
// Created by aaron on 2/6/17.
// periodicBC2d.cpp
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "process.hpp"
#include "tensor.hpp"
#include "element.hpp"

extern "C" void generate_dirichlet_bcs(job_t *job);
extern "C" void generate_node_number_override(job_t *job);

extern "C" void bc_init(job_t *job)
{
    // suppress unused argument warning
    //(void)job;
    bool doubleFilled = false;
    for (size_t b=0;b<job->num_bodies;b++){
        if (doubleFilled){
            break;
        }
        for (size_t i=0;i<job->bodies[b].p;i++){
            if (job->bodies[b].particles.x[i] > (job->Lx - job->hx/2.0)){
                //element double filled
                doubleFilled = true;
                break;
            } else if (job->bodies[b].particles.x[i] < (job->hx/2.0)){
                //element double filled
                doubleFilled = true;
                break;
            }
        }
    }

    if (doubleFilled){
        std::cout << "Boundary condition: WARNING - bounding elements appear double filled.\n";
    }

    return;
}

extern "C" int bc_validate(job_t *job)
{
    //return 0 for failure to validate properties in configuration file
    //not implemented as of 10/26/16
    return 1;
}

extern "C" void bc_time_varying(job_t *job)
{
    generate_dirichlet_bcs(job);
    generate_node_number_override(job);
    return;
}

/*----------------------------------------------------------------------------*/

void generate_dirichlet_bcs(job_t *job)
{
    // set dirichlet mask to all zeros
    for (size_t i = 0; i < job->num_nodes; i++) {
        for (size_t j = 0; j < NODAL_DOF; j++) {
            job->u_dirichlet_mask[NODAL_DOF * i + j] = 0;
        }
    }

    // set floor value (0 to Nx*Ny-1)
    for (size_t nx = 0; nx < job->Nx; nx++) {
        for (size_t ny = 0; ny < job->Ny; ny++) {
            job->u_dirichlet[NODAL_DOF * nx + NODAL_DOF * job->Nx * ny + ZDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * nx + NODAL_DOF * job->Nx * ny + ZDOF_IDX] = 1;

            job->u_dirichlet[NODAL_DOF * nx + NODAL_DOF * job->Nx * ny + XDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * nx + NODAL_DOF * job->Nx * ny + XDOF_IDX] = 1;

            job->u_dirichlet[NODAL_DOF * nx + NODAL_DOF * job->Nx * ny + YDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * nx + NODAL_DOF * job->Nx * ny + YDOF_IDX] = 1;
        }
    }

    double mx_t, my_t, mz_t;
    double m, m1, m2;
    double fx, fy, fz;
    // set y-wall value (0 to Nx and Nx*Ny to Nx*Ny + Nx)
    for (size_t nx = 0; nx < job->Nx; nx++) {
        for (size_t nz = 0; nz<job->Nz;nz++) {
            /*job->u_dirichlet[NODAL_DOF * nx + NODAL_DOF * (job->Nx * job->Ny * nz) + YDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * nx + NODAL_DOF * (job->Nx * job->Ny * nz) + YDOF_IDX] = 1;

            job->u_dirichlet[NODAL_DOF * (nx + job->Nx * (job->Ny-1)) + NODAL_DOF * (job->Nx * job->Ny * nz) + YDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * (nx + job->Nx * (job->Ny-1)) + NODAL_DOF * (job->Nx * job->Ny * nz) + YDOF_IDX] = 1;*/
            for (size_t b = 0; b < job->num_bodies; b++) {
                //set lhs to rhs-1
                mx_t = 0;
                my_t = 0;
                mz_t = 0;
                m = 0;
                fx = 0;
                fy = 0;
                fz = 0;

                mx_t += job->bodies[b].nodes.contact_mx_t[nx + job->Nx * job->Ny * nz];
                my_t += job->bodies[b].nodes.contact_my_t[nx + job->Nx * job->Ny * nz];
                mz_t += job->bodies[b].nodes.contact_mz_t[nx + job->Nx * job->Ny * nz];
                mx_t += job->bodies[b].nodes.contact_mx_t[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz];
                my_t += job->bodies[b].nodes.contact_my_t[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz];
                mz_t += job->bodies[b].nodes.contact_mz_t[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz];

                fx += job->bodies[b].nodes.contact_fx[nx + job->Nx * job->Ny * nz];
                fy += job->bodies[b].nodes.contact_fy[nx + job->Nx * job->Ny * nz];
                fz += job->bodies[b].nodes.contact_fz[nx + job->Nx * job->Ny * nz];
                fx += job->bodies[b].nodes.contact_fx[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz];
                fy += job->bodies[b].nodes.contact_fy[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz];
                fz += job->bodies[b].nodes.contact_fz[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz];

                m1 = job->bodies[b].nodes.m[nx + job->Nx * job->Ny * nz];
                m2 = job->bodies[b].nodes.m[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz];
                m = m1 + m2;

                job->bodies[b].nodes.contact_mx_t[nx + job->Nx * job->Ny * nz] = mx_t * m1 / m;
                job->bodies[b].nodes.contact_my_t[nx + job->Nx * job->Ny * nz] = my_t * m1 / m;
                job->bodies[b].nodes.contact_mz_t[nx + job->Nx * job->Ny * nz] = mz_t * m1 / m;
                job->bodies[b].nodes.contact_mx_t[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz] =
                        mx_t * m2 / m;
                job->bodies[b].nodes.contact_my_t[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz] =
                        my_t * m2 / m;
                job->bodies[b].nodes.contact_mz_t[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz] =
                        mz_t * m2 / m;

                //job->bodies[b].nodes.m[nx + job->Nx * job->Ny * nz] = m;
                //job->bodies[b].nodes.m[nx + job->Nx * (job->Ny-2) + job->Nx * job->Ny * nz] = m;

                job->bodies[b].nodes.contact_fx[nx + job->Nx * job->Ny * nz] = fx * m1 / m;
                job->bodies[b].nodes.contact_fy[nx + job->Nx * job->Ny * nz] = fy * m1 / m;
                job->bodies[b].nodes.contact_fz[nx + job->Nx * job->Ny * nz] = fz * m1 / m;
                job->bodies[b].nodes.contact_fx[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz] = fx * m2 / m;
                job->bodies[b].nodes.contact_fy[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz] = fy * m2 / m;
                job->bodies[b].nodes.contact_fz[nx + job->Nx * (job->Ny - 2) + job->Nx * job->Ny * nz] = fz * m2 / m;

                //set lhs+1 to rhs
                mx_t = 0;
                my_t = 0;
                mz_t = 0;
                m = 0;
                fx = 0;
                fy = 0;
                fz = 0;

                mx_t += job->bodies[b].nodes.contact_mx_t[nx + job->Nx * (job->Ny + 1) * nz];
                my_t += job->bodies[b].nodes.contact_my_t[nx + job->Nx * (job->Ny + 1) * nz];
                mz_t += job->bodies[b].nodes.contact_mz_t[nx + job->Nx * (job->Ny + 1) * nz];
                mx_t += job->bodies[b].nodes.contact_mx_t[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz];
                my_t += job->bodies[b].nodes.contact_my_t[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz];
                mz_t += job->bodies[b].nodes.contact_mz_t[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz];

                fx += job->bodies[b].nodes.contact_fx[nx + job->Nx * (job->Ny + 1) * nz];
                fy += job->bodies[b].nodes.contact_fy[nx + job->Nx * (job->Ny + 1) * nz];
                fz += job->bodies[b].nodes.contact_fz[nx + job->Nx * (job->Ny + 1) * nz];
                fx += job->bodies[b].nodes.contact_fx[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz];
                fy += job->bodies[b].nodes.contact_fy[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz];
                fz += job->bodies[b].nodes.contact_fz[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz];

                m1 = job->bodies[b].nodes.m[nx + job->Nx * job->Ny * nz];
                m2 = job->bodies[b].nodes.m[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz];
                m = m1 + m2;

                job->bodies[b].nodes.contact_mx_t[nx + job->Nx * (job->Ny + 1) * nz] = mx_t * m1 / m;
                job->bodies[b].nodes.contact_my_t[nx + job->Nx * (job->Ny + 1) * nz] = my_t * m1 / m;
                job->bodies[b].nodes.contact_mz_t[nx + job->Nx * (job->Ny + 1) * nz] = mz_t * m1 / m;
                job->bodies[b].nodes.contact_mx_t[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz] =
                        mx_t * m2 / m;
                job->bodies[b].nodes.contact_my_t[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz] =
                        my_t * m2 / m;
                job->bodies[b].nodes.contact_mz_t[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz] =
                        mz_t * m2 / m;

                //job->bodies[b].nodes.m[nx + job->Nx * job->Ny * nz ] = m;
                //job->bodies[b].nodes.m[nx + job->Nx * (job->Ny-1) + job->Nx * job->Ny * nz] = m;

                job->bodies[b].nodes.contact_fx[nx + job->Nx * (job->Ny + 1) * nz] = fx * m1 / m;
                job->bodies[b].nodes.contact_fy[nx + job->Nx * (job->Ny + 1) * nz] = fy * m1 / m;
                job->bodies[b].nodes.contact_fz[nx + job->Nx * (job->Ny + 1) * nz] = fz * m1 / m;
                job->bodies[b].nodes.contact_fx[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz] = fx * m2 / m;
                job->bodies[b].nodes.contact_fy[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz] = fy * m2 / m;
                job->bodies[b].nodes.contact_fz[nx + job->Nx * (job->Ny - 1) + job->Nx * job->Ny * nz] = fz * m2 / m;

                //wrap particle positions
                for (size_t i = 0; i < job->bodies[b].p; i++) {
                    if (job->bodies[b].particles.y[i] > (job->Ly - job->hy / 2.0)) {
                        job->bodies[b].particles.y[i] = job->bodies[b].particles.y[i] - (job->Ly - job->hy);
                    } else if (job->bodies[b].particles.y[i] < (job->hy / 2.0)) {
                        job->bodies[b].particles.y[i] = job->bodies[b].particles.y[i] + (job->Ly - job->hy);
                    }
                }
            }
        }
    }

    // set x-wall value (0 + Nx*n),(Nx-1 + Nx*n)
    for (size_t ny = 0; ny < job->Ny; ny++) {
        for (size_t nz = 0; nz<job->Nz;nz++) {
            /*job->u_dirichlet[NODAL_DOF * (job->Nx * ny) + NODAL_DOF * (job->Nx * job->Ny * nz) + XDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny) + NODAL_DOF * (job->Nx * job->Ny * nz) + XDOF_IDX] = 1;

            job->u_dirichlet[NODAL_DOF * (job->Nx * (ny + 1) - 1) + NODAL_DOF * (job->Nx * job->Ny * nz) + XDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * (ny + 1) - 1) + NODAL_DOF * (job->Nx * job->Ny * nz) + XDOF_IDX] = 1;*/
            for (size_t b = 0; b < job->num_bodies; b++) {
                mx_t = 0;
                my_t = 0;
                mz_t = 0;
                m = 0;
                fx = 0;
                fy = 0;
                fz = 0;

                mx_t += job->bodies[b].nodes.contact_mx_t[job->Nx * ny + job->Nx * job->Ny * nz];
                my_t += job->bodies[b].nodes.contact_my_t[job->Nx * ny + job->Nx * job->Ny * nz];
                mz_t += job->bodies[b].nodes.contact_mz_t[job->Nx * ny + job->Nx * job->Ny * nz];
                mx_t += job->bodies[b].nodes.contact_mx_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2];
                my_t += job->bodies[b].nodes.contact_my_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2];
                mz_t += job->bodies[b].nodes.contact_mz_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2];

                fx += job->bodies[b].nodes.contact_fx[job->Nx * ny + job->Nx * job->Ny * nz];
                fy += job->bodies[b].nodes.contact_fy[job->Nx * ny + job->Nx * job->Ny * nz];
                fz += job->bodies[b].nodes.contact_fz[job->Nx * ny + job->Nx * job->Ny * nz];
                fx += job->bodies[b].nodes.contact_fx[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2];
                fy += job->bodies[b].nodes.contact_fy[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2];
                fz += job->bodies[b].nodes.contact_fz[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2];

                m1 = job->bodies[b].nodes.m[job->Nx * ny + job->Nx * job->Ny * nz];
                m2 = job->bodies[b].nodes.m[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2];
                m = m1 + m2;

                job->bodies[b].nodes.contact_mx_t[job->Nx * ny + job->Nx * job->Ny * nz] = mx_t * m1 / m;
                job->bodies[b].nodes.contact_my_t[job->Nx * ny + job->Nx * job->Ny * nz] = my_t * m1 / m;
                job->bodies[b].nodes.contact_mz_t[job->Nx * ny + job->Nx * job->Ny * nz] = mz_t * m1 / m;
                job->bodies[b].nodes.contact_mx_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2] = mx_t * m2 / m;
                job->bodies[b].nodes.contact_my_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2] = my_t * m2 / m;
                job->bodies[b].nodes.contact_mz_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2] = mz_t * m2 / m;

                //job->bodies[b].nodes.m[job->Nx * ny + job->Nx * job->Ny *nz] = m;
                //job->bodies[b].nodes.m[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2] = m;

                job->bodies[b].nodes.contact_fx[job->Nx * ny + job->Nx * job->Ny * nz] = fx * m1 / m;
                job->bodies[b].nodes.contact_fy[job->Nx * ny + job->Nx * job->Ny * nz] = fy * m1 / m;
                job->bodies[b].nodes.contact_fz[job->Nx * ny + job->Nx * job->Ny * nz] = fz * m1 / m;
                job->bodies[b].nodes.contact_fx[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2] = fx * m2 / m;
                job->bodies[b].nodes.contact_fy[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2] = fy * m2 / m;
                job->bodies[b].nodes.contact_fz[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 2] = fz * m2 / m;

                //set lhs+1 to rhs
                mx_t = 0;
                my_t = 0;
                mz_t = 0;
                m = 0;
                fx = 0;
                fy = 0;
                fz = 0;

                mx_t += job->bodies[b].nodes.contact_mx_t[job->Nx * ny + job->Nx * job->Ny * nz + 1];
                my_t += job->bodies[b].nodes.contact_my_t[job->Nx * ny + job->Nx * job->Ny * nz + 1];
                mz_t += job->bodies[b].nodes.contact_mz_t[job->Nx * ny + job->Nx * job->Ny * nz + 1];
                mx_t += job->bodies[b].nodes.contact_mx_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1];
                my_t += job->bodies[b].nodes.contact_my_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1];
                mz_t += job->bodies[b].nodes.contact_mz_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1];

                fx += job->bodies[b].nodes.contact_fx[job->Nx * ny + job->Nx * job->Ny * nz + 1];
                fy += job->bodies[b].nodes.contact_fy[job->Nx * ny + job->Nx * job->Ny * nz + 1];
                fz += job->bodies[b].nodes.contact_fz[job->Nx * ny + job->Nx * job->Ny * nz + 1];
                fx += job->bodies[b].nodes.contact_fx[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1];
                fy += job->bodies[b].nodes.contact_fy[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1];
                fz += job->bodies[b].nodes.contact_fz[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1];

                m1 = job->bodies[b].nodes.m[job->Nx * ny + job->Nx * job->Ny * nz + 1];
                m2 = job->bodies[b].nodes.m[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1];
                m = m1 + m2;

                job->bodies[b].nodes.contact_mx_t[job->Nx * ny + job->Nx * job->Ny * nz + 1] = mx_t * m1 / m;
                job->bodies[b].nodes.contact_my_t[job->Nx * ny + job->Nx * job->Ny * nz + 1] = my_t * m1 / m;
                job->bodies[b].nodes.contact_mz_t[job->Nx * ny + job->Nx * job->Ny * nz + 1] = mz_t * m1 / m;
                job->bodies[b].nodes.contact_mx_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1] = mx_t * m2 / m;
                job->bodies[b].nodes.contact_my_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1] = my_t * m2 / m;
                job->bodies[b].nodes.contact_mz_t[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1] = mz_t * m2 / m;

                //job->bodies[b].nodes.m[job->Nx * ny + job->Nx * job->Ny *nz + 1] = m;
                //job->bodies[b].nodes.m[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1] = m;

                job->bodies[b].nodes.contact_fx[job->Nx * ny + job->Nx * job->Ny * nz + 1] = fx * m1 / m;
                job->bodies[b].nodes.contact_fy[job->Nx * ny + job->Nx * job->Ny * nz + 1] = fy * m1 / m;
                job->bodies[b].nodes.contact_fz[job->Nx * ny + job->Nx * job->Ny * nz + 1] = fz * m1 / m;
                job->bodies[b].nodes.contact_fx[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1] = fx * m2 / m;
                job->bodies[b].nodes.contact_fy[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1] = fy * m2 / m;
                job->bodies[b].nodes.contact_fz[job->Nx * (ny + 1) + job->Nx * job->Ny * nz - 1] = fz * m2 / m;

                //wrap particle positions
                for (size_t i = 0; i < job->bodies[b].p; i++) {
                    if (job->bodies[b].particles.x[i] > (job->Lx - job->hx / 2.0)) {
                        job->bodies[b].particles.x[i] = job->bodies[b].particles.x[i] - (job->Lx - job->hx);
                    } else if (job->bodies[b].particles.x[i] < (job->hx / 2.0)) {
                        job->bodies[b].particles.x[i] = job->bodies[b].particles.x[i] + (job->Lx - job->hx);
                    }
                }
            }
        }
    }

    return;
}

void generate_node_number_override(job_t *job)
{
    // set node_number_overide to node number (i.e. not periodic)
    for (size_t i = 0; i < job->num_nodes; i++) {
        for (size_t j = 0; j < NODAL_DOF; j++) {
            job->node_number_override[NODAL_DOF * i + j] = (NODAL_DOF * i + j);
        }
    }
    return;
}

/* Only zeros out entries for now... */
extern "C" void bc_momentum(job_t *job)
{
    for (size_t i = 0; i < job->num_nodes; i++) {
        for (size_t j = 0; j < NODAL_DOF; j++) {
            size_t n = NODAL_DOF * i + j;
            size_t m = job->node_number_override[n];
            if (job->u_dirichlet_mask[m] != 0) {
                /* only handle 0 displacement right now. */
                if (job->u_dirichlet[m] == 0) {
                    for (size_t b = 0; b < job->num_bodies; b++) {
                        if (j == XDOF_IDX) {
                            job->bodies[b].nodes.contact_mx_t[i] = 0;
                        } else if (j == YDOF_IDX) {
                            job->bodies[b].nodes.contact_my_t[i] = 0;
                        } else if (j == ZDOF_IDX) {
                            job->bodies[b].nodes.contact_mz_t[i] = 0;
                        }
                    }
                }
            }
        }
    }
    return;
}

/* Only zero out entries for now... */
extern "C" void bc_force(job_t *job)
{
    for (size_t i = 0; i < job->num_nodes; i++) {
        for (size_t j = 0; j < NODAL_DOF; j++) {
            size_t n = NODAL_DOF * i + j;
            size_t m = job->node_number_override[n];
            if (job->u_dirichlet_mask[m] != 0) {
                /* only handle 0 displacement right now. */
                if (job->u_dirichlet[m] == 0) {
                    for (size_t b = 0; b < job->num_bodies; b++) {
                        if (j == XDOF_IDX) {
                            job->bodies[b].nodes.contact_fx[i] = 0;
                        } else if (j == YDOF_IDX) {
                            job->bodies[b].nodes.contact_fy[i] = 0;
                        } else if (j == ZDOF_IDX) {
                            job->bodies[b].nodes.contact_fz[i] = 0;
                        }
                    }
                }
            }
        }
    }
    return;
}


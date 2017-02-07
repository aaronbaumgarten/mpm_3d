//
// Created by aaron on 1/31/17.
// syringeBC2d.cpp
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

double v = 0;

extern "C" void generate_dirichlet_bcs(job_t *job);
extern "C" void generate_node_number_override(job_t *job);

extern "C" void bc_init(job_t *job)
{
    if (job->boundary.num_fp64_props < 1) {
        std::cout << job->boundary.num_fp64_props << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 1 property defined (wall velocity).\n",
                __FILE__, __func__);
        exit(EXIT_FAILURE);
    }

    if (job->Nx < 12 || job->Ny < 70) {
        std::cout << job->Nx << "x" << job->Ny << std::endl;
        fprintf(stderr,
                "%s:%s: Need at 12x70 grid.\n",
                __FILE__, __func__);
        exit(EXIT_FAILURE);
    }

    v = job->boundary.fp64_props[0];
    printf("Boundary properties (v = %g).\n", v);
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

    // set floor value (0 to Nx and Nx*Ny to Nx*Ny + Nx)
    for (size_t nx = 0; nx < job->Nx; nx++) {
        job->u_dirichlet[NODAL_DOF * nx + YDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * nx + YDOF_IDX] = 1;

        job->u_dirichlet[NODAL_DOF * (nx + job->Nx*job->Ny) + YDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (nx + job->Nx*job->Ny) + YDOF_IDX] = 1;


        job->u_dirichlet[NODAL_DOF * nx + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * nx + XDOF_IDX] = 1;

        job->u_dirichlet[NODAL_DOF * (nx + job->Nx*job->Ny) + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (nx + job->Nx*job->Ny) + XDOF_IDX] = 1;
    }

    // set wall value (0 + Nx*n),(Nx-1 + Nx*n)
    for (size_t ny = 0; ny < job->Ny*job->Nz; ny++) {
        job->u_dirichlet[NODAL_DOF * (job->Nx * ny) + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny) + XDOF_IDX] = 1;

        job->u_dirichlet[NODAL_DOF * (job->Nx * (ny+1) - 1) + XDOF_IDX] = 0;
        job->u_dirichlet_mask[NODAL_DOF * (job->Nx * (ny+1) - 1) + XDOF_IDX] = 1;
    }

    //set syringe value
    //syringe width: 10 elements
    //syringe height: 70 elements
    //lock x-direction, 1 is -y velocity
    for (size_t ny = 0; ny < job->Ny*job->Nz; ny++){
        if ((job->Ny - ny%job->Ny) <= 51) {
            //upper syringe
            //left syringe wall
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny) + YDOF_IDX] = 1;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny) + YDOF_IDX] = 1;

            //right syringe wall
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny + 10) + XDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny + 10) + XDOF_IDX] = 1;
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny + 10) + YDOF_IDX] = 1;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny + 10) + YDOF_IDX] = 1;

            //second right syringe wall
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny + 11) + XDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny + 11) + XDOF_IDX] = 1;
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny + 11) + YDOF_IDX] = 1;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny + 11) + YDOF_IDX] = 1;
        } else if ((job->Ny - ny%job->Ny) <= 61){
            //buffer area
            //right syringe wall
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny + 10) + XDOF_IDX] = 1;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny + 10) + XDOF_IDX] = 1;

            //second right syringe wall
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny + 11) + XDOF_IDX] = 1;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny + 11) + XDOF_IDX] = 1;
        } else if ((job->Ny - ny%job->Ny) <= 71){
            //mouth
            //right syringe wall
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny + 10) + XDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny + 10) + XDOF_IDX] = 1;
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny + 10) + YDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny + 10) + YDOF_IDX] = 1;

            //second right syringe wall
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny + 11) + XDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny + 11) + XDOF_IDX] = 1;
            job->u_dirichlet[NODAL_DOF * (job->Nx * ny + 11) + YDOF_IDX] = 0;
            job->u_dirichlet_mask[NODAL_DOF * (job->Nx * ny + 11) + YDOF_IDX] = 1;
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
                } else if (job->u_dirichlet[m] == 1) {
                    for (size_t b = 0; b < job->num_bodies; b++) {
                        if (j == XDOF_IDX && job->bodies[b].nodes.contact_mx_t[i] > 0) {
                            job->bodies[b].nodes.contact_mx_t[i] = 0;
                        } else if (j == YDOF_IDX) {
                            job->bodies[b].nodes.contact_my_t[i] = job->bodies[b].nodes.m[i] * v;
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
                } else if (job->u_dirichlet[m] == 1) {
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


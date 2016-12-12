//
// Created by aaron on 11/12/16.
// borrows heavily from mpm_2d_legacy
// floorBC.cpp
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

void generate_dirichlet_bcs(job_t *job);
void generate_node_number_override(job_t *job);

void bc_init(job_t *job)
{
    // suppress unused argument warning
    (void)job;
    return;
}

int bc_validate(job_t *job)
{
    //return 0 for failure to validate properties in configuration file
    //not implemented as of 10/26/16
    return 1;
}

void bc_time_varying(job_t *job)
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
void bc_momentum(job_t *job)
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
                            job->bodies[b].node_contact_mx_t[i] = 0;
                        } else if (j == YDOF_IDX) {
                            job->bodies[b].node_contact_my_t[i] = 0;
                        } else if (j == ZDOF_IDX) {
                            job->bodies[b].node_contact_mz_t[i] = 0;
                        }
                    }
                }
            }
        }
    }
    return;
}

/* Only zero out entries for now... */
void bc_force(job_t *job)
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
                            job->bodies[b].node_contact_fx[i] = 0;
                        } else if (j == YDOF_IDX) {
                            job->bodies[b].node_contact_fy[i] = 0;
                        } else if (j == ZDOF_IDX) {
                            job->bodies[b].node_contact_fz[i] = 0;
                        }
                    }
                }
            }
        }
    }
    return;
}

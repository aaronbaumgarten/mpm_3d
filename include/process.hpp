//
// Created by aaron on 8/26/16.
// process.hpp
//

#ifndef MPM_3D_PROCESS_HPP
#define MPM_3D_PROCESS_HPP

#include <stdlib.h>
#include <vector>
#include "body.hpp"

class job_t{
public:
    //time tracking
    double t;
    double dt;
    size_t frame;
    int stepcount;

    //total count
    size_t num_bodies;
    size_t num_particles;
    size_t num_nodes;
    size_t num_elements;

    //thread coloring
    /*size_t num_colors;
    size_t *color_indices;
    size_t *particle_by_element_color_lengths;
    size_t **particle_by_element_color_lists;*/

    //grid
    size_t Nx, Ny, Nz;
    double Lx, Ly, Lz;
    double hx, hy, hz;

    //mpm type
    int use_cpdi;

    //node positions
    int *node_u_map;
    int *inv_node_u_map;

    //which DOF are bc controlled
    double *u_dirichlet;
    double *u_dirichlet_mask;

    //step number
    int step_number;
    double step_start_time;

    //object vector
    std::vector<Body> bodies;

    //threading
    //pthread_barrier_t *step_barrier;
    //pthread_barrier_t *serialize_barrier;
    //size_t num_threads;

    //functions
    job_t();
    int importNodesandParticles(const char*,const char*);
    void createBody(Body*,size_t,size_t,size_t);
};

class threadtask_t{
public:
    size_t id;
    static size_t num_threads;

    //responsibility offset to offset+blocksize
    size_t offset;
    size_t blocksize;

    //node responsibility
    size_t n_offset;
    size_t n_blocksize;

    //element responsibilty
    size_t e_offset;
    size_t e_blocksize;

    size_t stride;

    job_t *job;
};

#endif //MPM_3D_PROCESS_HPP

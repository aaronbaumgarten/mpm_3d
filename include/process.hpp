//
// Created by aaron on 8/26/16.
// process.hpp
//

#ifndef MPM_3D_PROCESS_HPP
#define MPM_3D_PROCESS_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include "body.hpp"
#include "element.hpp"
#include "boundary.hpp"

class job_t{
public:
    //time tracking
    double t;
    double dt;
    size_t frame;
    int stepcount;

    //implicit
    double linearStepSize;
    double newtonTOL;
    double dt_base;
    double dt_minimum;

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
    int use_3d;
    int use_implicit;
    int use_cpdi;

    //which DOF are bc controlled
    std::vector<double> u_dirichlet;
    std::vector<double> u_dirichlet_mask;
    std::vector<int> node_number_override;

    //step number
    int step_number;
    double step_start_time;

    //object vector
    std::vector<Body> bodies;

    //element object
    Elements elements;

    //boundary condition
    Boundary boundary;

    //threading
    //pthread_barrier_t *step_barrier;
    //pthread_barrier_t *serialize_barrier;
    //size_t num_threads;

    //functions
    job_t();
    inline void node_number_to_coords(double *, double *, double *, size_t, size_t, size_t, size_t, double, double, double);
    inline int ijkton_safe(int,int,int,int,int,int);

    //initialization
    //int importNodesandParticles(const char*,const char*);
    int importNodesandParticles(std::string ,std::string);
    int importNodesandParticles2D(std::string, std::string);
    void createBody(Body*,size_t,size_t,size_t,size_t);
    int assignMaterials();
    int assignMaterials(const char*, const char*);
    int assignBoundaryConditions();
    int assignBoundaryConditions(const char*);

    int mpmStepUSLExplicit();
    int mpmStepUSLImplicit();
    int mpmStepUSLExplicit2D();
    int mpmTrialStepUSLExplicit();
    //usl step
    int createMappings();
    void mapParticles2Grid();
    void addContactForces();
    void addContactForces2D();
    void addBoundaryConditions();
    void moveGridExplicit();
    void moveGridExplicit2D();
    void moveParticlesExplicit();
    void moveParticlesImplicit();
    void moveParticlesExplicit2D();
    void calculateStrainRate();
    void calculateStrainRate2D();
    void updateDensity();
    void updateTrialDensity();
    void updateStress();
    void updateTrialStress();
    //implicit calculations
    void mapTrialStress2Grid();
    void calculateImplicitResidual();
    void moveGridImplicitCG();
    void moveGridImplicitBiCGSTAB();

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

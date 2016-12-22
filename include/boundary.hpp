//
// Created by aaron on 10/26/16.
// boundary.hpp
//

#include <stdlib.h>
#include <functional>

#ifndef MPM_3D_BOUNDARY_HPP
#define MPM_3D_BOUNDARY_HPP

class Body;
class job_t;
class threadtask_t;

class Boundary{
public:
    std::string boundary_filename;
    int use_builtin;

    void *handle;

    void (*bc_init)(job_t*);
    void (*bc_validate)(job_t*);
    void (*bc_time_varying)(job_t*);

    void (*generate_dirichlet_bcs)(job_t*);
    void (*generate_node_number_override)(job_t*);
    void (*bc_momentum)(job_t*);
    void (*bc_force)(job_t*);


    double *fp64_props;
    int *int_props;
    size_t num_fp64_props;
    size_t num_int_props;

    Boundary();
    Boundary(std::string,size_t,size_t,double*,int*);
    ~Boundary();
    void setBoundary(std::string,size_t,size_t,double*,int*);
};
#endif //MPM_3D_BOUNDARY_HPP

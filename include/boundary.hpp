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
    char *boundary_filename;
    int use_builtin;

    std::function<void(job_t *)> bc_init;
    std::function<void(job_t *)> bc_validate;
    std::function<void(job_t *)> bc_time_varying;

    std::function<void(job_t *)> generate_dirichlet_bcs;
    std::function<void(job_t *)> generate_node_number_override;
    std::function<void(job_t *)> bc_momentum;
    std::function<void(job_t *)> bc_force;


    double *fp64_props;
    int *int_props;
    size_t num_fp64_props;
    size_t num_int_props;

    Boundary(){
        num_fp64_props = 0;
        num_int_props = 0;
        fp64_props = NULL;
        int_props = NULL;
        boundary_filename = NULL;

        /*material_init = NULL;
        calculate_stress = NULL;
        calculate_stress_threaded = NULL;*/
    }
};

namespace boundary {
    void bc_init(job_t *job);
    int bc_validate(job_t *job);
    void bc_time_varying(job_t *job);
    void generate_dirichlet_bcs(job_t *job);
    void generate_node_number_override(job_t *job);
    void bc_momentum(job_t *job);
    void bc_force(job_t *job);
}
#endif //MPM_3D_BOUNDARY_HPP

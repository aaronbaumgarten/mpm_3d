//
// Created by aaron on 9/6/16.
// material.hpp
//

#include <stdlib.h>
#include <functional>

#ifndef MPM_3D_MATERIAL_HPP
#define MPM_3D_MATERIAL_HPP

#endif //MPM_3D_MATERIAL_HPP

class Body;
class threadtask_t;

class Material{
public:
    char *material_filename;
    int use_builtin;

    //template <class bodyT, class taskT>
    /*void (*material_init)(Body *);
    void (*calculate_stress)(Body *, double);
    void (*calculate_stress_threaded)(threadtask_t *, Body *);*/

    std::function<void(Body *)> material_init;
    std::function<void(Body *, double)> calculate_stress;
    std::function<void(threadtask_t *, Body *)> calculate_stress_threaded;

    double *fp64_props;
    int *int_props;
    size_t num_fp64_props;
    size_t num_int_props;

    Material(){
        num_fp64_props = 0;
        num_int_props = 0;
        fp64_props = NULL;
        int_props = NULL;
        material_filename = NULL;

        /*material_init = NULL;
        calculate_stress = NULL;
        calculate_stress_threaded = NULL;*/
    }
};

namespace material1 {
    //template <class bodyT>
    //template <class taskT>
    void material_init(Body *);
    void calculate_stress(Body *, double);
    void calculate_stress_threaded(threadtask_t *, Body *);
}

namespace material2 {
    //template <class bodyT>
    //template <class taskT>
    void material_init(Body *);
    void calculate_stress(Body *, double);
    void calculate_stress_threaded(threadtask_t *, Body *);
}
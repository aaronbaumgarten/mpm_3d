//
// Created by aaron on 9/6/16.
// material.hpp
//

#include <stdlib.h>
#include <string>

#ifndef MPM_3D_MATERIAL_HPP
#define MPM_3D_MATERIAL_HPP

class Body;
class threadtask_t;

class Material{
public:
    std::string material_filename;
    int use_builtin;
    void *handle;

    void (*material_init)(Body*);
    void (*calculate_stress)(Body*, double);
    void (*calculate_stress_implicit)(Body*, double);
    void (*calculate_stress_threaded)(threadtask_t*,Body*,double);

    double *fp64_props;
    int *int_props;
    size_t num_fp64_props;
    size_t num_int_props;

    Material();
    Material(std::string,size_t,size_t,double*,int*);
    ~Material();
    void setMaterial(std::string,size_t,size_t,double*,int*);
};

/*namespace material1 {
    //template <class bodyT>
    //template <class taskT>
    void material_init(Body *);
    void calculate_stress(Body *, double);
    void calculate_stress_threaded(threadtask_t *, Body *, double);
    void calculate_stress_implicit(Body *, double);
}

namespace material2 {
    //template <class bodyT>
    //template <class taskT>
    void material_init(Body *);
    void calculate_stress(Body *, double);
    void calculate_stress_threaded(threadtask_t *, Body *, double);
    void calculate_stress_implicit(Body *, double);
}*/

#endif //MPM_3D_MATERIAL_HPP
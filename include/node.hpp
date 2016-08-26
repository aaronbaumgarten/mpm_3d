//
// Created by aaron on 8/26/16.
// node.hpp (heavily borrow from mpm-2d-legacy)
//

#ifndef MPM_3D_NODE_HPP
#define MPM_3D_NODE_HPP

inline size_t get_global_index(size_t body, size_t local_index, size_t blocksize){
    return (body * blocksize + local_index);
}

class Node{
public:
    //filled element neighbors
    int num_filled_element_neighbors;
    double mass_filled_element_neighbors;

    //mass
    double m;
    double body_m;

    //position
    double x;
    double y;
    double z;

    //displacement
    double ux;
    double uy;
    double uz;

    //velocity
    double x_t;
    double y_t;
    double z_t;

    //velocity difference
    double diff_x_t;
    double diff_y_t;
    double diff_z_t;

    //momentum
    double mx_t;
    double my_t;
    double mz_t;

    //force
    double fx;
    double fy;
    double fz;

    //density
    double rho;

    //marker for displacement/velocity update
    double velocity_update_flag;

    double sum_sqrt_m_neighbors;
    double max_m_neighbors;

    //body contact resolution
    double contact_x_t;
    double contact_y_t;
    double contact_z_t;

    double contact_fx;
    double contact_fy;
    double contact_fz;

    double real_contact_fx;
    double real_contact_fy;
    double real_contact_fz;

    double contact_normal_x;
    double contact_normal_y;
    double contact_normal_z;

    size_t num_bodies;

    size_t blocksize;

    //construcors
    Node() {}
    //destructors
    ~Node() {}
};

#endif //MPM_3D_NODE_HPP
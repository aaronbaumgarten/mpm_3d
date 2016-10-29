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
    //unique id
    size_t id;

    //filled element neighbors
    int num_filled_element_neighbors;
    double mass_filled_element_neighbors;

    //mass
    double* m;
    double* body_m;

    //position
    double* x;
    double* y;
    double* z;

    //displacement
    double* ux;
    double* uy;
    double* uz;

    //velocity
    double* x_t;
    double* y_t;
    double* z_t;

    //velocity difference
    double* diff_x_t;
    double* diff_y_t;
    double* diff_z_t;

    //momentum
    double* mx_t;
    double* my_t;
    double* mz_t;

    //force
    double* fx;
    double* fy;
    double* fz;

    //density
    double* rho;

    //marker for displacement/velocity update
    double velocity_update_flag;

    double sum_sqrt_m_neighbors;
    double max_m_neighbors;

    //body contact resolution
    double* contact_x_t;
    double* contact_y_t;
    double* contact_z_t;

    double* contact_mx_t;
    double* contact_my_t;
    double* contact_mz_t;

    double* contact_fx;
    double* contact_fy;
    double* contact_fz;

    double* real_contact_fx;
    double* real_contact_fy;
    double* real_contact_fz;

    double* contact_normal_x;
    double* contact_normal_y;
    double* contact_normal_z;

    size_t num_bodies;

    size_t blocksize;

    //construcors

    template <class bodyT>
    Node(bodyT* bd, size_t idIn):
            id(idIn),
            //mass
            m(&(bd->node_m[idIn])),

            //position
            x(&(bd->node_x[idIn])),
            y(&(bd->node_y[idIn])),
            z(&(bd->node_z[idIn])),

            //displacement
            ux(&(bd->node_ux[idIn])),
            uy(&(bd->node_uy[idIn])),
            uz(&(bd->node_uz[idIn])),

            //velocity
            x_t(&(bd->node_x_t[idIn])),
            y_t(&(bd->node_y_t[idIn])),
            z_t(&(bd->node_z_t[idIn])),

            //velocity difference
            diff_x_t(&(bd->node_diff_x_t[idIn])),
            diff_y_t(&(bd->node_diff_y_t[idIn])),
            diff_z_t(&(bd->node_diff_z_t[idIn])),

            //momentum
            mx_t(&(bd->node_mx_t[idIn])),
            my_t(&(bd->node_my_t[idIn])),
            mz_t(&(bd->node_mz_t[idIn])),

            //force
            fx(&(bd->node_fx[idIn])),
            fy(&(bd->node_fy[idIn])),
            fz(&(bd->node_fz[idIn])),

            //density
            rho(&(bd->node_rho[idIn])),

            //body contact resolution
            contact_mx_t(&(bd->node_contact_mx_t[idIn])),
            contact_my_t(&(bd->node_contact_my_t[idIn])),
            contact_mz_t(&(bd->node_contact_mz_t[idIn])),

            contact_x_t(&(bd->node_contact_x_t[idIn])),
            contact_y_t(&(bd->node_contact_y_t[idIn])),
            contact_z_t(&(bd->node_contact_z_t[idIn])),

            contact_fx(&(bd->node_contact_fx[idIn])),
            contact_fy(&(bd->node_contact_fy[idIn])),
            contact_fz(&(bd->node_contact_fz[idIn])),

            real_contact_fx(&(bd->node_real_contact_fx[idIn])),
            real_contact_fy(&(bd->node_real_contact_fy[idIn])),
            real_contact_fz(&(bd->node_real_contact_fz[idIn])),

            contact_normal_x(&(bd->node_contact_normal_x[idIn])),
            contact_normal_y(&(bd->node_contact_normal_y[idIn])),
            contact_normal_z(&(bd->node_contact_normal_z[idIn]))
    { }

    Node() {}
    //destructors
    //~Node() {}
};

#endif //MPM_3D_NODE_HPP
//
// Created by aaron on 8/26/16.
// node.hpp (heavily borrow from mpm-2d-legacy)
//

#include <eigen3/Eigen/Core>

#ifndef MPM_3D_NODE_HPP
#define MPM_3D_NODE_HPP

class Nodes{
public:
    //unique id
    size_t num_nodes;

    //mass
    Eigen::VectorXd m;

    //position
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd z;

    //displacement
    Eigen::VectorXd ux;
    Eigen::VectorXd uy;
    Eigen::VectorXd uz;

    //velocity
    Eigen::VectorXd x_t;
    Eigen::VectorXd y_t;
    Eigen::VectorXd z_t;

    //velocity difference
    Eigen::VectorXd diff_x_t;
    Eigen::VectorXd diff_y_t;
    Eigen::VectorXd diff_z_t;

    //momentum
    Eigen::VectorXd mx_t;
    Eigen::VectorXd my_t;
    Eigen::VectorXd mz_t;

    //force
    Eigen::VectorXd fx;
    Eigen::VectorXd fy;
    Eigen::VectorXd fz;

    //density
    Eigen::VectorXd rho;

    //body contact resolution
    Eigen::VectorXd contact_x_t;
    Eigen::VectorXd contact_y_t;
    Eigen::VectorXd contact_z_t;

    Eigen::VectorXd contact_mx_t;
    Eigen::VectorXd contact_my_t;
    Eigen::VectorXd contact_mz_t;

    Eigen::VectorXd contact_fx;
    Eigen::VectorXd contact_fy;
    Eigen::VectorXd contact_fz;

    Eigen::VectorXd contact_normal_x;
    Eigen::VectorXd contact_normal_y;
    Eigen::VectorXd contact_normal_z;
    
    //nodal initial momentum
    Eigen::VectorXd mx_t_k;
    Eigen::VectorXd my_t_k;
    Eigen::VectorXd mz_t_k;

    //nodal trial velocity
    Eigen::VectorXd x_t_trial;
    Eigen::VectorXd y_t_trial;
    Eigen::VectorXd z_t_trial;

    //nodal iterated velocity
    Eigen::VectorXd x_t_n;
    Eigen::VectorXd y_t_n;
    Eigen::VectorXd z_t_n;

    //nodal initial force
    //Eigen::VectorXd fx_k;
    //Eigen::VectorXd fy_k;
    //Eigen::VectorXd fz_k;

    //nodal final force
    Eigen::VectorXd fx_L;
    Eigen::VectorXd fy_L;
    Eigen::VectorXd fz_L;

    //nodal residuals
    Eigen::VectorXd Rx;
    Eigen::VectorXd Ry;
    Eigen::VectorXd Rz;

    //saved residuals
    Eigen::VectorXd Rvx;
    Eigen::VectorXd Rvy;
    Eigen::VectorXd Rvz;

    //directional residual derivative
    Eigen::VectorXd DhRx;
    Eigen::VectorXd DhRy;
    Eigen::VectorXd DhRz;

    //implicit algorithm
    Eigen::VectorXd wk;
    double ak;
    Eigen::VectorXd sk;
    Eigen::VectorXd rk;
    Eigen::VectorXd r0;
    Eigen::VectorXd qk;
    Eigen::VectorXd tk;
    Eigen::VectorXd hk;
    double ok;
    double rhok;
    double bk;
    Eigen::VectorXd pk;

    //minimum residuals
    Eigen::VectorXd skMin;
    Eigen::VectorXd rkMin;

    //construcors
    Nodes(size_t);
    Nodes() {}
    //destructors
    //~Node() {}

    //functions
    void addNode(double,double,double,size_t);
};

#endif //MPM_3D_NODE_HPP
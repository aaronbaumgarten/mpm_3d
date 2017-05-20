//
// Created by aaron on 5/5/17.
// nodes.hpp
//

#ifndef MPM_3D_NODES_HPP
#define MPM_3D_NODES_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

class Job;
class Serializer;
class Body;

class Nodes{
public:
    //nodes properties here
    Eigen::MatrixXd x;
    Eigen::MatrixXd u;
    Eigen::MatrixXd x_t;
    Eigen::MatrixXd diff_x_t;
    Eigen::VectorXd m;
    Eigen::MatrixXd mx_t;
    Eigen::MatrixXd f;
    Eigen::VectorXi active;

    //objects here

    //nodes object specific functions
    Nodes();
    int nodesInit(Job*, Body*);

    void nodesWriteFrame(Job*, Body*, Serializer*); //write frame data to serializer
    std::string nodesSaveState(Job*, Body*, Serializer*,std::string); //save data to file in serializer directory and return name
    int nodesLoadState(Job*, Body*, Serializer*,std::string); //load data from full path

    //other functions
};

#endif //MPM_3D_NODES_HPP

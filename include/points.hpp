//
// Created by aaron on 5/5/17.
// points.hpp
//

#ifndef MPM_3D_POINTS_HPP
#define MPM_3D_POINTS_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

class Job;
class Serializer;
class Body;

class Points{
public:
    //point-file
    std::string file;

    //points properties here
    Eigen::MatrixXd x;
    Eigen::MatrixXd u;
    Eigen::MatrixXd x_t;
    Eigen::VectorXd m;
    Eigen::VectorXd v;
    Eigen::MatrixXd mx_t;
    Eigen::MatrixXd b;
    Eigen::MatrixXd T;
    Eigen::MatrixXd L;
    Eigen::VectorXi active;
    Eigen::VectorXd extent;

    //objects here

    //points object specific functions
    Points();
    int pointsInit(Job*, Body*);
    int pointsReadFromFile(Job*, Body*, std::string);

    void pointsWriteFrame(Job*, Body*, Serializer*); //write frame data to serializer
    std::string pointsSaveState(Job*, Body*, Serializer*,std::string); //save data to file in serializer directory and return name
    int pointsLoadState(Job*, Body*, Serializer*,std::string); //load data from full path

    //other functions
};

#endif //MPM_3D_POINTS_HPP

//
// Created by aaron on 5/6/17.
// grid.hpp
//

#ifndef MPM_3D_GRID_HPP
#define MPM_3D_GRID_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>

class Job;
class Serializer;

class Grid{
public:
    //grid properties here
    std::string filename; //name of file
    std::string filepath; //directory of file for access
    std::vector<double> fp64_props; //double properties
    std::vector<int> int_props; //integer properties
    std::vector<std::string> str_props; //string properties
    void *handle; //.so file handle

    size_t node_count;
    size_t element_count;

    //objects here

    //grid object specific functions
    Grid();
    ~Grid();
    void gridSetPlugin(Job*, std::string, std::string, std::vector<double>, std::vector<int>, std::vector<std::string>); //assign .so plugin for functions
    void gridSetFnPointers(void*); //set function pointers to .so file handle

    void (*gridInit)(Job*); //initialize grid file

    void (*gridWriteFrame)(Job*, Serializer*); //write frame to serializer
    std::string (*gridSaveState)(Job*, Serializer*,std::string); //save state to serializer folder with given filename
    int (*gridLoadState)(Job*, Serializer*, std::string); //read state from given full path

    //other functions
    int (*gridWhichElement)(Job*, Eigen::VectorXd); //return element id
    Eigen::VectorXd (*gridNodeIDToPosition)(Job*, int); //return position of node
    void (*gridEvaluateShapeFnValue)(Job*, Eigen::VectorXd, std::vector<int>&, std::vector<double>&); //evaluate shape function values at given point (add ids and vals to vectors)
    void (*gridEvaluateShapeFnGradient)(Job*, Eigen::VectorXd, std::vector<int>&, std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>&); //evaluate shape function gradients at given point (add ids and vectors to vec/matrix)
};

#endif //MPM_3D_GRID_HPP

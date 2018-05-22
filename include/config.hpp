//
// Created by aaron on 5/10/18.
// config.hpp
//

#ifndef MPM_V3_CONFIG_HPP
#define MPM_V3_CONFIG_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "job.hpp"
#include "registry.hpp"
#include "mpm_objects.hpp"

class Configurator{
public:
    std::string file; //name of config file
    std::string mainpath; //directory of main function

    Registry<Serializer> serializer_registry;
    Registry<Driver> driver_registry;
    Registry<Solver> solver_registry;
    Registry<Grid> grid_registry;
    Registry<Body> body_registry;
    Registry<Points> points_registry;
    Registry<Nodes> nodes_registry;
    Registry<Contact> contact_registry;
    Registry<Material> material_registry;
    Registry<Boundary> boundary_registry;

    void init(std::string); //initialization
    void setMainPath(std::string); //set path to main
    void checkConfigFile(std::string); //check headers and report problems

    int configureJob(Job*); //configure job object;
};

#endif //MPM_V3_CONFIG_HPP

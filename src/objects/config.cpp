//
// Created by aaron on 5/14/17.
// config.cpp
//


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "job.hpp"
#include "config.hpp"

#include "serializer.hpp"
#include "driver.hpp"
#include "solver.hpp"

#include "grid.hpp"
#include "contact.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "points.hpp"

#include "material.hpp"
#include "boundary.hpp"

Config::Config() {
    filename = "";
    filepath = "";
}

int Config::configInit(std::string pathIN, std::string nameIN){
    filename = nameIN;
    filepath = pathIN;
    return 1;
}

int configConfigureJob(Job*){
//configure job object
    return -1;
}
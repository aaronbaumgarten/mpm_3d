//
// Created by aaron on 5/11/17.
// job.cpp
//

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "stringparser.hpp"

#include "job.hpp"

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

Job::Job():
        activeBodies(0),
        activeContacts(0),
        bodies(0),
        contacts(0)
{
    DIM = 3;
    XX = 0; XY = 1; XZ = 2;
    YX = 3; YY = 4; YZ = 5;
    ZX = 6; ZY = 7; ZZ = 8;

    X = 0; Y = 1; Z = 2;

    t = 0; dt = 1e-3;

    serializer = Serializer();
    driver = Driver();
    solver = Solver();

    grid = Grid();
}

int Job::jobInit(){

    std::cout << "Job properties (DIM = " << DIM << ", t = " << t << ", dt = " << dt << ")." << std::endl;
    std::cout << "Job Initialized." << std::endl;

    serializer.serializerInit(this);
    driver.driverInit(this);
    solver.solverInit(this);

    grid.gridInit(this);

    for (size_t i=0; i<bodies.size(); i++){
        bodies[i].bodyInit(this);
    }

    for (size_t i=0; i<contacts.size(); i++){
        contacts[i].contactInit(this, &(contacts[i]));
    }

    return 1;
}
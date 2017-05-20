//
// Created by aaron on 5/13/17.
// body.cpp
//


#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "job.hpp"

#include "serializer.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "points.hpp"

#include "material.hpp"
#include "boundary.hpp"

Body::Body():
{
    id = 0; //when config sets body (set to id)
    name = "default";
    activeMaterial = 0; //when config sets material (set to 1)
    activeBoundary = 0; //when config sets boundary (set to 1)

    points = Points();
    nodes = Nodes();
    material = Material();
    boundary = Boundary();
}

int Body::bodyInit(Job*){
    //body initialization stuff
    //kind of useless right now
    //use for stuff that can't be in constructor
    return 1;
}


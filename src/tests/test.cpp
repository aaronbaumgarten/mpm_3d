//
// Created by aaron on 9/12/16.
// test.cpp
//

#include <stdlib.h>
#include <math.h>

#include "test.hpp"
#include "process.hpp"

void testNodeandParticleSums(job_t *job) {
    double sumX = 0;
    for (size_t i = 0; i < job->num_nodes; i++) {
        sumX += job->bodies[0].nodes[i].x[0];
    }
    std::cout << "test: " << sumX << " =? 515150\n";
    sumX = 0;
    for (size_t i = 0; i < job->num_particles; i++) {
        sumX += job->bodies[0].particles[i].id;
    }
    std::cout << "test: " << sumX << " =? " << 7999 * 8000 / 2 << "\n";
}

void testParticleCorners(job_t *job) {
    job->bodies[0].particles[0].updateCorners(job);
    std::cout << "element number: " << job->bodies[0].particles[0].corner_elements[0] << "\n";
    std::cout << "particle x: " << job->bodies[0].particles[0].x[0] << "\n";
    for (size_t i = 0; i < 8; i++) {
        std::cout << "node x: " << job->bodies[0].nodes[job->bodies[0].elements[505050].nodeID[i]].x[0] << "\n";
    }
}
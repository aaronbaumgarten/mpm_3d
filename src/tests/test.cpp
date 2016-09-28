//
// Created by aaron on 9/12/16.
// test.cpp
//

#include <stdlib.h>
#include <math.h>

#include "test.hpp"
#include "process.hpp"

void testParticlePositions(job_t *job){
    std::cout << "particles[0].{x,y,x} : particles[1].{x,y,z}\n";
    std::cout << "{" << job->bodies[0].particles[0].x[0] << " " << job->bodies[0].particles[0].y[0] << " " << job->bodies[0].particles[0].z[0] << "} : ";
    std::cout << "{" << job->bodies[0].particles[1].x[0] << " " << job->bodies[0].particles[1].y[0] << " " << job->bodies[0].particles[1].z[0] << "}\n";
}

void testNodePositions(job_t *job){
    size_t N = cbrt(job->num_nodes);

    std::cout << "nodes[0].{x,y,x} : nodes[1].{x,y,z}\n";
    std::cout << "{" << job->bodies[0].nodes[0].x[0] << " " << job->bodies[0].nodes[0].y[0] << " " << job->bodies[0].nodes[0].z[0] << "} : ";
    std::cout << "{" << job->bodies[0].nodes[1].x[0] << " " << job->bodies[0].nodes[1].y[0] << " " << job->bodies[0].nodes[1].z[0] << "}\n";

    std::cout << "nodes[0].{x,y,x} : nodes[N].{x,y,z}\n";
    std::cout << "{" << job->bodies[0].nodes[0].x[0] << " " << job->bodies[0].nodes[0].y[0] << " " << job->bodies[0].nodes[0].z[0] << "} : ";
    std::cout << "{" << job->bodies[0].nodes[N].x[0] << " " << job->bodies[0].nodes[N].y[0] << " " << job->bodies[0].nodes[N].z[0] << "}\n";

    std::cout << "nodes[0].{x,y,x} : nodes[N*N].{x,y,z}\n";
    std::cout << "{" << job->bodies[0].nodes[0].x[0] << " " << job->bodies[0].nodes[0].y[0] << " " << job->bodies[0].nodes[0].z[0] << "} : ";
    std::cout << "{" << job->bodies[0].nodes[N*N].x[0] << " " << job->bodies[0].nodes[N*N].y[0] << " " << job->bodies[0].nodes[N*N].z[0] << "}\n";
}

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
    int elementID = job->bodies[0].particles[0].corner_elements[0];
    job->bodies[0].particles[0].updateCorners(job);
    std::cout << "element number: " << elementID << "\n";
    std::cout << "particle x: " << job->bodies[0].particles[0].x[0] << "\n";
    for (size_t i = 0; i < 8; i++) {
        std::cout << "node x: " << job->bodies[0].nodes[job->bodies[0].elements[elementID].nodeID[i]].x[0] << "\n";
    }
}

void testSipSize(job_t *job) {
    std::cout << job->bodies[0].n << ", " << job->bodies[0].p << " ?= " << job->bodies[0].Sip.rows() << ", " << job->bodies[0].Sip.cols() << "\n";
}
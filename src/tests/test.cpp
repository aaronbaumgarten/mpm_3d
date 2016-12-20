//
// Created by aaron on 9/12/16.
// test.cpp
//

#include <stdlib.h>
#include <math.h>
#include <Eigen/Core>

#include "test.hpp"
#include "process.hpp"

void testParticlePositions(job_t *job){
    std::cout << "particles[0].{x,y,x} : particles[1].{x,y,z}\n";
    std::cout << "{" << job->bodies[0].particles[0].x[0] << " " << job->bodies[0].particles[0].y[0] << " " << job->bodies[0].particles[0].z[0] << "} : ";
    std::cout << "{" << job->bodies[0].particles[1].x[0] << " " << job->bodies[0].particles[1].y[0] << " " << job->bodies[0].particles[1].z[0] << "}\n";
    return;
}

void testNodePositions(job_t *job){
    size_t N = cbrt(job->num_nodes);

    std::cout << "nodes[0].{x,y,x} : nodes[1].{x,y,z}\n";
    std::cout << "{" << job->bodies[0].nodes.x[0] << " " << job->bodies[0].nodes.y[0] << " " << job->bodies[0].nodes.z[0] << "} : ";
    std::cout << "{" << job->bodies[0].nodes.x[1] << " " << job->bodies[0].nodes.y[1] << " " << job->bodies[0].nodes.z[1] << "}\n";

    std::cout << "nodes[0].{x,y,x} : nodes[N].{x,y,z}\n";
    std::cout << "{" << job->bodies[0].nodes.x[0] << " " << job->bodies[0].nodes.y[0] << " " << job->bodies[0].nodes.z[0] << "} : ";
    std::cout << "{" << job->bodies[0].nodes.x[N] << " " << job->bodies[0].nodes.y[N] << " " << job->bodies[0].nodes.z[N] << "}\n";

    std::cout << "nodes[0].{x,y,x} : nodes[N*N].{x,y,z}\n";
    std::cout << "{" << job->bodies[0].nodes.x[0] << " " << job->bodies[0].nodes.y[0] << " " << job->bodies[0].nodes.z[0] << "} : ";
    std::cout << "{" << job->bodies[0].nodes.x[N*N] << " " << job->bodies[0].nodes.y[N*N] << " " << job->bodies[0].nodes.z[N*N] << "}\n";
    return;
}

void testNodeandParticleSums(job_t *job) {
    double sumX = 0;
    for (size_t i = 0; i < job->num_nodes; i++) {
        sumX += job->bodies[0].nodes.x[i];
    }
    std::cout << "test: " << sumX << " =? 515150\n";
    sumX = 0;
    for (size_t i = 0; i < job->num_particles; i++) {
        sumX += job->bodies[0].particles[i].id;
    }
    std::cout << "test: " << sumX << " =? " << 7999 * 8000 / 2 << "\n";
    return;
}

void testParticleCorners(job_t *job) {
    int elementID = job->bodies[0].particles[0].corner_elements[0];
    job->bodies[0].particles[0].updateCorners(job);
    std::cout << "element number: " << elementID << "\n";
    std::cout << "particle x: " << job->bodies[0].particles[0].x[0] << "\n";
    for (size_t i = 0; i < 8; i++) {
        std::cout << "node x: " << job->bodies[0].nodes.x[job->elements.nodeID(elementID,i)] << "\n";
    }
    return;
}

void testPhiSize(job_t *job) {
    std::cout << job->bodies[0].n << ", " << job->bodies[0].p << " ?= " << job->bodies[0].Phi.rows() << ", " << job->bodies[0].Phi.cols() << "\n";
    return;
}

void testElementCorners(job_t *job) {
    int elementID = job->bodies[0].particles[0].corner_elements[0];
    job->bodies[0].particles[0].updateCorners(job);
    std::cout << "element number: " << elementID << "\n";
    std::cout << "particle x,y,z: " << job->bodies[0].particles[0].corner[0][0] << ", " << job->bodies[0].particles[0].corner[0][1] << ", " << job->bodies[0].particles[0].corner[0][2] << "\n";
    for (size_t i = 0; i < 8; i++) {
        std::cout << "node x,y,z: " << job->bodies[0].nodes.x[job->elements.nodeID(elementID,i)] << " , ";
        std::cout << job->bodies[0].nodes.y[job->elements.nodeID(elementID,i)] << " , ";
        std::cout << job->bodies[0].nodes.z[job->elements.nodeID(elementID,i)] << "\n";
    }
    return;
}

void testMappingP2G(job_t *job) {
    size_t n = floor(job->bodies[0].n/2);
    std::cout << "node[515151].m: " << job->bodies[0].nodes.m[n] << "\n";
    job->mapParticles2Grid();
    std::cout << "node[515151].m: " << job->bodies[0].nodes.m[n] << "\n";
    return;
}

void testMappingGradient(job_t *job){
    Eigen::VectorXd cstOnes(job->bodies[0].n);
    Eigen::VectorXd pvec(job->bodies[0].p);

    for (size_t i=0;i<job->bodies[0].n;i++){
        cstOnes[i] = i;//i;
    }

    pvec = job->bodies[0].gradPhiX.transpose() * cstOnes;
    std::cout << "non-zero gradients: 0 =? " << pvec.nonZeros() << "\n";
    std::cout << pvec << "\n";
    return;
}
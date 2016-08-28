//
// Created by aaron on 8/27/16.
// process.cpp
//

#include <iostream>
#include <stdlib.h>

#include "process.hpp"

int job_t::readParticles(const char *filename){
    //open particle file and parse header
    FILE *fp;
    int r;

    size_t numParticles;
    size_t numParticles1;
    size_t numParticles2;

    char s[16384]; //#from mpm-2d-legacy

    fp = fopen(filename, "r");
    if (fp == NULL){
        std::cout << "Cannot parse particle file: " << filename << "\n";
        return -1;
    }

    if (NULL==fgets(s, sizeof(s)/sizeof(char), fp)){
        std::cout << "Cannot parse particle file: " << filename << "\n";
        return -1;
    }

    r = sscanf(s,"%i %i %i",&numParticles,&numParticles1,&numParticles2);
    if (r!=3){
        std::cout << "Cannot parse particle file: " << filename << "\n";
        std::cout << "Expected 3 particle counts at file header." << "\n";
        std::cout << "Got: " << r << "\n";
        return -1;
    }

    std::cout << "Header parsed!" << "\n";
    std::cout << "Number of particles: " << numParticles << "\n";
}
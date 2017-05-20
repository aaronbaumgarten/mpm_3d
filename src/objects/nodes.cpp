//
// Created by aaron on 5/19/17.
// nodes.cpp
//

//
// Created by aaron on 5/14/17.
// points.cpp
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>
#include <sstream>
#include <Eigen/Core>
#include <ctime>

#include "job.hpp"

#include "serializer.hpp"

#include "body.hpp"
#include "nodes.hpp"

#include "stringparser.hpp"

Nodes::Nodes():
        x(0,0),
        u(0,0),
        x_t(0,0),
        diff_x_t(0,0),
        m(0,1),
        mx_t(0,0),
        f(0,0),
        active(0,1)
{}

int Nodes::nodesInit(Job* job, Body* body){
    //point initialization stuff
    //kind of useless right now
    //use for stuff that can't be in constructor
    return 1;
}

void Nodes::nodesWriteFrame(Job* job, Body* body, Serializer* serializer){
    //serializer will use x-position to create format for file
    //serializer.serializerWriteVector(&x, "position")
    serializer->serializerWriteVector(u,"displacement");
    serializer->serializerWriteVector(x_t,"velocity");
    serializer->serializerWriteScalar(m,"mass");
    serializer->serializerWriteVector(mx_t,"momentum");
    serializer->serializerWriteVector(f,"force");
    serializer->serializerWriteScalar(active,"active");

    return;
}

std::string Nodes::nodesSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath) {
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);

    //create filename
    std::ostringstream s;
    s << "mpm_v2." << body->name << "." << body->id << ".nodes." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    std::string filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    //write data
    if (ffile.is_open()) {
        size_t len = x.size();

        ffile << "# mpm_v2 Points\n";
        ffile << "count\n" << len << "\n\n";

        ffile << "x\n";
        ffile << "{";
        job->jobVectorToFile(x,ffile);
        ffile << "}\n\n";

        ffile << "u\n";
        ffile << "{";
        job->jobVectorToFile(u,ffile);
        ffile << "}\n\n";

        ffile << "x_t\n";
        ffile << "{";
        job->jobVectorToFile(x_t,ffile);
        ffile << "}\n\n";

        ffile << "m\n";
        ffile << "{";
        job->jobScalarToFile(m,ffile);
        ffile << "}\n\n";

        ffile << "mx_t\n";
        ffile << "{";
        job->jobVectorToFile(mx_t,ffile);
        ffile << "}\n\n";

        ffile << "f\n";
        job->jobVectorToFile(f,ffile);
        ffile << "}\n\n";

        ffile << "active\n";
        ffile << "{";
        job->jobScalarToFile(active,ffile);
        ffile << "}\n\n";

        ffile.close();
    } else {
        std::cout << "Unable to open .\n";
        return "ERR";
    }
    return filename;
}

int Nodes::nodesLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){

    std::string line; //read line

    std::ifstream fin(fullpath); //file to load from

    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin,line)) {
            //remove spaces and comments from line
            line = StringParser::stringRemoveComments(line);
            line = StringParser::stringRemoveSpaces(line);

            //check if line gives 'count' of points
            if (line.compare("count")) {
                std::getline(fin, line);
                line = StringParser::stringRemoveComments(line);
                line = StringParser::stringRemoveSpaces(line);

                size_t len = std::stoi(line); //count\n #\n

                x = job->jobVectorArray<double>(len);
                u = job->jobVectorArray<double>(len);
                x_t = job->jobVectorArray<double>(len);
                m.resize(len);
                mx_t = job->jobVectorArray<double>(len);
                f = job->jobVectorArray<double>(len);
                active.resize(len);

                x.setZero();
                u.setZero();
                x_t.setZero();
                m.setZero();
                mx_t.setZero();
                f.setZero();
                active.setZero();

                std::vector<std::string> nodeFields = {"x","u","x_t","m","mx_t","f","active"};
                while(std::getline(fin, line)){
                    line = StringParser::stringRemoveComments(line);
                    line = StringParser::stringRemoveSpaces(line);
                    if (line.size() > 0){
                        switch (StringParser::stringFindStringID(nodeFields,line)){
                            case 0:
                                //x
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorFromFile(x,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"x\". Got: " << line << std::endl;
                                }
                                break;
                            case 1:
                                //u
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorFromFile(u,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"u\". Got: " << line << std::endl;
                                }
                                break;
                            case 2:
                                //x_t
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorFromFile(x_t,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"x_t\". Got: " << line << std::endl;
                                }
                                break;
                            case 3:
                                //m
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobScalarFromFile(m,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"m\". Got: " << line << std::endl;
                                }
                                break;
                            case 4:
                                //mx_t
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorFromFile(mx_t,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"mx_t\". Got: " << line << std::endl;
                                }
                                break;
                            case 5:
                                //f
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorFromFile(f,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"b\". Got: " << line << std::endl;
                                }
                                break;
                            case 6:
                                //active
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorFromFile(active,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"active\". Got: " << line << std::endl;
                                }
                                break;
                            default:
                                std::cerr << "Unknown field title: " << line << std::endl;
                        }
                    }
                }
            }
        }
        //close file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    return 1;
}
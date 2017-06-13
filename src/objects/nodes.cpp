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
#include <grid.hpp>

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

    size_t len = job->grid.node_count;

    //initialize vectors to length of grid object
    if (x.rows() != len) {
        x = job->jobVectorArray<double>(len);
        u = job->jobVectorArray<double>(len);
        x_t = job->jobVectorArray<double>(len);
        diff_x_t = job->jobVectorArray<double>(len);
        m.resize(len);
        mx_t = job->jobVectorArray<double>(len);
        f = job->jobVectorArray<double>(len);
        active.resize(len);
    }

    x.setZero();
    u.setZero();
    x_t.setZero();
    diff_x_t.setZero();

    m.setZero();
    mx_t.setZero();
    f.setZero();
    active.setOnes();

    for (size_t i=0;i<len;i++){
        x.row(i) = job->grid.gridNodeIDToPosition(job,i);
    }

    active.setOnes(); //all nodes active

    return 1;
}

void Nodes::nodesWriteFrame(Job* job, Body* body, Serializer* serializer){
    //serializer will use x-position to create format for file
    //serializer.serializerWriteVectorArray(&x, "position")
    serializer->serializerWriteVectorArray(u,"displacement");
    serializer->serializerWriteVectorArray(x_t,"velocity");
    serializer->serializerWriteVectorArray(diff_x_t,"velocity");
    serializer->serializerWriteScalarArray(m,"mass");
    serializer->serializerWriteVectorArray(mx_t,"momentum");
    serializer->serializerWriteVectorArray(f,"force");
    //need to make double
    Eigen::VectorXd tmpVec = active.cast<double>();
    serializer->serializerWriteScalarArray(tmpVec,"active");

    return;
}

std::string Nodes::nodesSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath) {
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename
    std::ostringstream s;
    s << "mpm_v2." << body->name << "." << body->id << ".nodes." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    //write data
    if (ffile.is_open()) {
        size_t len = x.rows();

        ffile << "# mpm_v2 Nodes\n";
        ffile << "count\n" << len << "\n\n";

        ffile << "x\n";
        ffile << "{\n";
        job->jobVectorArrayToFile(x, ffile);
        ffile << "}\n\n";

        ffile << "u\n";
        ffile << "{\n";
        job->jobVectorArrayToFile(u, ffile);
        ffile << "}\n\n";

        ffile << "x_t\n";
        ffile << "{\n";
        job->jobVectorArrayToFile(x_t, ffile);
        ffile << "}\n\n";

        ffile << "diff_x_t\n";
        ffile << "{\n";
        job->jobVectorArrayToFile(diff_x_t, ffile);
        ffile << "}\n\n";

        ffile << "m\n";
        ffile << "{\n";
        job->jobScalarArrayToFile(m, ffile);
        ffile << "}\n\n";

        ffile << "mx_t\n";
        ffile << "{\n";
        job->jobVectorArrayToFile(mx_t, ffile);
        ffile << "}\n\n";

        ffile << "f\n";
        ffile << "{\n";
        job->jobVectorArrayToFile(f, ffile);
        ffile << "}\n\n";

        ffile << "active\n";
        ffile << "{\n";
        job->jobScalarArrayToFile(active, ffile);
        ffile << "}\n\n";

        ffile.close();
    } else {
        std::cout << "Unable to open .\n";
        return "ERR";
    }

    std::cout << "Nodes Saved: [" << body->name << "]." << std::endl;

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
            if (line.compare("count")==0) {
                std::getline(fin, line);
                line = StringParser::stringRemoveComments(line);
                line = StringParser::stringRemoveSpaces(line);

                size_t len = std::stoi(line); //count\n #\n

                x = job->jobVectorArray<double>(len);
                u = job->jobVectorArray<double>(len);
                x_t = job->jobVectorArray<double>(len);
                diff_x_t = job->jobVectorArray<double>(len);
                m.resize(len);
                mx_t = job->jobVectorArray<double>(len);
                f = job->jobVectorArray<double>(len);
                active.resize(len);

                x.setZero();
                u.setZero();
                x_t.setZero();
                diff_x_t.setZero();
                m.setZero();
                mx_t.setZero();
                f.setZero();
                active.setZero();

                std::vector<std::string> nodeFields = {"x","u","x_t","diff_x_t","m","mx_t","f","active","}"};
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
                                    job->jobVectorArrayFromFile(x, fin);
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
                                    job->jobVectorArrayFromFile(u, fin);
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
                                    job->jobVectorArrayFromFile(x_t, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"x_t\". Got: " << line << std::endl;
                                }
                                break;
                            case 3:
                                //diff_x_t
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorArrayFromFile(diff_x_t, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"diff_x_t\". Got: " << line << std::endl;
                                }
                                break;
                            case 4:
                                //m
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobScalarArrayFromFile(m, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"m\". Got: " << line << std::endl;
                                }
                                break;
                            case 5:
                                //mx_t
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorArrayFromFile(mx_t, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"mx_t\". Got: " << line << std::endl;
                                }
                                break;
                            case 6:
                                //f
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorArrayFromFile(f, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"b\". Got: " << line << std::endl;
                                }
                                break;
                            case 7:
                                //active
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobScalarArrayFromFile(active, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"active\". Got: " << line << std::endl;
                                }
                                break;
                            case 8:
                                //}
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

    std::cout << "Nodes Loaded: [" << body->name << "]." << std::endl;

    return 1;
}
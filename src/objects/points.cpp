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
#include "points.hpp"

#include "stringparser.hpp"

Points::Points():
        x(0,0),
        u(0,0),
        x_t(0,0),
        m(0,1),
        v(0,1),
        mx_t(0,0),
        b(0,0),
        T(0,0),
        L(0,0),
        active(0,1),
        extent(0,1)
{
    file = "";
}

int Points::pointsInit(Job* job, Body* body){
    //extent initialization
    if(job->DIM == 1){
        for (size_t i=0;i<v.cols();i++){
            extent[i] = 0.5 * v[i];
        }
    } else if (job->DIM == 2){
        for (size_t i=0;i<v.cols();i++){
            extent[i] = 0.5 * std::sqrt(v[i]);
        }
    } else if (job->DIM == 3){
        for (size_t i = 0; i < v.cols(); i++) {
            extent[i] = 0.5 * std::cbrt(v[i]);
        }
    }
    return 1;
}

int Points::pointsReadFromFile(Job *job, Body *body, std::string fileIN) {
    //first line lists number of points in file
    //subsequent lines contain m, v, x, x_t, active
    file = fileIN;

    std::string line;
    std::ifstream fin(file);
    std::stringstream ss;

    if (fin.is_open()){
        std::getline(fin,line);
        size_t len = std::stoi(line);

        x = job->jobVectorArray<double>(len);
        u = job->jobVectorArray<double>(len);
        x_t = job->jobVectorArray<double>(len);
        m.resize(len);
        v.resize(len);
        mx_t = job->jobVectorArray<double>(len);
        b = job->jobVectorArray<double>(len);
        T = job->jobTensorArray<double>(len);
        L = job->jobTensorArray<double>(len);
        active.resize(len);
        extent.resize(len);

        x.setZero();
        u.setZero();
        x_t.setZero();
        m.setZero();
        v.setZero();
        mx_t.setZero();
        b.setZero();
        T.setZero();
        L.setZero();
        active.setZero();
        extent.setZero();

        for (size_t i=0; i<len; i++){
            std::getline(fin,line);
            ss = std::stringstream(line);
            ss >> m[i];
            ss >> v[i];
            for (size_t pos=0; pos < job->DIM; pos++){
                ss >> x(i,pos);
            }
            for (size_t pos=0; pos < job->DIM; pos++){
                ss >> x_t(i,pos);
            }
            if(!(ss >> active[i])){
                std::cerr << "ERROR: Unable to read line: " << file << std::endl;
                return 0;
            }
        }

    } else {
        std::cerr << "ERROR: Unable to open file: " << file << std::endl;
        return 0;
    }

    return 1;
}

void Points::pointsWriteFrame(Job* job, Body* body, Serializer* serializer){
    //serializer will use x-position to create format for file
    //serializer.serializerWriteVector(&x, "position")
    serializer->serializerWriteVector(u,"displacement");
    serializer->serializerWriteVector(x_t,"velocity");
    serializer->serializerWriteScalar(m,"mass");
    serializer->serializerWriteScalar(v,"volume");
    serializer->serializerWriteVector(mx_t,"momentum");
    serializer->serializerWriteVector(b,"body_force");
    serializer->serializerWriteTensor(T,"cauchy_stress");
    serializer->serializerWriteTensor(L,"velocity_gradient");
    serializer->serializerWriteScalar(active,"active");
    serializer->serializerWriteScalar(extent,"extent");

    return;
}

std::string Points::pointsSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath) {
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);

    //create filename
    std::ostringstream s;
    s << "mpm_v2."  << body->name << "." << body->id << ".points." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    std::string filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    //write data
    if (ffile.is_open()) {
        size_t len = x.rows();

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

        ffile << "v\n";
        ffile << "{";
        job->jobScalarToFile(v,ffile);
        ffile << "}\n\n";

        ffile << "mx_t\n";
        ffile << "{";
        job->jobVectorToFile(mx_t,ffile);
        ffile << "}\n\n";

        ffile << "b\n";
        job->jobVectorToFile(b,ffile);
        ffile << "}\n\n";

        ffile << "T\n";
        ffile << "{";
        job->jobTensorToFile(T,ffile);
        ffile << "}\n\n";

        ffile << "L\n";
        ffile << "{";
        job->jobTensorToFile(L,ffile);
        ffile << "}\n\n";

        ffile << "active\n";
        ffile << "{";
        job->jobScalarToFile(active,ffile);
        ffile << "}\n\n";

        ffile << "extent\n";
        ffile << "{";
        job->jobScalarToFile(extent,ffile);
        ffile << "}\n\n";

        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filename << "\" !\n";
        return "ERR";
    }
    return filename;
}

int Points::pointsLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){

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
                v.resize(len);
                mx_t = job->jobVectorArray<double>(len);
                b = job->jobVectorArray<double>(len);
                T = job->jobTensorArray<double>(len);
                L = job->jobTensorArray<double>(len);
                active.resize(len);
                extent.resize(len);

                x.setZero();
                u.setZero();
                x_t.setZero();
                m.setZero();
                v.setZero();
                mx_t.setZero();
                b.setZero();
                T.setZero();
                L.setZero();
                active.setZero();
                extent.setZero();

                std::vector<std::string> pointFields = {"x","u","x_t","m","v","mx_t","b","T","L","active","extent"};
                while(std::getline(fin, line)){
                    line = StringParser::stringRemoveComments(line);
                    line = StringParser::stringRemoveSpaces(line);
                    if (line.size() > 0){
                        switch (StringParser::stringFindStringID(pointFields,line)){
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
                                //v
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobScalarFromFile(v,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"v\". Got: " << line << std::endl;
                                }
                                break;
                            case 5:
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
                            case 6:
                                //b
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorFromFile(b,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"b\". Got: " << line << std::endl;
                                }
                                break;
                            case 7:
                                //T
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobTensorFromFile(T,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"T\". Got: " << line << std::endl;
                                }
                                break;
                            case 8:
                                //L
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobTensorFromFile(L,fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"L\". Got: " << line << std::endl;
                                }
                                break;
                            case 9:
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
                            case 10:
                                //extent
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorFromFile(extent,fin);
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

    std::cout << "Points Loaded: [" << body->id << "]." << std::endl;

    return 1;
}
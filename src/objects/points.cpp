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
        v0(0,1),
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
    //v0 initialization
    v0 = v;

    //extent initialization
    if(job->DIM == 1){
        for (size_t i=0;i<v.rows();i++){
            extent[i] = 0.5 * v[i];
        }
    } else if (job->DIM == 2){
        for (size_t i=0;i<v.rows();i++){
            extent[i] = 0.5 * std::sqrt(v[i]);
        }
    } else if (job->DIM == 3){
        for (size_t i = 0; i < v.rows(); i++) {
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
        v0.resize(len);
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
        v0.setZero();
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
    //serializer.serializerWriteVectorArray(&x, "position")
    serializer->serializerWriteVectorArray(u,"displacement");
    serializer->serializerWriteVectorArray(x_t,"velocity");
    serializer->serializerWriteScalarArray(m,"mass");
    serializer->serializerWriteScalarArray(v,"volume");
    serializer->serializerWriteVectorArray(mx_t,"momentum");
    serializer->serializerWriteVectorArray(b,"body_force");
    serializer->serializerWriteTensorArray(T,"cauchy_stress");
    serializer->serializerWriteTensorArray(L,"velocity_gradient");
    //need to make double
    Eigen::VectorXd tmpVec = active.cast<double>();
    serializer->serializerWriteScalarArray(tmpVec,"active");
    serializer->serializerWriteScalarArray(extent,"extent");

    //pressure
    if (job->DIM == 1) {
        tmpVec = -1.0*T;
    } else if (job->DIM == 2) {
        tmpVec = -1.0/2.0 * (T.col(0) + T.col(3));
    } else if (job->DIM == 3) {
        tmpVec = -1.0/3.0 * (T.col(job->XX) + T.col(job->YY) + T.col(job->ZZ));
    }
    serializer->serializerWriteScalarArray(tmpVec,"pressure");

    tmpVec = m.array() / v.array();
    serializer->serializerWriteScalarArray(tmpVec,"density");

    return;
}

std::string Points::pointsSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath) {
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename
    std::stringstream s;
    s << "mpm_v2."  << body->name << "." << body->id << ".points." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    //write data
    if (ffile.is_open()) {
        size_t len = x.rows();

        ffile << "# mpm_v2 Points\n";
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

        ffile << "m\n";
        ffile << "{\n";
        job->jobScalarArrayToFile(m, ffile);
        ffile << "}\n\n";

        ffile << "v\n";
        ffile << "{\n";
        job->jobScalarArrayToFile(v, ffile);
        ffile << "}\n\n";

        ffile << "v0\n";
        ffile << "{\n";
        job->jobScalarArrayToFile(v0, ffile);
        ffile << "}\n\n";

        ffile << "mx_t\n";
        ffile << "{\n";
        job->jobVectorArrayToFile(mx_t, ffile);
        ffile << "}\n\n";

        ffile << "b\n";
        ffile << "{\n";
        job->jobVectorArrayToFile(b, ffile);
        ffile << "}\n\n";

        ffile << "T\n";
        ffile << "{\n";
        job->jobTensorArrayToFile(T, ffile);
        ffile << "}\n\n";

        ffile << "L\n";
        ffile << "{\n";
        job->jobTensorArrayToFile(L, ffile);
        ffile << "}\n\n";

        ffile << "active\n";
        ffile << "{\n";
        job->jobScalarArrayToFile(active, ffile);
        ffile << "}\n\n";

        ffile << "extent\n";
        ffile << "{\n";
        job->jobScalarArrayToFile(extent, ffile);
        ffile << "}\n\n";

        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Points Saved: [" << body->name << "]." << std::endl;

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
            if (line.compare("count") == 0) {
                std::getline(fin, line);
                line = StringParser::stringRemoveComments(line);
                line = StringParser::stringRemoveSpaces(line);

                size_t len = std::stoi(line); //count\n #\n

                x = job->jobVectorArray<double>(len);
                u = job->jobVectorArray<double>(len);
                x_t = job->jobVectorArray<double>(len);
                m.resize(len);
                v.resize(len);
                v0.resize(len);
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
                v0.setZero();
                mx_t.setZero();
                b.setZero();
                T.setZero();
                L.setZero();
                active.setZero();
                extent.setZero();

                std::vector<std::string> pointFields = {"x","u","x_t","m","v","v0","mx_t","b","T","L","active","extent","}"};
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
                            case 4:
                                //v
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobScalarArrayFromFile(v, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"v\". Got: " << line << std::endl;
                                }
                                break;
                            case 5:
                                //v0
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobScalarArrayFromFile(v0, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"v0\". Got: " << line << std::endl;
                                }
                                break;
                            case 6:
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
                            case 7:
                                //b
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobVectorArrayFromFile(b, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"b\". Got: " << line << std::endl;
                                }
                                break;
                            case 8:
                                //T
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobTensorArrayFromFile(T, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"T\". Got: " << line << std::endl;
                                }
                                break;
                            case 9:
                                //L
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobTensorArrayFromFile(L, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"L\". Got: " << line << std::endl;
                                }
                                break;
                            case 10:
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
                            case 11:
                                //extent
                                std::getline(fin, line);
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);
                                if (line.compare("{") == 0){
                                    job->jobScalarArrayFromFile(extent, fin);
                                } else {
                                    std::cerr << "Expected \"{\" symbol after \"extent\". Got: " << line << std::endl;
                                }
                                break;
                            case 12:
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

    std::cout << "Points Loaded: [" << body->name << "]." << std::endl;

    return 1;
}
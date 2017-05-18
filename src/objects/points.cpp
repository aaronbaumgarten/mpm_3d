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

#include "stringparse.hpp"

Points::Points():
        x(0,0),
        u(0,0),
        x_t(0,0),
        m(0),
        mx_t(0,0),
        b(0,0),
        T(0,0),
        L(0,0),
        active(0)
{}

int Points::pointsInit(Job* job, Body* body){
    //point initialization stuff
    //kind of useless right now
    //use for stuff that can't be in constructor
    return 1;
}

void Points::pointsWriteFrame(Job* job, Body* body, Serializer* serializer){
    //serializer will use x-position to create format for file
    //serializer.serializerWriteVector(&x, "position")
    serializer->serializerWriteVector(&u,"displacement");
    serializer->serializerWriteVector(&x_t,"velocity");
    serializer->serializerWriteScalar(&m,"mass");
    serializer->serializerWriteVector(&mx_t,"momentum");
    serializer->serializerWriteVector(&b,"body_force");
    serializer->serializerWriteTensor(&T,"cauchy_stress");
    serializer->serializerWriteTensor(&T,"velocity_gradient");
    serializer->serializerWriteScalar(&active,"active");

    return;
}

std::string Points::pointsSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath) {
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);

    //create filename
    std::ostringstream s;
    s << "mpm_v2.points." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
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
        for (size_t i = 0; i < len; i++) {
            ffile << x(i,0) << " " << x(i,1) << " " << x(i,2) << "\n";
        }
        ffile << "}\n\n";

        ffile << "u\n";
        ffile << "{";
        for (size_t i = 0; i < len; i++) {
            ffile << u(i,0) << " " << u(i,1) << " " << u(i,2) << "\n";
        }
        ffile << "}\n\n";

        ffile << "x_t\n";
        ffile << "{";
        for (size_t i = 0; i < len; i++) {
            ffile << x_t(i,0) << " " << x_t(i,1) << " " << x_t(i,2) << "\n";
        }
        ffile << "}\n\n";

        ffile << "m\n";
        ffile << "{";
        for (size_t i = 0; i < len; i++) {
            ffile << m(i) << "\n";
        }
        ffile << "}\n\n";

        ffile << "mx_t\n";
        ffile << "{";
        for (size_t i = 0; i < len; i++) {
            ffile << mx_t(i,0) << " " << mx_t(i,1) << " " << mx_t(i,2) << "\n";
        }
        ffile << "}\n\n";

        ffile << "b\n";
        ffile << "{";
        for (size_t i = 0; i < len; i++) {
            ffile << b(i,0) << " " << b(i,1) << " " << b(i,2) << "\n";
        }
        ffile << "}\n\n";

        ffile << "T\n";
        ffile << "{";
        for (size_t i = 0; i < len; i++) {
            ffile << T(i,0) << " " << T(i,1) << " " << T(i,2) << " ";
            ffile << T(i,3) << " " << T(i,4) << " " << T(i,5) << " ";
            ffile << T(i,6) << " " << T(i,7) << " " << T(i,8) << "\n";
        }
        ffile << "}\n\n";

        ffile << "L\n";
        ffile << "{";
        for (size_t i = 0; i < len; i++) {
            ffile << L(i,0) << " " << L(i,1) << " " << L(i,2) << " ";
            ffile << L(i,3) << " " << L(i,4) << " " << L(i,5) << " ";
            ffile << L(i,6) << " " << L(i,7) << " " << L(i,8) << "\n";
        }
        ffile << "}\n\n";

        ffile << "active\n";
        ffile << "{";
        for (size_t i = 0; i < len; i++) {
            ffile << active(i) << "\n";
        }
        ffile << "}\n\n";

        ffile.close();
    } else {
        std::cout << "Unable to open .\n";
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
            line = StringParse::stringRemoveComments(line);
            line = StringParse::stringRemoveSpaces(line);

            //check if line gives 'count' of points
            if (line.compare("count")) {
                std::getline(fin, line);
                line = StringParse::stringRemoveComments(line);
                line = StringParse::stringRemoveSpaces(line);

                size_t len = std::stoi(line); //count\n #\n

                x = job->jobVectorArray<double>(len);
                u = job->jobVectorArray<double>(len);
                x_t = job->jobVectorArray<double>(len);
                m = Eigen::VectorXd(len);
                mx_t = job->jobVectorArray<double>(len);
                b = job->jobVectorArray<double>(len);
                T = job->jobTensorArray<double>(len);
                L = job->jobTensorArray<double>(len);
                active = Eigen::VectorXi(len);

                x.setZero();
                u.setZero();
                x_t.setZero();
                m.setZero();
                mx_t.setZero();
                b.setZero();
                T.setZero();
                L.setZero();
                active.setZero();
/////////////////////////////////////////////////////////////////////////


                //check if line matches desired header
                if (line.compare("job") == 0) {
                    std::getline(fin, line);
                    line = this->removeComments(line);
                    line = this->removeSpaces(line);
                    //check if next line is open bracket
                    if (line.compare("{") == 0) {
                        while (std::getline(fin, line)) {
                            line = this->removeComments(line);
                            line = this->removeSpaces(line);
                            if (line.compare("}") == 0) {
                                //stop reading headers
                                break;
                            }
                            lvec = this->splitString(line, '=');
                            if (lvec.size() > 1) {
                                propName = lvec[0];
                                propValue = lvec[1];

                                //check if property name matches job inputs
                                switch (this->findStringID(this->jobParams, propName)) {
                                    //cases for setting job values need to process strings
                                    //{"dt","use_3d","use_implicit","use_cpdi","newtonTOL","linearStepSize","use_smoothing"}
                                    case 0 :
                                        job->dt = std::stod(propValue);
                                        break;
                                    case 1 :
                                        job->use_3d = std::stoi(propValue);
                                        break;
                                    case 2 :
                                        job->use_implicit = std::stoi(propValue);
                                        break;
                                    case 3 :
                                        job->use_cpdi = std::stoi(propValue);
                                        break;
                                    case 4 :
                                        job->newtonTOL = std::stod(propValue);
                                        break;
                                    case 5 :
                                        job->linearStepSize = std::stod(propValue);
                                        break;
                                    case 6 :
                                        job->use_smoothing = std::stoi(propValue);
                                        break;
                                    default:
                                        std::cout << "Parameter \"" << propName <<
                                                  "\" not recognized for \"job\" configuration." << std::endl;
                                }
                            }
                        }
                    }
                    std::cout << "Job configured." << std::endl;
                    //read only first job entry (stop reading doc)
                    break;
                }
            }
        }
        //close file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << this->configFile << std::endl;
        return 0;
    }


    return -1;
}
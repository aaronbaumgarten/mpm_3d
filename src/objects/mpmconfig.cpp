//
// Created by aaron on 12/26/16.
// mpmconfig.cpp
//
#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include "process.hpp"
#include "mpmconfig.hpp"

MPMconfig::MPMconfig(){
    configFile = "default.cfg";
    //knownHeaders = {"job","input","output","material","boundary","contact"};
    jobParams = {"dt","use_3d","use_implicit","use_cpdi","newtonTOL","linearStepSize"};
    inputParams = {"particle-file","grid-file"};
    materialParams = {"bodies","material-file","properties","int-properties"};
    boundaryParams = {"boundary-file","properties","int-properties"};
    contactParams = {"bodies","contact-file","properties","int-properties"};
    outputParams = {"frame-rate","simulation-time","output-file","output-dir"};
}

MPMconfig::MPMconfig(std::string fileIn) {
    configFile = fileIn;
    jobParams = {"dt","use_3d","use_implicit","use_cpdi","newtonTOL","linearStepSize"};
    inputParams = {"particle-file","grid-file"};
    materialParams = {"bodies","material-file","properties","int-properties"};
    boundaryParams = {"boundary-file","properties","int-properties"};
    contactParams = {"bodies","contact-file","properties","int-properties"};
    outputParams = {"frame-rate","simulation-time","output-file","output-dir"};
}

void MPMconfig::setConfigFile(std::string fileIn){
    this->configFile = fileIn;
    return;
}

std::string MPMconfig::removeSpaces(std::string s){
    s.erase(std::remove_if(s.begin(), s.end(), isspace), s.end());
    return s;
}

std::vector<std::string> MPMconfig::splitString(std::string s,char delim){
    //credit to Evan Teran
    //found on stack exchange
    std::stringstream ss;
    std::string tmp;
    std::vector<std::string> elems;
    ss.str(s);
    while (std::getline(ss, tmp, delim)){
        elems.push_back(tmp);
    }
    return elems;
}

std::string MPMconfig::removeComments(std::string s) {
    std::vector<std::string> svec;
    svec = this->splitString(s,'#');
    if (svec.size()>0) {
        return svec[0];
    } else {
        return std::string();
    }
}

std::string MPMconfig::removeBraces(std::string s){
    std::vector<std::string> svec;
    std::stringstream ss;
    //std::cout << s << std::endl;
    svec = this->splitString(s,'{');
    for(size_t i=0;i<svec.size();i++){
        ss << svec[i];
    }
    //std::cout << ss.str() << std::endl;
    svec = this->splitString(ss.str(),'}');
    ss.str("");
    ss.clear();
    //std::cout << ss.str() << std::endl;
    for(size_t i=0;i<svec.size();i++){
        ss << svec[i];
    }
    //std::cout << ss.str() << std::endl;
    return ss.str();
}

int MPMconfig::findStringID(std::vector<std::string> svec, std::string s){
    for (size_t i=0;i<svec.size();i++){
        if(svec[i].compare(s) == 0){
            return i;
        }
    }
    return -1;
}

int MPMconfig::checkConfigFile(std::string filename) {
    std::vector<std::string> headers = {"job","input","material","boundary","contact","output"};
    std::vector<int> check = {0,0,0,0,0,0};
    size_t id;
    std::string line;
    std::ifstream fin(filename);
    if (fin.is_open()){
        while (std::getline(fin,line)){
            line = this->removeComments(line);
            line = this->removeSpaces(line);
            id = this->findStringID(headers,line);
            if (id != -1){
                check[id] = 1;
            }
        }
        for (size_t i=0;i<headers.size();i++){
            if (check[i]==0){
                std::cout << "Could not find \"" << headers[i] << "\" in configuration file." << std::endl;
            }
        }
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << filename << std::endl;
    }
}

int MPMconfig::configJob(job_t * job){
    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::string propName;
    std::string propValue;
    std::vector<double> fp64_props;
    std::vector<int> int_props;

    //declare fstream from filename
    std::ifstream fin(this->configFile);
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin,line)){
            //remove spaces and comments from line
            line = this->removeComments(line);
            line = this->removeSpaces(line);
            //check if line matches desired header
            if (line.compare("job") == 0){
                std::getline(fin,line);
                line = this->removeComments(line);
                line = this->removeSpaces(line);
                //check if next line is open bracket
                if (line.compare("{") == 0){
                    while (std::getline(fin,line)){
                        line = this->removeComments(line);
                        line = this->removeSpaces(line);
                        if (line.compare("}") == 0) {
                            //stop reading headers
                            break;
                        }
                        lvec = this->splitString(line,'=');
                        if (lvec.size()>1) {
                            propName = lvec[0];
                            propValue = lvec[1];

                            //check if property name matches job inputs
                            switch (this->findStringID(this->jobParams, propName)) {
                                //cases for setting job values need to process strings
                                //{"dt","use_3d","use_implicit","use_cpdi","newtonTOL","linearStepSize"}
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
        //close file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << this->configFile << std::endl;
        return 0;
    }

    return 1;
}

int MPMconfig::configInput(job_t * job){
    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::string propName;
    std::string propValue;
    std::string nfilename;
    std::string pfilename;

    //declare fstream from filename
    std::ifstream fin(this->configFile);
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin,line)){
            //remove spaces and comments from line
            line = this->removeComments(line);
            line = this->removeSpaces(line);
            //check if line matches desired header
            if (line.compare("input") == 0){
                std::getline(fin,line);
                line = this->removeComments(line);
                line = this->removeSpaces(line);
                //check if next line is open bracket
                if (line.compare("{") == 0){
                    while (std::getline(fin,line)){
                        line = this->removeComments(line);
                        line = this->removeSpaces(line);
                        if (line.compare("}") == 0) {
                            //stop reading headers
                            break;
                        }

                        lvec = this->splitString(line,'=');
                        if (lvec.size()>1) {
                            propName = lvec[0];
                            propValue = lvec[1];

                            //check if property name matches input inputs
                            switch (this->findStringID(this->inputParams, propName)) {
                                //cases for setting input values need to process strings
                                //{"particle-file","grid-file"}
                                case 0 :
                                    pfilename = propValue;
                                    break;
                                case 1 :
                                    nfilename = propValue;
                                    break;
                                default:
                                    std::cout << "Parameter \"" << propName <<
                                    "\" not recognized for \"input\" configuration." << std::endl;
                            }
                        }
                    }
                }
                if (nfilename.size() != 0 && pfilename.size() != 0){
                    std::cout << "Input configured. particles: " << pfilename << " nodes: " << nfilename << std::endl;
                    if (job->use_3d == 1){
                        fin.close();
                        return job->importNodesandParticles(nfilename,pfilename);
                    } else {
                        fin.close();
                        return job->importNodesandParticles2D(nfilename,pfilename);
                    }
                } else {
                    std::cout << "Cannot configure inputs without particle-file or grid-file. Exiting." << std::endl;
                    fin.close();
                    return 0;
                }
                break;
            }
        }
        //close file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << this->configFile << std::endl;
        return 0;
    }

    return 1;
}

int MPMconfig::configMaterial(job_t * job){
    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::string propName;
    std::string propValue;
    std::string mfilename;
    std::vector<int> bodyIDs;
    std::vector<double> fp64_props;
    std::vector<int> int_props;

    //declare fstream from filename
    std::ifstream fin(this->configFile);
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin,line)){
            //remove spaces and comments from line
            line = this->removeComments(line);
            line = this->removeSpaces(line);
            //check if line matches desired header
            if (line.compare("material") == 0){
                std::getline(fin,line);
                line = this->removeComments(line);
                line = this->removeSpaces(line);
                //check if next line is open bracket
                if (line.compare("{") == 0){
                    while (std::getline(fin,line)){
                        line = this->removeComments(line);
                        line = this->removeSpaces(line);
                        if (line.compare("}") == 0) {
                            //stop reading headers
                            break;
                        }

                        line = this->removeBraces(line);
                        lvec = this->splitString(line,'=');
                        if (lvec.size()>1) {
                            propName = lvec[0];
                            propValue = lvec[1];

                            std::vector<std::string> ssplit = this->splitString(propValue, ',');
                            //check if property name matches material inputs
                            switch (this->findStringID(this->materialParams, propName)) {
                                //cases for setting material values need to process strings
                                //{"bodies","material-file","properties","int-properties"}
                                case 0 :
                                    //read body ids into integer vector
                                    for (size_t i = 0; i < ssplit.size(); i++) {
                                        bodyIDs.push_back(std::stoi(ssplit[i]));
                                    }
                                    break;
                                case 1 :
                                    mfilename = propValue;
                                    break;
                                case 2 :
                                    for (size_t i = 0; i < ssplit.size(); i++) {
                                        fp64_props.push_back(std::stod(ssplit[i]));
                                    }
                                    break;
                                case 3 :
                                    for (size_t i = 0; i < ssplit.size(); i++) {
                                        int_props.push_back(std::stoi(ssplit[i]));
                                    }
                                    break;
                                default:
                                    std::cout << "Parameter \"" << propName <<
                                    "\" not recognized for \"material\" configuration." << std::endl;
                            }
                        }
                    }
                }
                if (mfilename.size() != 0 && bodyIDs.size() != 0){
                    std::cout << "Material configured. material: " << mfilename << std::endl;
                    for (size_t i=0;i<bodyIDs.size();i++){
                        job->assignMaterial(mfilename,bodyIDs[i],fp64_props,int_props);
                    }
                } else {
                    std::cout << "Cannot configure material without material-file or bodies." << std::endl;
                }
                //reset arrays and vectors
                mfilename.clear();
                bodyIDs.clear();
                fp64_props.clear();
                int_props.clear();
            }
        }
        //close file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << this->configFile << std::endl;
        return 0;
    }

    for (size_t b=0; b<job->num_bodies; b++){
        if(job->bodies[b].material.material_filename.size() == 0){
            std::cout << "Failed to configure all materials. Body " << b << " unassigned. Exiting." << std::endl;
            return 0;
        }
    }

    return 1;
}

int MPMconfig::configBoundary(job_t * job){
    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::string propName;
    std::string propValue;
    std::string bfilename;
    std::vector<double> fp64_props;
    std::vector<int> int_props;

    //declare fstream from filename
    std::ifstream fin(this->configFile);
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin,line)){
            //remove spaces and comments from line
            line = this->removeComments(line);
            line = this->removeSpaces(line);
            //check if line matches desired header
            if (line.compare("boundary") == 0){
                std::getline(fin,line);
                line = this->removeComments(line);
                line = this->removeSpaces(line);
                //check if next line is open bracket
                if (line.compare("{") == 0){
                    while (std::getline(fin,line)){
                        line = this->removeComments(line);
                        line = this->removeSpaces(line);
                        if (line.compare("}") == 0) {
                            //stop reading headers
                            break;
                        }

                        line = this->removeBraces(line);
                        lvec = this->splitString(line,'=');
                        if (lvec.size()>1) {
                            propName = lvec[0];
                            propValue = lvec[1];

                            std::vector<std::string> ssplit = this->splitString(propValue, ',');
                            //check if property name matches boundary inputs
                            switch (this->findStringID(this->boundaryParams, propName)) {
                                //cases for setting boundary values need to process strings
                                //{"boundary-file","properties","int-properties"}
                                case 0 :
                                    bfilename = propValue;
                                    break;
                                case 1 :
                                    for (size_t i = 0; i < ssplit.size(); i++) {
                                        fp64_props.push_back(std::stod(ssplit[i]));
                                    }
                                    break;
                                case 2 :
                                    for (size_t i = 0; i < ssplit.size(); i++) {
                                        int_props.push_back(std::stoi(ssplit[i]));
                                    }
                                    break;
                                default:
                                    std::cout << "Parameter \"" << propName <<
                                    "\" not recognized for \"boundary\" configuration." << std::endl;
                            }
                        }
                    }
                }
                if (bfilename.size() != 0){
                    std::cout << "Boundary configured. boundary: " << bfilename << std::endl;
                    job->assignBoundaryConditions(bfilename,fp64_props,int_props);
                } else {
                    std::cout << "Cannot configure boundary without boundary-file." << std::endl;
                }
                //reset arrays and vectors
                bfilename.clear();
                fp64_props.clear();
                int_props.clear();
            }
        }
        //close file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << this->configFile << std::endl;
        return 0;
    }

    return 1;
}

int MPMconfig::configContact(job_t * job){
    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::string propName;
    std::string propValue;
    std::string cfilename;
    std::vector<int> bodyIDs;
    std::vector<double> fp64_props;
    std::vector<int> int_props;

    //declare fstream from filename
    std::ifstream fin(this->configFile);
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin,line)){
            //remove spaces and comments from line
            line = this->removeComments(line);
            line = this->removeSpaces(line);
            //check if line matches desired header
            if (line.compare("contact") == 0){
                std::getline(fin,line);
                line = this->removeComments(line);
                line = this->removeSpaces(line);
                //check if next line is open bracket
                if (line.compare("{") == 0){
                    while (std::getline(fin,line)){
                        line = this->removeComments(line);
                        line = this->removeSpaces(line);
                        if (line.compare("}") == 0) {
                            //stop reading headers
                            break;
                        }

                        line = this->removeBraces(line);
                        lvec = this->splitString(line,'=');
                        if (lvec.size()>1) {
                            propName = lvec[0];
                            propValue = lvec[1];

                            std::vector<std::string> ssplit = this->splitString(propValue, ',');
                            //check if property name matches contact inputs
                            switch (this->findStringID(this->contactParams, propName)) {
                                //cases for setting contact values need to process strings
                                //{"bodies","contact-file","properties","int-properties"}
                                case 0 :
                                    //read body ids into integer vector
                                    for (size_t i = 0; i < ssplit.size(); i++) {
                                        bodyIDs.push_back(std::stoi(ssplit[i]));
                                    }
                                    break;
                                case 1 :
                                    cfilename = propValue;
                                    break;
                                case 2 :
                                    for (size_t i = 0; i < ssplit.size(); i++) {
                                        fp64_props.push_back(std::stod(ssplit[i]));
                                    }
                                    break;
                                case 3 :
                                    for (size_t i = 0; i < ssplit.size(); i++) {
                                        int_props.push_back(std::stoi(ssplit[i]));
                                    }
                                    break;
                                default:
                                    std::cout << "Parameter \"" << propName <<
                                    "\" not recognized for \"contact\" configuration." << std::endl;
                            }
                        }
                    }
                }
                if (cfilename.size() != 0 && bodyIDs.size() == 2){
                    std::cout << "Contact configured. contact: " << cfilename << std::endl;
                    job->assignContact(cfilename,bodyIDs,fp64_props,int_props);
                } else {
                    std::cout << "Cannot configure contact without contact-file or bodies." << std::endl;
                }
                //reset arrays and vectors
                cfilename.clear();
                bodyIDs.clear();
                fp64_props.clear();
                int_props.clear();
            }
        }
        //close file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << this->configFile << std::endl;
        return 0;
    }

    return 1;
}

int MPMconfig::configOutput(job_t *job, MPMio *mpmout){
    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::string propName;
    std::string propValue;
    std::vector<double> fp64_props;
    std::vector<int> int_props;

    //declare fstream from filename
    std::ifstream fin(this->configFile);
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin,line)){
            //remove spaces and comments from line
            line = this->removeComments(line);
            line = this->removeSpaces(line);
            //check if line matches desired header
            if (line.compare("output") == 0){
                std::getline(fin,line);
                line = this->removeComments(line);
                line = this->removeSpaces(line);
                //check if next line is open bracket
                if (line.compare("{") == 0){
                    while (std::getline(fin,line)){
                        line = this->removeComments(line);
                        line = this->removeSpaces(line);
                        if (line.compare("}") == 0) {
                            //stop reading headers
                            break;
                        }

                        line = this->removeBraces(line);
                        lvec = this->splitString(line,'=');
                        if (lvec.size()>1) {
                            propName = lvec[0];
                            propValue = lvec[1];

                            //check if property name matches output inputs
                            switch (this->findStringID(this->outputParams, propName)) {
                                //cases for setting output values need to process strings
                                //{"frame-rate","simulation-time","output-file","output-dir"}
                                case 0 :
                                    mpmout->setSampleRate(std::stod(propValue));
                                    break;
                                case 1 :
                                    mpmout->setSampleTime(std::stod(propValue));
                                    break;
                                case 2 :
                                    mpmout->defineFrameFile(propValue);
                                    break;
                                case 3 :
                                    mpmout->defineFrameDirectory(propValue);
                                    break;
                                default:
                                    std::cout << "Parameter \"" << propName <<
                                    "\" not recognized for \"output\" configuration." << std::endl;
                            }
                        }
                    }
                }
                std::cout << "Output configured." << std::endl;
                //read only first job entry (stop reading doc)
                break;
            }
        }
        mpmout->setJob(job);
        //close file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << this->configFile << std::endl;
        return 0;
    }

    return 1;
}
//
// Created by aaron on 5/14/17.
// config.cpp
//


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "stringparser.hpp"

#include "job.hpp"
#include "config.hpp"

#include "serializer.hpp"
#include "driver.hpp"
#include "solver.hpp"

#include "grid.hpp"
#include "contact.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "points.hpp"

#include "material.hpp"
#include "boundary.hpp"

Config::Config() {
    file= "";
    mainpath = "";
}

int Config::configInit(std::string fileIN){
    file= fileIN;
    return 1;
}

void Config::configSetMainPath(std::string program){
    std::vector<std::string> svec;
    svec = StringParser::stringSplitString(program,'/');
    std::string filepath = "";
    for (size_t i=0; i<(svec.size()-1);i++){
        filepath += svec[i];
        filepath += "/";
    }
    mainpath = filepath;
    return;
}

void Config::configCheckConfigFile(std::string filename) {
    std::vector<std::string> headers = {"job","serializer","driver","solver",
                                        "body","contact","grid"};
    std::vector<int> check = {0,0,0,0,0,0,0};
    size_t id;
    std::string line;
    std::ifstream fin(filename);
    if (fin.is_open()){
        while (std::getline(fin,line)){
            line = StringParser::stringRemoveComments(line);
            line = StringParser::stringRemoveSpaces(line);
            id = StringParser::stringFindStringID(headers,line);
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
    return;
}

int Config::configConfigureJob(Job *job){
    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::vector<std::string> params;
    std::string propName;
    std::string propValue;
    std::vector<double> fp64_props;
    std::vector<int> int_props;
    size_t id; //for bodies and contacts

    //section headers
    std::vector<std::string> headers = {"job","serializer","driver","solver",
                                        "body","contact","grid"};

    //declare fstream from filename
    std::ifstream fin(file);

    //read config file and configure
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin,line)){
            //remove spaces and comments from line
            line = StringParser::stringRemoveComments(line);
            line = StringParser::stringRemoveSpaces(line);

            //check if line matches header
            if (line.size()>0){
                switch (StringParser::stringFindStringID(headers,line)) {
                    case 0:
                        //job
                        params = {"dt","t","DIM"};
                        line = StringParser::stringRemoveComments(line);
                        line = StringParser::stringRemoveSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = StringParser::stringRemoveBraces(line);
                                lvec = StringParser::stringSplitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];
                                    switch (StringParser::stringFindStringID(params, propName)) {
                                        case 0:
                                            //dt
                                            job->dt = std::stod(propValue);
                                            break;
                                        case 1:
                                            //t
                                            job->t = std::stod(propValue);
                                            break;
                                        case 3:
                                            //DIM
                                            job->DIM = std::stod(propValue);
                                            if (job->DIM > 4 || job->DIM < 1){
                                                std::cerr << "DIM ERROR!" << std::endl;
                                                return 0;
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                        }
                        break;
                    case 1:
                        //serializer
                        params = {"filepath","filename","properties","int-properties"};
                        line = StringParser::stringRemoveComments(line);
                        line = StringParser::stringRemoveSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = StringParser::stringRemoveBraces(line);
                                lvec = StringParser::stringSplitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = StringParser::stringSplitString(propValue, ',');
                                    switch (StringParser::stringFindStringID(params, propName)) {
                                        case 0:
                                            //filepath
                                            job->serializer.filepath = propValue;
                                            break;
                                        case 1:
                                            //filename
                                            job->serializer.filename = propValue;
                                            break;
                                        case 3:
                                            //props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->serializer.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 4:
                                            //int-props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->serializer.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            job->serializer.serializerSetPlugin(job,job->serializer.filepath,
                                                                job->serializer.filename,
                                                                job->serializer.fp64_props,
                                                                job->serializer.int_props);
                            std::cout << "Serializer Configured." << std::endl;
                        }
                        break;
                    case 2:
                        //driver
                        params = {"filepath","filename","properties","int-properties"};
                        line = StringParser::stringRemoveComments(line);
                        line = StringParser::stringRemoveSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = StringParser::stringRemoveBraces(line);
                                lvec = StringParser::stringSplitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = StringParser::stringSplitString(propValue, ',');
                                    switch (StringParser::stringFindStringID(params, propName)) {
                                        case 0:
                                            //filepath
                                            job->driver.filepath = propValue;
                                            break;
                                        case 1:
                                            //filename
                                            job->driver.filename = propValue;
                                            break;
                                        case 3:
                                            //props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->driver.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 4:
                                            //int-props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->driver.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            job->driver.driverSetPlugin(job,job->driver.filepath,
                                                                job->driver.filename,
                                                                job->driver.fp64_props,
                                                                job->driver.int_props);
                            std::cout << "Driver Configured." << std::endl;
                        }
                        break;
                    case 3:
                        //solver
                        params = {"filepath","filename","properties","int-properties"};
                        line = StringParser::stringRemoveComments(line);
                        line = StringParser::stringRemoveSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = StringParser::stringRemoveBraces(line);
                                lvec = StringParser::stringSplitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = StringParser::stringSplitString(propValue, ',');
                                    switch (StringParser::stringFindStringID(params, propName)) {
                                        case 0:
                                            //filepath
                                            job->solver.filepath = propValue;
                                            break;
                                        case 1:
                                            //filename
                                            job->solver.filename = propValue;
                                            break;
                                        case 3:
                                            //props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->solver.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 4:
                                            //int-props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->solver.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            job->solver.solverSetPlugin(job,job->solver.filepath,
                                                        job->solver.filename,
                                                        job->solver.fp64_props,
                                                        job->solver.int_props);
                            std::cout << "Solver Configured." << std::endl;
                        }
                        break;
                    case 4:
                        //body
                        job->bodies.push_back(Body());
                        job->activeBodies.push_back(0); //set body to inactive
                        id = job->bodies.size() - 1;
                        job->bodies[id].id = id;

                        params = {"name", "point-file",
                                  "material-filepath","material-filename",
                                  "material-props","material-int-props",
                                  "boundary-filepath","boundary-filename",
                                  "boundary-props","boundary-int-props"};
                        line = StringParser::stringRemoveComments(line);
                        line = StringParser::stringRemoveSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = StringParser::stringRemoveBraces(line);
                                lvec = StringParser::stringSplitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = StringParser::stringSplitString(propValue, ',');
                                    switch (StringParser::stringFindStringID(params, propName)) {
                                        case 0:
                                            //name
                                            job->bodies[id].name = propValue;
                                            break;
                                        case 1:
                                            //point-file
                                            job->bodies[id].points.file = propValue;
                                            break;
                                        case 2:
                                            //material-filepath
                                            job->bodies[id].material.filepath = propValue;
                                            break;
                                        case 3:
                                            //material-filename
                                            job->bodies[id].material.filename = propValue;
                                            break;
                                        case 4:
                                            //material-props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->bodies[id].material.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 5:
                                            //material-int-props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->bodies[id].material.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 6:
                                            //boundary-filepath
                                            job->bodies[id].boundary.filepath = propValue;
                                            break;
                                        case 7:
                                            //boundary-filename
                                            job->bodies[id].boundary.filename = propValue;
                                            break;
                                        case 8:
                                            //boundary-props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->bodies[id].boundary.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 9:
                                            //boundary-int-props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->bodies[id].boundary.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            if (!job->bodies[id].points.file.empty()){
                                job->bodies[id].points.pointsReadFromFile(job,
                                                                          &(job->bodies[id]),
                                                                          job->bodies[id].points.file);
                                job->activeBodies[id] = 1;
                            } else {
                                std::cerr << "No point file for body: " << job->bodies[id].name << "!" << std::endl;
                            }
                            if (!job->bodies[id].boundary.filename.empty()) {
                                job->bodies[id].boundary.boundarySetPlugin(job, &(job->bodies[id]),
                                                                           job->bodies[id].boundary.filepath,
                                                                           job->bodies[id].boundary.filename,
                                                                           job->bodies[id].boundary.fp64_props,
                                                                           job->bodies[id].boundary.int_props);
                                job->bodies[id].activeBoundary = 1;
                            }
                            if (!job->bodies[id].material.filename.empty()) {
                                job->bodies[id].material.materialSetPlugin(job, &(job->bodies[id]),
                                                                           job->bodies[id].material.filepath,
                                                                           job->bodies[id].material.filename,
                                                                           job->bodies[id].material.fp64_props,
                                                                           job->bodies[id].material.int_props);
                                job->bodies[id].activeMaterial = 1;
                            }

                            std::cout << "Body [" << id << ", " << job->bodies[id].name << "] Configured." << std::endl;
                        }
                        break;
                    case 5:
                        //contact
                        job->contacts.push_back(Contact());
                        job->activeContacts.push_back(0); //set body to inactive
                        id = job->contacts.size() - 1;
                        job->contacts[id].id = id;

                        params = {"filepath","filename","properties","int-properties"};
                        line = StringParser::stringRemoveComments(line);
                        line = StringParser::stringRemoveSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = StringParser::stringRemoveBraces(line);
                                lvec = StringParser::stringSplitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = StringParser::stringSplitString(propValue, ',');
                                    switch (StringParser::stringFindStringID(params, propName)) {
                                        case 0:
                                            //filepath
                                            job->contacts[id].filepath = propValue;
                                            break;
                                        case 1:
                                            //filename
                                            job->contacts[id].filename = propValue;
                                            break;
                                        case 3:
                                            //props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->contacts[id].fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 4:
                                            //int-props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->contacts[id].int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            job->contacts[id].contactSetPlugin(job,job->contacts[id].filepath,
                                                        job->contacts[id].filename,
                                                        job->contacts[id].fp64_props,
                                                        job->contacts[id].int_props);
                            job->activeContacts[id] = 1;
                            std::cout << "Contact [" << id << ", " << job->contacts[id].name << "] Configured." << std::endl;
                        }
                        break;
                    case 6:
                        //grid
                        params = {"filepath","filename","properties","int-properties"};
                        line = StringParser::stringRemoveComments(line);
                        line = StringParser::stringRemoveSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = StringParser::stringRemoveComments(line);
                                line = StringParser::stringRemoveSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = StringParser::stringRemoveBraces(line);
                                lvec = StringParser::stringSplitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = StringParser::stringSplitString(propValue, ',');
                                    switch (StringParser::stringFindStringID(params, propName)) {
                                        case 0:
                                            //filepath
                                            job->grid.filepath = propValue;
                                            break;
                                        case 1:
                                            //filename
                                            job->grid.filename = propValue;
                                            break;
                                        case 3:
                                            //props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->grid.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 4:
                                            //int-props
                                            for (size_t i = 0; i < lvec.size(); i++) {
                                                job->grid.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            job->grid.gridSetPlugin(job,job->grid.filepath,
                                                        job->grid.filename,
                                                        job->grid.fp64_props,
                                                        job->grid.int_props);

                            std::cout << "Grid Configured." << std::endl;
                        }
                        break;
                    default:
                        std::cerr << "Header \"" << line << "\" not recognized." << std::endl;
                }
            }
        }
        //close file
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << file << std::endl;
        return 0;
    }

    return 1;
}
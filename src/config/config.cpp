//
// Created by aaron on 5/17/18.
// config.cpp
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "parser.hpp"
#include "job.hpp"
#include "config.hpp"

#include "mpm_objects.hpp"
#include "registry.hpp"

#include <Eigen/Core>
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

/*----------------------------------------------------------------------------*/
//
void Configurator::init(std::string fileIN){
    serializer_registry = Registry<Serializer>();
    driver_registry = Registry<Driver>();
    solver_registry = Registry<Solver>();

    grid_registry = Registry<Grid>();
    body_registry = Registry<Body>();
    points_registry = Registry<Points>();
    nodes_registry = Registry<Nodes>();
    contact_registry = Registry<Contact>();

    material_registry = Registry<Material>();
    boundary_registry = Registry<Boundary>();

    file = fileIN;
    return;
}


/*----------------------------------------------------------------------------*/
//
void Configurator::setMainPath(std::string program){
    std::vector<std::string> svec;
    svec = Parser::splitString(program,'/');
    std::string filepath = "";
    for (int i=0; i<(svec.size()-1);i++){
        filepath += svec[i];
        filepath += "/";
    }
    mainpath = filepath;
    return;
}

/*----------------------------------------------------------------------------*/
//
void Configurator::checkConfigFile(std::string filename) {
    std::vector<std::string> headers = {"job","serializer","driver","solver","grid","body","contact"};
    std::vector<int> check = {0,0,0,0,0,0,0};
    int id;
    std::string line;
    std::ifstream fin(filename);
    if (fin.is_open()){
        while (std::getline(fin,line)){
            line = Parser::removeComments(line);
            line = Parser::removeSpaces(line);
            id = Parser::findStringID(headers,line);
            if (id != -1){
                check[id] = 1;
            }
        }
        for (int i=0;i<headers.size();i++){
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


/*----------------------------------------------------------------------------*/
//
int Configurator::configureJob(Job* job){
    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::vector<std::string> params;
    std::string propName;
    std::string propValue;
    int id; //for bodies and contacts
    std::string filename, name; //for points and body


    //section headers
    std::vector<std::string> headers = {"job","serializer","driver","solver","grid","body","contact"};

    //declare fstream from filename
    std::ifstream fin(file);

    //read config file and configure
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin,line)){
            //remove spaces and comments from line
            line = Parser::removeComments(line);
            line = Parser::removeSpaces(line);

            //check if line matches header
            if (line.size()>0){
                //temporary objects for constructing job objects
                MPMObject tmp, point_tmp, node_tmp, material_tmp, boundary_tmp;

                //switch for headers
                switch (Parser::findStringID(headers,line)) {
                    case 0:
                        //job
                        params = {"dt","t","TYPE"};
                        std::getline(fin,line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    //end of block
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //dt
                                            job->dt = std::stod(propValue);
                                            break;
                                        case 1:
                                            //t
                                            job->t = std::stod(propValue);
                                            break;
                                        case 2:
                                            //JOB_TYPE
                                            job->assignJobType(std::stoi(propValue));
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
                        params = {"class","properties","int-properties","str-properties"};
                        std::getline(fin,line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            job->serializer = serializer_registry.get_object(propValue);
                                            break;
                                        case 1:
                                            //props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->serializer->fp64_props.push_back(std::stod(lvec[i]));
                                                tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 2:
                                            //int-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->serializer->int_props.push_back(std::stoi(lvec[i]));
                                                tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 3:
                                            //str-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->serializer->str_props.push_back(lvec[i]);
                                                tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            if (job->serializer){
                                //if serializer is valid
                                job->serializer->fp64_props = tmp.fp64_props;
                                job->serializer->int_props = tmp.int_props;
                                job->serializer->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: Serializer object needs valid \'class\' defined. Exiting." << std::endl;
                                exit(0);
                            }
                            std::cout << "Serializer Configured: " << job->serializer->object_name << std::endl;
                        }
                        break;
                    case 2:
                        //driver
                        params = {"class","properties","int-properties","str-properties"};
                        std::getline(fin,line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            job->driver = driver_registry.get_object(propValue);
                                            break;
                                        case 1:
                                            //props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->driver->fp64_props.push_back(std::stod(lvec[i]));
                                                tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 2:
                                            //int-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->driver->int_props.push_back(std::stoi(lvec[i]));
                                                tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 3:
                                            //str-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->driver->str_props.push_back(lvec[i]);
                                                tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            if (job->driver){
                                //if driver is valid
                                job->driver->fp64_props = tmp.fp64_props;
                                job->driver->int_props = tmp.int_props;
                                job->driver->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: Driver object needs valid \'class\' defined. Exiting." << std::endl;
                                exit(0);
                            }
                            std::cout << "Driver Configured: " << job->driver->object_name  << std::endl;
                        }
                        break;
                    case 3:
                        //solver
                        params = {"class","properties","int-properties","str-properties"};
                        std::getline(fin,line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            job->solver = solver_registry.get_object(propValue);
                                            break;
                                        case 1:
                                            //props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->solver->fp64_props.push_back(std::stod(lvec[i]));
                                                tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 2:
                                            //int-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->solver->int_props.push_back(std::stoi(lvec[i]));
                                                tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 3:
                                            //str-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->solver->str_props.push_back(lvec[i]);
                                                tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            if (job->solver){
                                //if solver is valid
                                job->solver->fp64_props = tmp.fp64_props;
                                job->solver->int_props = tmp.int_props;
                                job->solver->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: Solver object needs valid \'class\' defined. Exiting." << std::endl;
                                exit(0);
                            }
                            std::cout << "Solver Configured: " << job->solver->object_name << std::endl;
                        }
                        break;
                    case 4:
                        //grid
                        params = {"class","properties","int-properties","str-properties"};
                        std::getline(fin,line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            job->grid = grid_registry.get_object(propValue);
                                            break;
                                        case 1:
                                            //props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->grid->fp64_props.push_back(std::stod(lvec[i]));
                                                tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 2:
                                            //int-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->grid->int_props.push_back(std::stoi(lvec[i]));
                                                tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 3:
                                            //str-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->grid->str_props.push_back(lvec[i]);
                                                tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            if (job->grid){
                                job->grid->fp64_props = tmp.fp64_props;
                                job->grid->int_props = tmp.int_props;
                                job->grid->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: Grid object needs valid \'class\' defined. Exiting." << std::endl;
                                exit(0);
                            }
                            std::cout << "Grid Configured: " << job->grid->object_name << std::endl;
                        }
                        break;
                    case 5:
                        //body
                        job->bodies.push_back(std::unique_ptr<Body>(nullptr)); //
                        job->activeBodies.push_back(0); //set body to inactive
                        id = job->bodies.size() - 1;
                        //job->bodies[id].id = id;

                        params = {"class", "name", "point-file", "properties", "int-properties", "str-properties",
                                  "point-class", "point-props", "point-int-props", "point-str-props",
                                  "node-class", "node-props", "node-int-props", "node-str-props",
                                  "material-class", "material-props", "material-int-props", "material-str-props",
                                  "boundary-class", "boundary-props", "boundary-int-props", "boundary-str-props"};
                        std::getline(fin,line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            job->bodies[id] = body_registry.get_object(propValue);
                                            job->bodies[id]->id = id;
                                        case 1:
                                            //name
                                            //job->bodies[id]->name = propValue;
                                            name = propValue;
                                            break;
                                        case 2:
                                            //point-file
                                            //job->bodies[id]->points->file = propValue;
                                            filename = propValue;
                                            break;
                                        case 3:
                                            //properties
                                            for (int i = 0; i < lvec.size(); i++){
                                                tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 4:
                                            //int-properties
                                            for (int i = 0; i < lvec.size(); i++){
                                                tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 5:
                                            //str-properties
                                            for (int i = 0; i < lvec.size(); i++){
                                                tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        case 6:
                                            //point-class
                                            //job->bodies[id]->points = points_registry.get_object(propValue);
                                            point_tmp.object_name = propValue;
                                            break;
                                        case 7:
                                            //point-props
                                            for (int i=0;i<lvec.size();i++){
                                                point_tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 8:
                                            //point-int-props
                                            for (int i=0;i<lvec.size();i++){
                                                point_tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 9:
                                            //point-str-props
                                            for (int i=0;i<lvec.size();i++){
                                                point_tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        case 10:
                                            //node-class
                                            node_tmp.object_name = propValue;
                                            break;
                                        case 11:
                                            //node-props
                                            for (int i=0;i<lvec.size();i++){
                                                node_tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 12:
                                            //node-int-props
                                            for (int i=0;i<lvec.size();i++){
                                                node_tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 13:
                                            //node-str-props
                                            for (int i=0;i<lvec.size();i++){
                                                node_tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        case 14:
                                            //material-class
                                            //job->bodies[id]->material = material_registry.get_object(propValue);
                                            material_tmp.object_name = propValue;
                                            break;
                                        case 15:
                                            //material-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->bodies[id]->material->fp64_props.push_back(std::stod(lvec[i]));
                                                material_tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 16:
                                            //material-int-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->bodies[id]->material->int_props.push_back(std::stoi(lvec[i]));
                                                material_tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 17:
                                            //material-str-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->bodies[id]->material->str_props.push_back(lvec[i]);
                                                material_tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        case 18:
                                            //boundary-class
                                            //job->bodies[id]->boundary = boundary_registry.get_object(propValue);
                                            boundary_tmp.object_name = propValue;
                                            break;
                                        case 19:
                                            //boundary-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->bodies[id]->boundary->fp64_props.push_back(std::stod(lvec[i]));
                                                boundary_tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 20:
                                            //boundary-int-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->bodies[id]->boundary->int_props.push_back(std::stoi(lvec[i]));
                                                boundary_tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 21:
                                            //boundary-str-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->bodies[id]->boundary->str_props.push_back(lvec[i]);
                                                boundary_tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            if (job->bodies[id]){
                                //if body is valid
                                //job->activeBodies[id] = 1;
                                job->bodies[id]->name = name;
                                job->bodies[id]->fp64_props = tmp.fp64_props;
                                job->bodies[id]->int_props = tmp.int_props;
                                job->bodies[id]->str_props = tmp.str_props;

                                if (!point_tmp.object_name.empty()){
                                    //if a point object is defined
                                    job->bodies[id]->points = points_registry.get_object(point_tmp.object_name);
                                    if (job->bodies[id]->points){
                                        job->bodies[id]->points->file = filename;
                                        job->bodies[id]->points->fp64_props = point_tmp.fp64_props;
                                        job->bodies[id]->points->int_props = point_tmp.int_props;
                                        job->bodies[id]->points->str_props = point_tmp.str_props;

                                        //setup job and points
                                        job->bodies[id]->points->readFromFile(job, job->bodies[id].get(), filename);
                                        job->activeBodies[id] = 1;
                                    } else {
                                        //if object is not created, this is fatal
                                        std::cerr << "ERROR: Points object needs valid \'class\' defined. Exiting." << std::endl;
                                        exit(0);
                                    }
                                } else {
                                    std::cerr << "No Points!" << std::endl;
                                }

                                if (!node_tmp.object_name.empty()){
                                    //if a point object is defined
                                    job->bodies[id]->nodes = nodes_registry.get_object(node_tmp.object_name);
                                    if (job->bodies[id]->nodes){
                                        job->bodies[id]->nodes->fp64_props = node_tmp.fp64_props;
                                        job->bodies[id]->nodes->int_props = node_tmp.int_props;
                                        job->bodies[id]->nodes->str_props = node_tmp.str_props;
                                    } else {
                                        //if object is not created, this is fatal
                                        std::cerr << "ERROR: Nodes object needs valid \'class\' defined. Exiting." << std::endl;
                                        exit(0);
                                    }
                                } else {
                                    std::cerr << "No Nodes!" << std::endl;
                                }

                                if (!material_tmp.object_name.empty()){
                                    //if a point object is defined
                                    job->bodies[id]->material = material_registry.get_object(material_tmp.object_name);
                                    if (job->bodies[id]->material){
                                        job->bodies[id]->material->fp64_props = material_tmp.fp64_props;
                                        job->bodies[id]->material->int_props = material_tmp.int_props;
                                        job->bodies[id]->material->str_props = material_tmp.str_props;

                                        job->bodies[id]->activeMaterial = 1;
                                    } else {
                                        //if object is not created, this is fatal
                                        std::cerr << "ERROR: Material object needs valid \'class\' defined. Exiting." << std::endl;
                                        exit(0);
                                    }
                                } else {
                                    std::cerr << "No Material!" << std::endl;
                                }

                                if (!boundary_tmp.object_name.empty()){
                                    //if a point object is defined
                                    job->bodies[id]->boundary = boundary_registry.get_object(boundary_tmp.object_name);
                                    if (job->bodies[id]->boundary){
                                        job->bodies[id]->boundary->fp64_props = boundary_tmp.fp64_props;
                                        job->bodies[id]->boundary->int_props = boundary_tmp.int_props;
                                        job->bodies[id]->boundary->str_props = boundary_tmp.str_props;

                                        job->bodies[id]->activeBoundary = 1;
                                    } else {
                                        //if object is not created, this is fatal
                                        std::cerr << "ERROR: Boundary object needs valid \'class\' defined. Exiting." << std::endl;
                                        exit(0);
                                    }
                                } else {
                                    std::cerr << "No Points!" << std::endl;
                                }
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: Body object needs valid \'class\' defined. Exiting." << std::endl;
                                exit(0);
                            }

                            std::cout << "Body [" << id << ", " << job->bodies[id]->name << "] Configured: " << job->bodies[id]->object_name << std::endl;
                        }
                        break;
                    case 6:
                        //contact
                        job->contacts.push_back(std::unique_ptr<Contact>(nullptr));
                        job->activeContacts.push_back(0); //set body to inactive
                        id = job->contacts.size() - 1;
                        //job->contacts[id]->id = id;

                        params = {"name","class","properties","int-properties","str-properties"};
                        std::getline(fin,line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin,line)){
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line,'=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //name
                                            //job->contacts[id]->name = propValue;
                                            name = propValue;
                                            break;
                                        case 1:
                                            //class
                                            job->contacts[id] = contact_registry.get_object(propValue);
                                            job->contacts[id]->id = id;
                                            break;
                                        case 2:
                                            //props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->contacts[id]->fp64_props.push_back(std::stod(lvec[i]));
                                                tmp.fp64_props.push_back(std::stod(lvec[i]));
                                            }
                                            break;
                                        case 3:
                                            //int-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->contacts[id]->int_props.push_back(std::stoi(lvec[i]));
                                                tmp.int_props.push_back(std::stoi(lvec[i]));
                                            }
                                            break;
                                        case 4:
                                            //str-props
                                            for (int i = 0; i < lvec.size(); i++) {
                                                //job->contacts[id]->str_props.push_back(lvec[i]);
                                                tmp.str_props.push_back(lvec[i]);
                                            }
                                            break;
                                        default:
                                            std::cerr << "Parameter \"" << propName << "\" not recognized." << std::endl;
                                    }
                                }
                            }
                            if (job->contacts[id]){
                                //if body is valid
                                job->activeContacts[id] = 1;
                                job->contacts[id]->name = name;
                                job->contacts[id]->fp64_props = tmp.fp64_props;
                                job->contacts[id]->int_props = tmp.int_props;
                                job->contacts[id]->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: Contact object needs valid \'class\' defined. Exiting." << std::endl;
                                exit(0);
                            }
                            std::cout << "Contact [" << id << ", " << job->contacts[id]->name << "] Configured: " << job->contacts[id]->object_name << std::endl;
                        }
                        break;
                    default:
                        std::cerr << "Header \"" << line << "\" not recognized." << std::endl;
                }
            }
        }
        //close file
        fin.close();

        std::cout << "Job Configured.\n" << std::endl;
    } else {
        std::cout << "ERROR: Unable to open file: " << file << std::endl;
        return 0;
    }

    return 1;
}
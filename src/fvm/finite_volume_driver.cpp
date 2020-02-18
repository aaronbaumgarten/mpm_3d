//
// Created by aaron on 12/12/19.
// finite_volume_driver.cpp
//

//
// Created by aaron on 5/10/18.
// default_driver.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <time.h>

#include "job.hpp"
#include "mpm_objects.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"

#include "fvm_objects.hpp"
#include "registry.hpp"

/*----------------------------------------------------------------------------*/

void FiniteVolumeDriver::checkConfigFile(std::string filename) {
    std::vector<std::string> headers = {"serializer", "grid", "body", "material", "solver"};
    std::vector<int> check = {0,0,0,0,0};
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

void FiniteVolumeDriver::init(Job* job){
    int GRID_DIM;
    //assign grid dimension from job type
    if (job->JOB_TYPE == job->JOB_1D){
        GRID_DIM = 1;
    } else if (job->JOB_TYPE == job->JOB_2D){
        GRID_DIM = 2;
    } else if (job->JOB_TYPE == job->JOB_3D){
        GRID_DIM = 3;
    } else if (job->JOB_TYPE == job->JOB_2D_OOP){
        GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
    } else if (job->JOB_TYPE == job->JOB_AXISYM){
        GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
    } else {
        std::cerr << "Job doesn't have defined type for input " << job->JOB_TYPE << "." << std::endl;
    }

    //check size of properties passed to driver object
    if (fp64_props.size() < 1 || str_props.size() < 1) {
        std::cout << fp64_props.size() << ", " << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined (stop_time, file_name).\n",
                __FILE__, __func__);
        exit(0);
    } else if (fp64_props.size() < 1+GRID_DIM) {

        //store stop_time
        stop_time = fp64_props[0];

        //let gravity be zeroo
        gravity = KinematicVector(job->JOB_TYPE);

        //store filename
        file = str_props[0];

    } else {
        stop_time = fp64_props[0];
        gravity = KinematicVector(job->JOB_TYPE);
        for (int pos=0; pos<GRID_DIM; pos++){
            gravity[pos] = fp64_props[pos+1];
        }
        file = str_props[0];
    }

    //define ORDER of element-wise flux reconstructions
    if (int_props.size() >= 1){
        ORDER = int_props[0];
    }

    //ensure ORDER is valid
    if (ORDER < 1){
        std::cout << "ERROR: Order selected is less than 1, that doesn't seem right." << std::endl;
        ORDER = 1;
    }

    //set simulation type, i.e. underlying fluid assumptions
    if (int_props.size() >= 2){
        TYPE = int_props[1];
    }

    //ensure TYPE is valid
    if (TYPE != THERMAL && TYPE != ISOTHERMAL && TYPE != INCOMPRESSIBLE){
        std::cout << "ERROR: Simulation type " << TYPE << " is not defined. Exiting." << std::endl;
        exit(0);
    }

    //print grid properties
    std::cout << "Driver properties (stop_time = " << stop_time << ", gravity = < ";
    for (int pos=0; pos<GRID_DIM; pos++){
        std::cout << gravity[pos] << " ";
    }
    std::cout << ">, file = " << file << ", ORDER = " << ORDER << ")." << std::endl;

    /*---------------------------------------*/
    //read in finite volume configuration file
    /*---------------------------------------*/

    //headers
    std::vector<std::string> headers = {"serializer", "grid", "body", "material", "solver"};

    //declare variables and strings and vectors
    std::string line;
    std::vector<std::string> lvec;
    std::vector<std::string> params;
    std::string propName;
    std::string propValue;

    //declare fstream from filename
    std::ifstream fin(file);

    //registries
    Registry<FiniteVolumeSerializer> serializer_registry = Registry<FiniteVolumeSerializer>();
    Registry<FiniteVolumeGrid> grid_registry = Registry<FiniteVolumeGrid>();
    Registry<FiniteVolumeBody> body_registry = Registry<FiniteVolumeBody>();
    Registry<FiniteVolumeMaterial> material_registry = Registry<FiniteVolumeMaterial>();
    Registry<FiniteVolumeSolver> solver_registry = Registry<FiniteVolumeSolver>();

    //read config file and configure
    if (fin.is_open()) {
        //if open, read lines
        while (std::getline(fin, line)) {
            //remove spaces and comments from line
            line = Parser::removeComments(line);
            line = Parser::removeSpaces(line);

            //check if line matches header
            if (line.size() > 0) {
                //temporary objects for constructing job objects
                MPMObject tmp;

                line = Parser::removeBraces(line);
                line = Parser::removeQuotes(line);
                lvec = Parser::splitString(line, '=');

                //switch for headers
                switch (Parser::findStringID(headers, lvec[0])) {
                    case 0:
                        //serializer
                        params = {"class", "properties", "int-properties", "str-properties"};
                        std::getline(fin, line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin, line)) {
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line, '=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            serializer = serializer_registry.get_object(propValue);
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
                                            std::cerr << "Parameter \"" << propName << "\" not recognized."
                                                      << std::endl;
                                    }
                                }
                            }
                            if (serializer) {
                                //if serializer is valid
                                serializer->fp64_props = tmp.fp64_props;
                                serializer->int_props = tmp.int_props;
                                serializer->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: FiniteVolumeSerializer object needs valid \'class\' defined. Exiting."
                                          << std::endl;
                                exit(0);
                            }
                            std::cout << "FiniteVolumeSerializer Configured: " << serializer->object_name << std::endl;
                        }
                        break;
                    case 1:
                        //grid
                        params = {"class", "properties", "int-properties", "str-properties"};
                        std::getline(fin, line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin, line)) {
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line, '=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            fluid_grid = grid_registry.get_object(propValue);
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
                                            std::cerr << "Parameter \"" << propName << "\" not recognized."
                                                      << std::endl;
                                    }
                                }
                            }
                            if (fluid_grid) {
                                //if serializer is valid
                                fluid_grid->fp64_props = tmp.fp64_props;
                                fluid_grid->int_props = tmp.int_props;
                                fluid_grid->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: FiniteVolumeGrid object needs valid \'class\' defined. Exiting."
                                          << std::endl;
                                exit(0);
                            }
                            std::cout << "FiniteVolumeGrid Configured: " << fluid_grid->object_name << std::endl;
                        }
                        break;
                    case 2:
                        //body
                        params = {"class", "properties", "int-properties", "str-properties"};
                        std::getline(fin, line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin, line)) {
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line, '=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            fluid_body = body_registry.get_object(propValue);
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
                                            std::cerr << "Parameter \"" << propName << "\" not recognized."
                                                      << std::endl;
                                    }
                                }
                            }
                            if (fluid_body) {
                                //if serializer is valid
                                fluid_body->fp64_props = tmp.fp64_props;
                                fluid_body->int_props = tmp.int_props;
                                fluid_body->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: FiniteVolumeBody object needs valid \'class\' defined. Exiting."
                                          << std::endl;
                                exit(0);
                            }
                            std::cout << "FiniteVolumeBody Configured: " << fluid_body->object_name << std::endl;
                        }
                        break;
                    case 3:
                        //material
                        params = {"class", "properties", "int-properties", "str-properties"};
                        std::getline(fin, line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin, line)) {
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line, '=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            fluid_material = material_registry.get_object(propValue);
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
                                            std::cerr << "Parameter \"" << propName << "\" not recognized."
                                                      << std::endl;
                                    }
                                }
                            }
                            if (fluid_material) {
                                //if serializer is valid
                                fluid_material->fp64_props = tmp.fp64_props;
                                fluid_material->int_props = tmp.int_props;
                                fluid_material->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: FiniteVolumeMaterial object needs valid \'class\' defined. Exiting."
                                          << std::endl;
                                exit(0);
                            }
                            std::cout << "FiniteVolumeMaterial Configured: " << fluid_material->object_name << std::endl;
                        }
                        break;
                    case 4:
                        //solver
                        params = {"class", "properties", "int-properties", "str-properties"};
                        std::getline(fin, line);
                        line = Parser::removeComments(line);
                        line = Parser::removeSpaces(line);
                        if (line.compare("{") == 0) {
                            while (std::getline(fin, line)) {
                                line = Parser::removeComments(line);
                                line = Parser::removeSpaces(line);

                                if (line.compare("}") == 0) {
                                    break;
                                }

                                line = Parser::removeBraces(line);
                                line = Parser::removeQuotes(line);
                                lvec = Parser::splitString(line, '=');
                                if (lvec.size() > 1) {
                                    propName = lvec[0];
                                    propValue = lvec[1];

                                    lvec = Parser::splitString(propValue, ',');
                                    switch (Parser::findStringID(params, propName)) {
                                        case 0:
                                            //class
                                            solver = solver_registry.get_object(propValue);
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
                                            std::cerr << "Parameter \"" << propName << "\" not recognized."
                                                      << std::endl;
                                    }
                                }
                            }
                            if (solver) {
                                //if serializer is valid
                                solver->fp64_props = tmp.fp64_props;
                                solver->int_props = tmp.int_props;
                                solver->str_props = tmp.str_props;
                            } else {
                                //if object is not created, this is fatal
                                std::cerr << "ERROR: FiniteVolumeSolver object needs valid \'class\' defined. Exiting."
                                          << std::endl;
                                exit(0);
                            }
                            std::cout << "FiniteVolumeSolver Configured: " << solver->object_name << std::endl;
                        }
                        break;
                }
            }
        }
    }

    //close file
    fin.close();

    //wait to initialize FVM objects until after run() is called (ensures MPM is properly configured)
    std::cout << "Driver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void FiniteVolumeDriver::run(Job* job) {
    //initialize FVM objects
    serializer->init(job, this);
    fluid_grid->init(job, this);
    fluid_body->init(job, this);
    fluid_material->init(job, this);
    solver->init(job, this);

    //set counters to zero
    int stepCount = 0;
    int frameCount = 0;

    struct timespec timeStart, timeFrame, timeFinish;
    clock_gettime(CLOCK_MONOTONIC, &timeStart);
    timeFrame = timeStart;
    //clock_t clockSim = clock();
    //clock_t clockFrame = clock();
    double tSim = 0;
    double tFrame = 0;

    //initialize gravity
    generateGravity(job);
    applyGravity(job);

    //run simulation until stop_time
    while (job->t <= stop_time){
        //run solver
        solver->step(job,this);

        //std::cout << "Step Completed [" << ++stepCount << "]." << std::flush;
        if (job->serializer->writeFrame(job) == 1) {
            //call fvm serializer to write frame as well:
            serializer->writeFrame(job,this);

            //successful frame written
            //tFrame = (double)(clock() - clockFrame)/CLOCKS_PER_SEC;
            //tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
            clock_gettime(CLOCK_MONOTONIC,&timeFinish);
            tFrame = (timeFinish.tv_sec - timeFrame.tv_sec) + (timeFinish.tv_nsec - timeFrame.tv_nsec)/1000000000.0;
            tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
            timeFrame = timeFinish;
            printf("\33[2K");
            std::cout << "Frame Written [" << ++frameCount << "]. Time/Frame [" << tFrame << " s]. Elapsed Time [" << tSim << " s]." << std::flush;
        }
        std::cout << "\r";
        job->t += job->dt;
    }
    //tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
    clock_gettime(CLOCK_MONOTONIC,&timeFinish);
    tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
    std::cout << std::endl << std::endl << "Simulation Complete. Elapsed Time [" << tSim << "s]." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/

void FiniteVolumeDriver::generateGravity(Job* job) {
    //nothing to do here
    return;
}

/*----------------------------------------------------------------------------*/

void FiniteVolumeDriver::applyGravity(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (int i=0;i<job->bodies[b]->points->b.size();i++){
            job->bodies[b]->points->b(i) = gravity;
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/

std::string FiniteVolumeDriver::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "error";
}

/*----------------------------------------------------------------------------*/

int FiniteVolumeDriver::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

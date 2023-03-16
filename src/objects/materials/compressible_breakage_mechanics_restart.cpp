//
// Created by aaron on 5/15/18.
// isolin.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <dlfcn.h>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "materials.hpp"

/*----------------------------------------------------------------------------*/
//initialize assuming that general properties have been assigned correctly
//fp64_props etc. have been filled by configuration object
void CompressibleBreakageMechanicsRestart::init(Job* job, Body* body){

    // 0. Open Restart File

    if (str_props.size() < 1){
        std::cerr << "CompressibleBreakageMechanicsRestart requires 1 input file in order to restart. " << str_props.size() << " given. Exiting." << std::endl;
        exit(0);
    } else {
        point_file = str_props[0];
        std::cout << "CompressibleBreakageMechanicsRestart restarting simulation using data from:\n";
        std::cout << "    Points: " << point_file << "\n";
    }

    // 1. Call Parent Initializer
    CompressibleBreakageMechanicsSand::init(job,body);

    // 2. Read in Restart Data
    // Open Point File
    std::ifstream pin(point_file);
    if (!pin.is_open()){
        std::cerr << "ERROR: Unable to open point restart file: " << point_file << ". Exiting." << std::endl;
        exit(0);
    }

    // Write Update to Console
    std::cout << "Restarting CompressibleBreakageMechanicsSand Material." << std::endl;

    // Read Point Data
    std::string line;
    std::vector<std::string> svec;
    int result = -1;
    int plen = 0;
    if (pin.is_open()){
        // Look For "POINTS" Label
        while (std::getline(pin,line)){
            //Split Line Using " "
            svec = Parser::splitString(line,' ');

            //Check for POINTS Label
            if (svec.size() > 0 && svec[0].compare("POINTS") == 0){
                //POINTS Label Found!
                //Next Value is Number of Points
                plen = std::stod(svec[1]);

                //Exit While Loop
                break;
            }
        }

        // Read in Data
        if (plen == body->points->x.size()){

            // Write Update to Console
            std::cout << "   Reading Point Data." << std::endl;

            // Read Through Lines Until "POINT_DATA"
            while (std::getline(pin,line)){
                //Split Line Using " "
                svec = Parser::splitString(line,' ');

                //Check for POINTS Label
                if (svec.size() > 0 && svec[0].compare("POINT_DATA") == 0){
                    //POINT_DATA Label Found!
                    //Exit While Loop
                    break;
                }
            }

            // Write Update to Console
            std::cout << "    Found POINT_DATA.\n";

            // Keywords
            std::vector<std::string> scalar_tags = {"B","phiS","phiP","Es"};
            std::vector<std::string> vector_tags = {};
            std::vector<std::string> tensor_tags = {"F","Be"};

            // Read Through Lines Until SCALARS, VECTORS, or TENSORS Found
            while (std::getline(pin,line)){
                //Split Line Using " "
                svec = Parser::splitString(line,' ');

                //Check for Data Labels
                if (svec.size() > 0 && svec[0].compare("SCALARS") == 0){
                    //Found SCALARS Label
                    //Next Value is Data Name
                    result = -1;
                    result = Parser::findStringID(scalar_tags, svec[1]);

                    //Read Next Line
                    std::getline(pin,line);

                    //Read Data Into Correct Array
                    for (int i=0; i<plen; i++) {
                        //Read Data Line
                        std::getline(pin,line);

                        //SCALARS Only Have 1 Value
                        switch (result) {
                            case 0:
                                //B
                                B(i) = std::stod(line);
                                break;
                            case 1:
                                //phiS
                                phiS(i) = std::stod(line);
                                break;
                            case 2:
                                //phiP
                                phiP(i) = std::stoi(line);
                                break;
                            case 3:
                                //Es
                                Es(i) = std::stoi(line);
                                break;
                            default:
                                //do nothing
                                break;
                        }
                    }

                    if (result >= 0 && result < scalar_tags.size()) {
                        //Write Update to Console
                        std::cout << "    Read In " << scalar_tags[result] << " Data...\n";
                    }

                } else if (svec.size() > 0 && svec[0].compare("VECTORS") == 0){
                    //Found VECTORS Label
                    //Next Value is Data Name
                    result = -1;
                    result = Parser::findStringID(vector_tags, svec[1]);

                    //Read Data Into Correct Array
                    for (int i=0; i<plen; i++) {
                        //Read Data Line
                        std::getline(pin,line);

                        //VECTORS Have 3 Values
                        svec = Parser::splitString(line,' ');
                        for (int j=0; j<3; j++) {
                            switch (result) {
                                default:
                                    //do nothing
                                    break;
                            }
                        }
                    }

                    if (result >= 0 && result < vector_tags.size()) {
                        //Write Update to Console
                        std::cout << "    Read In " << vector_tags[result] << " Data...\n";
                    }

                } else if (svec.size() > 0 && svec[0].compare("TENSORS") == 0){
                    //Found TENSORS Label
                    //Next Value is Data Name
                    result = -1;
                    result = Parser::findStringID(tensor_tags, svec[1]);

                    //Read Data Into Correct Array
                    for (int i=0; i<plen; i++) {
                        //Read Data Line
                        std::getline(pin,line);

                        //TENSORS Have 3 Sets of 3 Values
                        svec = Parser::splitString(line,' ');
                        for (int j=0; j<3; j++) {
                            switch (result) {
                                case 0:
                                    //F
                                    F(i,j) = std::stod(svec[j]);
                                    break;
                                case 1:
                                    //Be
                                    Be(i,j) = std::stod(svec[j]);
                                    break;
                                default:
                                    //do nothing
                                    break;
                            }
                        }
                        std::getline(pin,line);
                        svec = Parser::splitString(line,' ');
                        for (int j=0; j<3; j++) {
                            switch (result) {
                                case 0:
                                    //F
                                    F(i,j+3) = std::stod(svec[j]);
                                    break;
                                case 1:
                                    //Be
                                    Be(i,j+3) = std::stod(svec[j]);
                                    break;
                                default:
                                    //do nothing
                                    break;
                            }
                        }
                        std::getline(pin,line);
                        svec = Parser::splitString(line,' ');
                        for (int j=0; j<3; j++) {
                            switch (result) {
                                case 0:
                                    //F
                                    F(i,j+6) = std::stod(svec[j]);
                                    break;
                                case 1:
                                    //Be
                                    Be(i,j+6) = std::stod(svec[j]);
                                    break;
                                default:
                                    //do nothing
                                    break;
                            }
                        }

                        //TENSORS Also Has Additional Line Break
                        std::getline(pin,line);
                    }

                    if (result >= 0 && result < tensor_tags.size()) {
                        //Write Update to Console
                        std::cout << "    Read In " << tensor_tags[result] << " Data...\n";
                    }
                }
            }

        } else {
            std::cerr << "ERROR: Unable to find POINTS label in " << point_file << ". Exiting." << std::endl;
            exit(0);
        }

        //Write Update to Console
        std::cout << "    Complete.\n";

        //Close File
        pin.close();
    }

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

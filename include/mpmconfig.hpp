//
// Created by aaron on 12/26/16.
// mpmconfig.hpp
//

#ifndef MPM_3D_MPMCONFIG_HPP
#define MPM_3D_MPMCONFIG_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include "process.hpp"
#include "mpmio.hpp"

class MPMconfig{
public:
    //config filename
    std::string configFile;
    std::string mainPath;

    //input parameters
    std::vector<std::string> jobParams;
    std::vector<std::string> inputParams;
    std::vector<std::string> materialParams;
    std::vector<std::string> boundaryParams;
    std::vector<std::string> contactParams;
    std::vector<std::string> outputParams;

    //constructors
    MPMconfig();
    MPMconfig(std::string);

    //functions
    void setConfigFile(std::string);
    void setMainPath(std::string);
    int checkConfigFile(std::string);
    std::string removeSpaces(std::string);
    std::string removeComments(std::string);
    std::string removeBraces(std::string);
    std::vector<std::string> splitString(std::string,char);
    int findStringID(std::vector<std::string>,std::string);
    int configJob(job_t*);
    int configInput(job_t*);
    int configMaterial(job_t*);
    int configBoundary(job_t*);
    int configContact(job_t*);
    int configOutput(job_t*,MPMio*);
};

#endif //MPM_3D_MPMCONFIG_HPP

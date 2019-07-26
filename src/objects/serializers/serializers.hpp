//
// Created by aaron on 5/10/18.
// serializers.hpp
//

#ifndef MPM_V3_SERIALIZERS_HPP
#define MPM_V3_SERIALIZERS_HPP

#include "mpm_objects.hpp"
#include <iostream>
#include <fstream>
#include <Eigen/Core>
#include <iomanip>
#include <sys/stat.h>

/*
 * IN THIS FILE, DEFINE SERIALIZER OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
class Serializer : public MPMObject{
public:
    //path to main program
    std::string mainpath;

    //set mainpath
    void setMainPath(std::string program){
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

    //functions which must be implemented by every serializer
    virtual void init(Job*) = 0;                    //initialize from Job object
    virtual std::string saveState(Job*) = 0;        //save entire job to file
    virtual int loadState(Job*, std::string) = 0;   //load from file
    virtual int writeFrame(Job*) = 0;               //write frame for job object (return 1 on success)
    virtual void writeDefaultPointHeader(Job*, Body*, std::ofstream&, int) = 0; //write default header to file
    virtual void writeDefaultNodeHeader(Job*, Body*, std::ofstream&, int) = 0; //write default header to file
    virtual void writeScalarArray(Eigen::VectorXd&, std::string) = 0;   //write scalar array to file (named by string)
    virtual void writeVectorArray(MPMVectorArray&, std::string) = 0;    //write vector array to file
    virtual void writeTensorArray(MPMTensorArray&, std::string) = 0;    //write tensor array to file
};
 */


/*----------------------------------------------------------------------------*/
//default serializer for vtk files
class DefaultVTK : public Serializer {
public:
    //default constructor
    DefaultVTK(){
        object_name = "DefaultVTK"; //set this here
        mainpath = "";              //for now, set to nothing
    }

    //file storing and writing variables
    std::string frameDirectory, outputDirectory, outputName;

    //job writing and sampling
    int sampledFrames;
    double sampleRate, t_last_frame;

    //frame writing variables
    std::string pfilename, nfilename;
    std::ofstream pfile, nfile;
    int plen, nlen;
    Body* currentBody;

    virtual void init(Job* job);
    virtual std::string saveState(Job* job);
    virtual int loadState(Job* job, std::string fullpath);
    virtual int writeFrame(Job*);

    virtual void writeDefaultPointHeader(Job* job, Body* body, std::ofstream& pfile, int SPEC);
    virtual void writeDefaultNodeHeader(Job* job, Body* body, std::ofstream& nfile, int SPEC);
    virtual void writeScalarArray(Eigen::VectorXd& scalarArray, std::string scalarName);
    virtual void writeVectorArray(MPMVectorArray& vectorArray, std::string vectorName);
    virtual void writeTensorArray(MPMTensorArray& tensorArray, std::string tensorName);
};

class MinimalVTK : public DefaultVTK {
public:
    //default constructor
    MinimalVTK(){
        object_name = "MinimalVTK";
        mainpath = "";
    }

    int writeFrame(Job*);
};

class SlicePointsVTK : public DefaultVTK {
public:
    //default constructor
    SlicePointsVTK(){
        object_name = "SlicePointsVTK";
        mainpath = "";
    }

    int stride = 1;
    void init(Job*);
    int writeFrame(Job*);
};

#endif //MPM_V3_SERIALIZERS_HPP

//
// Created by aaron on 12/23/19.
// fvm_serializers.hpp
//

#ifndef MPM_V3_FVM_SERIALIZERS_HPP
#define MPM_V3_FVM_SERIALIZERS_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <fstream>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "mpm_objects.hpp"
#include "fvm_objects.hpp"

/*
//responsible for writing frames
class FiniteVolumeSerializer : public MPMObject{
public:
    static const int VTK = 1; //in case i ever make another type of serializer

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
    virtual void init(Job*, FiniteVolumeDriver*) = 0;                   //initialize from Job object
    virtual int writeFrame(Job*, FiniteVolumeDriver*) = 0;              //write frame for job object (return 1 on success)
    virtual void writeScalarArray(Eigen::VectorXd&, std::string) = 0;   //write scalar array to file (named by string)
    virtual void writeVectorArray(MPMVectorArray&, std::string) = 0;    //write vector array to file
    virtual void writeTensorArray(MPMTensorArray&, std::string) = 0;    //write tensor array to file
};
 */

class FVMDefaultVTK: public FiniteVolumeSerializer {
public:
    FVMDefaultVTK(){
        object_name = "FVMDefaultVTK";
        mainpath = "";
    }

    //frame counter
    int sampledFrames = 0;

    //file storing and writing variables
    std::string frameDirectory, outputName, filename;
    std::ofstream file;

    //functions which must be implemented by every serializer
    virtual void init(Job* job, FiniteVolumeDriver* driver);                   //initialize from Job object
    virtual int writeFrame(Job* job, FiniteVolumeDriver* driver);              //write frame for job object (return 1 on success)
    virtual void writeScalarArray(Eigen::VectorXd& scalarArray, std::string scalarName);   //write scalar array to file (named by string)
    virtual void writeVectorArray(MPMVectorArray& vectorArray, std::string vectorName);    //write vector array to file
    virtual void writeTensorArray(MPMTensorArray& tensorArray, std::string tensorName);    //write tensor array to file
};

#endif //MPM_V3_FVM_SERIALIZERS_HPP

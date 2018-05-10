//
// Created by aaron on 5/9/18.
// mpm_objects.hpp
//

#ifndef MPM_V3_MPMOBJECT_HPP
#define MPM_V3_MPMOBJECT_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include "parser.hpp"
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

class Job; //forward declare job object
class Body; //forward declare body object
class Points; //forward declare points object
class Nodes; //forward declare nodes object
class Boundary; //forward declare boundary
class Material; //forward declare material

/*------------------------------------------------------------------------*/
//base class for all registry objects
class MPMObject{
public:
    std::vector<double> fp64_props;     //double properties
    std::vector<int> int_props;         //integer properties
    std::vector<std::string> str_props; //string properties

    virtual ~MPMObject() = default;
};

/*------------------------------------------------------------------------*/
//serializer class
//responsible for writing frames and reading/writing save files
class Serializer : public MPMObject{
public:
    //path to main program
    std::string mainpath;

    //set mainpath
    void setMainPath(std::string program){
        std::vector<std::string> svec;
        svec = Parser::splitString(program,'/');
        std::string filepath = "";
        for (size_t i=0; i<(svec.size()-1);i++){
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
    virtual void writeScalarArray(Eigen::VectorXd&, std::string) = 0;   //write scalar array to file (named by string)
    virtual void writeVectorArray(MPMVectorArray&, std::string) = 0;    //write vector array to file
    virtual void writeTensorArray(MPMTensorArray&, std::string) = 0;    //write tensor array to file
};


/*------------------------------------------------------------------------*/
//driver class
//responsible for running mpm problem (essentially the 'while loop')
class Driver : public MPMObject{
public:
    //functions which must be implemented by every driver
    virtual void init(Job*) = 0;                                        //initialize from Job
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save to file (in given directory) and return filename
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load from file
    virtual void run(Job*) = 0;                                         //run mpm according to problem
    virtual void generateGravity(Job*) = 0;                             //generate gravity
    virtual void applyGravity(Job*) = 0;                                //apply gravity
};


/*----------------------------------------------------------------------------*/
//solver class
//responsible for one forward mpm time-step (everything in the 'while loop')
class Solver : public MPMObject{
public:
    //functions that must be implemented by every solver (not many)
    virtual void init(Job*) = 0;                                        //initialize from Job
    virtual void step(Job*) = 0;                                        //perform single mpm step
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save to file (in given directory) and return filename
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load from file
};

/*----------------------------------------------------------------------------*/
//body class
//represents a single continuum body
class Body : public MPMObject{
public:
    int id;             //numerical id of body in simulation
    std::string name;   //name of body in simulation
    int activeMaterial; //is the continuum body material defined?
    int activeBoundary; //is the continuum body boundary defined?

    std::unique_ptr<Points> points;
    std::unique_ptr<Nodes> nodes;
    std::unique_ptr<Material> material;
    std::unique_ptr<Boundary> boundary;

    virtual void init(Job*) = 0;                                        //initialiaze from Job
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save to file (in given directory)
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load data from full path
};


/*----------------------------------------------------------------------------*/
//point class
//carries material point representation of continuum fields
class Points : public MPMObject{
public:
    std::string file;

    KinematicVectorArray x, u, x_t, mx_t, b;    //state vectors
    KinematicTensorArray L;                     //velocity gradient
    MaterialTensorArray T;                      //cauchy stress
    Eigen::VectorXd m, v, v0;                   //weight measures
    Eigen::VectorXi active;                     //active?

    virtual void init(Job*, Body*) = 0;                         //initialize from Job and Body
    virtual void readFromFile(Job*, Body*, std::string) = 0;    //construct points from given file

    virtual void writeFrame(Job*, Body*, Serializer*) = 0;                      //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;   //save to file (in given directory)
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;           //load data from full path
};

/*----------------------------------------------------------------------------*/
//node class
//carries grid node representation of continuum fields
class Nodes : public MPMObject{
public:
    KinematicVectorArray x, u, x_t, a, mx_t, f; //state vectors
    Eigen::VectorXd m, V;                       //weight measures
    Eigen::VectorXi active;                     //active?

    virtual void init(Job*, Body*) = 0;                                         //initialize from Job and Body
    virtual void writeFrame(Job*, Body*, Serializer*) = 0;                      //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;   //save to file (in given directory)
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;           //load from full path
};

/*----------------------------------------------------------------------------*/
//material class
//responsible for constitutive update on continuum body
class Material : public MPMObject{
public:
    static const int UPDATE = 1;    //update internal variables
    static const int TRIAL = 0;     //don't update internal variables

    virtual void init(Job*, Body*) = 0;                                     //initialize from Job and Body
    virtual void calculateStress(Job*, Body*, int) = 0;                     //calculate stress update
    virtual void assignStress(Job*, Body*, MaterialTensor&, int, int) = 0;  //assign stress to given id
    virtual void assignPressure(Job*, Body*, double, int, int) = 0;         //assign pressure to given id

    virtual void writeFrame(Job*, Body*, Serializer*) = 0;          //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*) = 0;    //save state to file
    virtual int loadState(Job*, Body*, Serializer*) = 0;            //load from full path
};

/*----------------------------------------------------------------------------*/
//boundary file
//responsible for boundary conditions on continuum body
class Boundary : public MPMObject{
    virtual void init(Job*, Body*) = 0;             //initialize from Job and Body
    virtual void generateRules(Job*, Body*) = 0;    //generate boundary rules
    virtual void applyRules(Job*, Body*) = 0;       //apply boundary rules

    virtual void writeFrame(Job*, Body*, Serializer*) = 0;                      //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;   //save state to file
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;           //load from fullpath
};

class Contact; class Grid;


#endif //MPM_V3_MPMOBJECT_HPP

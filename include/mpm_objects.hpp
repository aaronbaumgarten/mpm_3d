//
// Created by aaron on 5/9/18.
// mpm_objects.hpp
//

#ifndef MPM_V3_MPMOBJECT_HPP
#define MPM_V3_MPMOBJECT_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <fstream>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

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
    std::string object_name;            //name of object (in registry)
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
    int activeMaterial = 0; //is the continuum body material defined?
    int activeBoundary = 0; //is the continuum body boundary defined?

    MPMScalarSparseMatrix S = MPMScalarSparseMatrix(0,0); //S_ip maps ith node to pth point
    KinematicVectorSparseMatrix gradS = KinematicVectorSparseMatrix(0,0); //gradS_ip maps ith node gradient to pth point

    std::unique_ptr<Points> points;
    std::unique_ptr<Nodes> nodes;
    std::unique_ptr<Material> material;
    std::unique_ptr<Boundary> boundary;

    virtual void init(Job*) = 0;                                        //initialiaze from Job
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save to file (in given directory)
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load data from full path
    virtual void generateMap(Job*, int) = 0;                            //generate S and gradS

    virtual void generateLoads(Job*) = 0;   //arbitrary loading during simulation
    virtual void applyLoads(Job*) = 0;
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
    virtual void generateMap(Job*, Body*, int) = 0;             //generate S and gradS
    virtual void updateIntegrators(Job*, Body*) = 0;            //update integrators (extent, etc.)

    virtual void writeHeader(Job*, Body*, Serializer*, std::ofstream&, int) = 0;
    virtual void writeFrame(Job*, Body*, Serializer*) = 0;                      //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;   //save to file (in given directory)
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;           //load data from full path

    virtual void generateLoads(Job*, Body*) = 0;   //arbitrary loading during simulation
    virtual void applyLoads(Job*, Body*) = 0;
};

/*----------------------------------------------------------------------------*/
//node class
//carries grid node representation of continuum fields
class Nodes : public MPMObject{
public:
    KinematicVectorArray x, u, x_t, a, mx_t, f, diff_x_t; //state vectors
    Eigen::VectorXd m, V;                       //weight measures
    Eigen::VectorXi active;                     //active?

    virtual void init(Job*, Body*) = 0;                                         //initialize from Job and Body
    virtual void writeFrame(Job*, Body*, Serializer*) = 0;                      //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;   //save to file (in given directory)
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;           //load from full path

    virtual void generateLoads(Job*, Body*) = 0;   //arbitrary loading during simulation
    virtual void applyLoads(Job*, Body*) = 0;
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
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;    //save state to file
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;            //load from full path
};

/*----------------------------------------------------------------------------*/
//boundary file
//responsible for boundary conditions on continuum body
class Boundary : public MPMObject{
public:
    virtual void init(Job*, Body*) = 0;             //initialize from Job and Body
    virtual void generateRules(Job*, Body*) = 0;    //generate boundary rules
    virtual void applyRules(Job*, Body*) = 0;       //apply boundary rules

    virtual void writeFrame(Job*, Body*, Serializer*) = 0;                      //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;   //save state to file
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;           //load from fullpath
};


/*----------------------------------------------------------------------------*/
//contact class
//responsible for inter-body interactions
class Contact : public MPMObject{
public:
    static const int IMPLICIT = 1; //implicit contact rules
    static const int EXPLICIT = 0; //explicit contact rules

    int id;             //numerical id of contact in simulation
    std::string name;   //name of contact in simulation

    virtual void init(Job*) = 0;            //initialize from Job
    virtual void generateRules(Job*) = 0;   //generate contact rules
    virtual void applyRules(Job*, int) = 0; //apply contact rules

    virtual void writeFrame(Job*, Serializer*) = 0;                     //send frame data to Serializer
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save state to file
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load from file
};

/*----------------------------------------------------------------------------*/
//grid class
//responsible for the definition of background grid
class Grid : public MPMObject{
public:
    int node_count;      //number of nodes which define grid
    int element_count;   //number of elements in grid

    int GRID_DIM = -1; //dimension of grid (might be different than simulation dimension)

    virtual void init(Job*) = 0; //initialize from Job

    virtual void writeHeader(Job*, Body*, Serializer*, std::ofstream&, int) = 0; //write cell types
    virtual void writeFrame(Job*, Serializer*) = 0;                     //send frame data to Serializer
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save to file
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load from file

    virtual int whichElement(Job*, KinematicVector&) = 0;                                                       //return element given position
    virtual bool inDomain(Job*, KinematicVector&) = 0;                                                          //check if location is in domain
    virtual KinematicVector nodeIDToPosition(Job*, int) = 0;                                                    //return position of given node
    virtual void evaluateBasisFnValue(Job*, KinematicVector&, std::vector<int>&, std::vector<double>&) = 0;     //fill id and value for nodal weights of position
    virtual void evaluateBasisFnGradient(Job*, KinematicVector&, std::vector<int>&, KinematicVectorArray&) = 0; //fill id and gradient for nodal gradient of position
    virtual double nodeVolume(Job*, int) = 0;                                                                   //return node volume of id
    virtual double elementVolume(Job*, int) = 0;                                                                //return element volume of id
    virtual int nodeTag(Job*, int) = 0;                                                                         //return 'tag' of node id
};


#endif //MPM_V3_MPMOBJECT_HPP

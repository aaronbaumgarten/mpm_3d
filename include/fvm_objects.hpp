//
// Created by aaron on 12/12/19.
// fvm_objects.hpp
//

#ifndef MPM_V3_FVM_OBJECTS_HPP
#define MPM_V3_FVM_OBJECTS_HPP

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

class FiniteVolumeDriver;

/*----------------------------------------------------------------------------*/
//finite solver class
//responsible for one forward time-step
class FiniteVolumeSolver : public MPMObject{
public:
    //functions that must be implemented by every finite volume solver (not many)
    virtual void init(Job*, FiniteVolumeDriver*) = 0;                                        //initialize from Job
    virtual void step(Job*, FiniteVolumeDriver*) = 0;                                        //perform single mpm step
};

/*------------------------------------------------------------------------*/
//finite volume grid class
//responsible for tracking element, face definitions and cross integrals with mpm grid
class FiniteVolumeGrid: public MPMObject{
public:

    //boundary condition types
    static const int DIRICHLET = 0;
    static const int NEUMANN = 1;
    static const int PERIODIC = 2;
    static const int NEUMANN_DAMPING = 3;
    static const int SUPERSONIC_INLET = 4;

    //Harten entropy correction scale
    static constexpr double delta = 0.1;

    //grid definions
    int face_count;      //number of faces which define grid
    int node_count;      //number of nodes which define grid
    int element_count;   //number of elements in grid

    int GRID_DIM = -1; //dimension of grid (might be different than simulation dimension)

    MPMScalarSparseMatrix M = MPMScalarSparseMatrix(0,0); //M_ji maps jth element to ith mpm node
    KinematicVectorSparseMatrix gradM = KinematicVectorSparseMatrix(0,0); //gradM_ji maps jth element to ith node gradient

    virtual void init(Job*, FiniteVolumeDriver*) = 0; //initialize from mpm Job

    virtual void writeHeader(std::ofstream&, int) = 0; //write cell types

    virtual double getElementVolume(int) = 0;      //return element volume of id
    virtual int getElementTag(int) = 0;            //return 'tag' of id
    virtual double getFaceArea(int) = 0;           //return face area of id
    virtual int getFaceTag(int) = 0;               //return 'tag of face id
    virtual KinematicVector getFaceNormal(Job*, int) = 0; //return normal vector associated with given id

    virtual std::vector<int> getElementFaces(int) = 0; //return faces associated with given element id
    virtual std::array<int,2> getOrientedElementsByFace(int) = 0; //return elements associated with given face id A -|-> B
    virtual std::vector<int> getElementNeighbors(int) = 0; //return list of elements in neighborhood
    virtual KinematicVector getElementCentroid(Job*, int) = 0; //return element centroid
    virtual KinematicVector getFaceCentroid(Job*, int) = 0;

    virtual void generateMappings(Job*, FiniteVolumeDriver*) = 0; //generate M_ji and gradM_ji maps
    virtual void constructMomentumField(Job*, FiniteVolumeDriver*) = 0; //construct momentum from FV body
    virtual void constructDensityField(Job*, FiniteVolumeDriver*) = 0; //construct momentum from FV body
    virtual KinematicTensorArray getVelocityGradients(Job*, FiniteVolumeDriver*) = 0; //return velocity gradient in each element

    //functions to compute element-wise fluxes of field variables using reconstructed velocity field
    virtual Eigen::VectorXd calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, Eigen::VectorXd&) = 0;
    //virtual KinematicVectorArray calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, KinematicVectorArray&, int) = 0;
    //virtual KinematicTensorArray calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, KinematicTensorArray&, int) = 0;
    //virtual MaterialVectorArray calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, MaterialVectorArray&, int) = 0;
    //virtual MaterialTensorArray calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, MaterialTensorArray&, int) = 0;

    //functions to compute element mass flux
    virtual Eigen::VectorXd calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver) = 0;

    //functions to compute element momentum fluxes
    virtual KinematicVectorArray calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver) = 0;
};


/*----------------------------------------------------------------------------*/
//finite volume body class
//represents a single continuum body filling fluid domain
class FiniteVolumeBody : public MPMObject{
public:
    //initialize from job and driver
    virtual void init(Job*, FiniteVolumeDriver*) = 0;

    //container for fluid fields
    KinematicTensorArray p_x;      //momentum gradient
    KinematicVectorArray p, rho_x; //momentum and density gradient
    Eigen::VectorXd rho, P, theta; //density, pressure, and temperature
    MaterialTensorArray tau;       //shear stress
};

/*----------------------------------------------------------------------------*/
//finite volume material class
//represents a single continuum body filling fluid domain
class FiniteVolumeMaterial : public MPMObject{
public:
    //initialize from job and driver
    virtual void init(Job*, FiniteVolumeDriver*) = 0;

    //functions to calculate fluid fields
    virtual MaterialTensor getStress(Job*, FiniteVolumeDriver*, KinematicVector x, KinematicTensor L, double rho, double theta) = 0;
    virtual MaterialTensor getShearStress(Job*, FiniteVolumeDriver*, KinematicVector x, KinematicTensor L, double rho, double theta) = 0;
    virtual double getPressure(Job*, FiniteVolumeDriver*, KinematicVector x, double rho, double theta) = 0;
    virtual double getSpeedOfSound(Job*, FiniteVolumeDriver*, KinematicVector x, double rho, double theta) = 0;
    virtual void calculateElementPressures(Job*, FiniteVolumeDriver*) = 0;
    virtual void calculateElementShearStresses(Job*, FiniteVolumeDriver*) = 0;
    virtual void solveMaterialEquations(Job*, FiniteVolumeDriver*) = 0;
};

/*------------------------------------------------------------------------*/
//serializer class
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

/*------------------------------------------------------------------------*/
//driver class
//responsible for running mpm-fvm problem
class FiniteVolumeDriver: public Driver {
public:
    FiniteVolumeDriver(){
        object_name = "FiniteVolumeDriver"; //set object name here
    }

    int GRID_DIM = 3;

    //functions which must be implemented by every driver
    virtual void init(Job*);                                        //initialize from Job
    virtual std::string saveState(Job*, Serializer*, std::string);  //save to file (in given directory) and return filename
    virtual int loadState(Job*, Serializer*, std::string);          //load from file
    virtual void run(Job*);                                         //run mpm according to problem
    virtual void generateGravity(Job*);                             //generate gravity
    virtual void applyGravity(Job*);                                //apply gravity

    //objects for running finite volume method
    std::unique_ptr<FiniteVolumeSolver> solver;
    std::unique_ptr<FiniteVolumeGrid> fluid_grid;
    std::unique_ptr<FiniteVolumeBody> fluid_body;
    std::unique_ptr<FiniteVolumeMaterial> fluid_material;
    std::unique_ptr<FiniteVolumeSerializer> serializer;

    //function to check input file
    virtual void checkConfigFile(std::string);

    //internal data structures
    double stop_time;
    KinematicVector gravity;
    std::string file;

    //order of finite volume reconstruction
    int order = 2;
};

/*------------------------------------------------------------------------*/

#endif //MPM_V3_FVM_OBJECTS_HPP

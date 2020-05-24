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
    virtual void writeFrame(Job* job, FiniteVolumeDriver* driver){
        //do nothing
        return;
    }
};

/*------------------------------------------------------------------------*/
//finite volume grid class
//responsible for tracking element, face definitions and cross integrals with mpm grid
class FiniteVolumeGrid: public MPMObject{
public:

    //boundary conditions
    static const int VELOCITY_INLET         = 0;    //dirichlet
    static const int VELOCITY_TEMP_INLET    = 1;    //dirichlet
    static const int VELOCITY_DENSITY_INLET = 2;    //dirichlet
    static const int PRESSURE_INLET         = 3;
    static const int PRESSURE_OUTLET        = 4;
    static const int DAMPED_OUTLET          = 5;
    static const int ADIABATIC_WALL         = 6;    //dirichlet
    static const int THERMAL_WALL           = 7;    //dirichlet
    static const int SYMMETRIC_WALL         = 8;    //dirichlet
    static const int SUPERSONIC_INLET       = 9;    //dirichlet
    static const int SUPERSONIC_OUTLET      = 10;
    static const int PERIODIC               = 11;
    static const int DAMPED_WALL            = 12;
    static const int STAGNATION_INLET       = 13;

    //Harten entropy correction scale
    static constexpr double delta = 0.1;

    //grid definions
    int face_count;      //number of faces which define grid
    int node_count;      //number of nodes which define grid
    int element_count;   //number of elements in grid
    int int_quad_count;  //number of interior quadrature points
    int ext_quad_count;  //number of boundary/face quadrature points
    int qpe, qpf;        //number of quadrature points per element/face

    int GRID_DIM = -1; //dimension of grid (might be different than simulation dimension)

    MPMScalarSparseMatrix M = MPMScalarSparseMatrix(0,0); //M_ji maps jth element to ith mpm node
    MPMScalarSparseMatrix Q = MPMScalarSparseMatrix(0,0); //Q_ji maps ith mpm node to jth quadrature point
    KinematicVectorSparseMatrix gradQ = KinematicVectorSparseMatrix(0,0); //gradQ_ji maps ith node gradient to jth quadrature point

    virtual void init(Job*, FiniteVolumeDriver*) = 0; //initialize from mpm Job

    virtual void writeHeader(std::ofstream&, int) = 0; //write cell types

    virtual double getElementVolume(int) = 0;      //return element volume of id
    virtual int getElementTag(int) = 0;            //return 'tag' of id
    virtual double getFaceArea(int) = 0;           //return face area of id
    virtual int getFaceTag(int) = 0;               //return 'tag of face id
    virtual KinematicVector getFaceNormal(Job*, int) = 0; //return normal vector associated with given id

    virtual std::vector<int> getElementQuadraturePoints(int e) = 0; //return vector with global quadrature indices
    virtual std::vector<int> getFaceQuadraturePoints(int f) = 0;    //return vector with global quadrature indices

    virtual std::vector<int> getElementFaces(int) = 0; //return faces associated with given element id
    virtual std::array<int,2> getOrientedElementsByFace(int) = 0; //return elements associated with given face id A -|-> B
    virtual std::vector<int> getElementNeighbors(int) = 0; //return list of elements in neighborhood
    virtual KinematicVector getElementCentroid(Job*, int) = 0; //return element centroid
    virtual KinematicVector getFaceCentroid(Job*, int) = 0;
    virtual KinematicVector getQuadraturePosition(Job* job, int q) = 0;
    virtual double getQuadratureWeight(int q) = 0;

    virtual void generateMappings(Job*, FiniteVolumeDriver*) = 0; //generate M, Q, gradQ maps
    virtual void mapMixturePropertiesToQuadraturePoints(Job*, FiniteVolumeDriver*) = 0; //map porosity and solid velocity to quadrature points
    virtual void constructMomentumField(Job*, FiniteVolumeDriver*) = 0; //construct momentum from FV body
    virtual void constructDensityField(Job*, FiniteVolumeDriver*) = 0; //construct momentum from FV body
    virtual void constructEnergyField(Job*, FiniteVolumeDriver*) = 0; //construct energy field from FV body
    virtual void constructPorosityField(Job*, FiniteVolumeDriver*) = 0; //construct porosity field from FV body
    virtual KinematicTensorArray getVelocityGradients(Job*, FiniteVolumeDriver*) = 0; //return velocity gradient in each element

    //functions to compute element mass flux
    virtual Eigen::VectorXd calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver) = 0;

    //functions to compute element momentum fluxes
    virtual KinematicVectorArray calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver) = 0;

    //function to compute element energy fluxes
    virtual Eigen::VectorXd calculateElementEnergyFluxes(Job* job, FiniteVolumeDriver* driver) = 0;

    //function to calculate interphase force
    virtual KinematicVectorArray calculateInterphaseForces(Job* job, FiniteVolumeDriver* driver) = 0;
    virtual KinematicVectorArray calculateBuoyantForces(Job* job, FiniteVolumeDriver* driver) = 0;
    virtual KinematicVectorArray calculateDragForces(Job* job, FiniteVolumeDriver* driver) = 0;

    virtual KinematicVectorArray calculateCorrectedDragForces(Job *job,
                                                              FiniteVolumeDriver *driver,
                                                              const Eigen::VectorXd &K_n) = 0; // <- drag coefficient at each quad point

    virtual Eigen::VectorXd getCorrectedDragCoefficients(Job* job, FiniteVolumeDriver* driver) = 0; // <- function to collect those coefficients

    //interphase force functions when f_d_e != M*f_d_i
    virtual void calculateSplitIntegralInterphaseForces(Job* job,
                                                        FiniteVolumeDriver* driver,
                                                        KinematicVectorArray& f_i,
                                                        KinematicVectorArray& f_e) = 0;

    virtual void calculateSplitIntegralBuoyantForces(Job* job,
                                                     FiniteVolumeDriver* driver,
                                                     KinematicVectorArray& f_i,
                                                     KinematicVectorArray& f_e) = 0;

    virtual void calculateSplitIntegralDragForces(Job* job,
                                                  FiniteVolumeDriver* driver,
                                                  KinematicVectorArray& f_i,
                                                  KinematicVectorArray& f_e) = 0;

    virtual void calculateSplitIntegralCorrectedDragForces(Job* job,
                                                            FiniteVolumeDriver* driver,
                                                            KinematicVectorArray& f_i,
                                                            KinematicVectorArray& f_e,
                                                            const Eigen::VectorXd &K_n) = 0;

    //functions to calculate interphase energy terms
    virtual Eigen::VectorXd calculateInterphaseEnergyFlux(Job* job, FiniteVolumeDriver* driver) = 0;
    virtual Eigen::VectorXd calculateInterphaseEnergyFluxUsingNodeBasedDrag(Job* job,
                                                                            FiniteVolumeDriver* driver,
                                                                            const KinematicVectorArray &f_d) = 0;
    virtual Eigen::VectorXd calculateInterphaseEnergyFluxUsingElementBasedForce(Job* job,
                                                                                FiniteVolumeDriver* driver,
                                                                                const KinematicVectorArray &f_e) = 0;
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
    KinematicVectorArray p, rho_x, rhoE_x; //momentum and density gradient
    Eigen::VectorXd rho, rhoE, P, theta; //density, pressure, and temperature
    MaterialTensorArray tau;       //shear stress

    //container for solid phase fields in mixture
    KinematicVectorArray v_s;   //solid phase velocity field (defined on MPM grid)
    Eigen::VectorXd n;          //mixture porosity field (defined on MPM grid)

    Eigen::VectorXd n_e;                 //mixture porosity field (defined on finite volumes)
    KinematicVectorArray n_e_x;        //mixture porosity gradient (defined on finite volumes)
};

/*----------------------------------------------------------------------------*/
//finite volume material class
//represents a single continuum body filling fluid domain
class FiniteVolumeMaterial : public MPMObject{
public:
    static const int REGULAR_DRAG = 0;
    static const int CORRECTED_DRAG = 1;

    //initialize from job and driver
    virtual void init(Job*, FiniteVolumeDriver*) = 0;

    //functions to calculate fluid fields
    virtual MaterialTensor getStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n) = 0;
    virtual MaterialTensor getShearStress(Job*, FiniteVolumeDriver*, const KinematicTensor& L, double rho, const KinematicVector& p, double rhoE, double n) = 0;
    virtual double getPressure(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n) = 0;
    virtual double getTemperature(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n) = 0;
    virtual double getSpeedOfSound(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& p, double rhoE, double n) = 0;
    virtual void calculateElementPressures(Job*, FiniteVolumeDriver*) = 0;
    virtual void calculateElementShearStresses(Job*, FiniteVolumeDriver*) = 0;
    virtual void calculateElementTemperatures(Job*, FiniteVolumeDriver*) = 0;

    //fluid equations of state
    virtual double getDensityFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n) = 0;
    virtual double getInternalEnergyFromPressureAndTemperature(Job*, FiniteVolumeDriver*, double pressure, double theta, double n) = 0;
    virtual double getPressureFromDensityAndTemperature(Job*, FiniteVolumeDriver*, double rho, double theta, double n) = 0;
    virtual KinematicVector getHeatFlux(Job*, FiniteVolumeDriver*, double rho, double theta, const KinematicVector& theta_x, double n) = 0;
    virtual double getSpeedOfSoundFromEnthalpy(Job *, FiniteVolumeDriver *, double rho, const KinematicVector &p, double rhoH, double n) = 0;
    virtual double getPressureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Pt, double M) = 0;
    virtual double getTemperatureFromStagnationProperties(Job*, FiniteVolumeDriver*, double Tt, double M) = 0;

    //mixture model functions
    virtual int calculatePorosity(Job*, FiniteVolumeDriver*) = 0; //return 1 if mixture problem, return 0 if FVM problem only
    virtual int updateSolidPhaseVelocity(Job*, FiniteVolumeDriver*) = 0; //return 1 if mixture problem, return 0 if FVM problem only
    virtual KinematicVector getInterphaseDrag(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n) = 0;
    virtual double getInterphaseDragCoefficient(Job*, FiniteVolumeDriver*, double rho, const KinematicVector& v_f, const KinematicVector& v_s, double n, int SPEC = REGULAR_DRAG) = 0;
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

    //ORDER of finite volume reconstruction
    int ORDER = 2;

    //TYPE of finite volume method
    static const int THERMAL = 0;
    static const int ISOTHERMAL = 1;
    static const int INCOMPRESSIBLE = 2;
    int TYPE = 1;

};

/*------------------------------------------------------------------------*/

#endif //MPM_V3_FVM_OBJECTS_HPP

//
// Created by aaron on 12/23/19.
// fvm_bodies.hpp
//

#ifndef MPM_V3_FVM_DRIVERS_HPP
#define MPM_V3_FVM_DRIVERS_HPP

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

#include "mpm_objects.hpp"
#include "fvm_objects.hpp"
#include "fvm_artificial_viscosity.hpp"

/*
 * class FiniteVolumeDriver: public Driver {
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

    //new function for returning body forces at specified position in space
    virtual KinematicVector getFluidLoading(Job *job, const KinematicVector &x);
    virtual KinematicVector getSolidLoading(Job *job, const KinematicVector &x);

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
 */

class FVMVariableStepDriver : public FiniteVolumeDriver{
public:
    FVMVariableStepDriver(){
        object_name = "FVMVariableStepDriver"; //set object name here
    }

    bool USE_ARTIFICIAL_VISCOSITY = false;
    double dt0; //baseline time increment
    double lambda; //element quality (1.0 = perfect scaling, 0.0 = terrible scaling)
    double eta; //fluid viscosity

    //initialize from job
    virtual void init(Job* job);
    virtual void run(Job* job);

    //parallel functions
    virtual void applyArtificialViscosityFluxes(Job* job, std::vector<ArtificialViscosityCalculator::fluxVector> &fluxVectors);
    static void calculateSubsetOfAVFluxes(Job* job,
                                               FiniteVolumeDriver* driver,
                                               std::vector<ArtificialViscosityCalculator::fluxVector> &fluxVectors,
                                               int f_begin, int f_end,
                                               volatile bool &done);
};

class FVMNumericalDampingDriver : public FiniteVolumeDriver{
public:
    FVMNumericalDampingDriver(){
        object_name = "FVMNumericalDampingDriver"; //set object name here
    }

    double t0 = 1.0; //damping decay rate
    double C0 = 0.0; //damping scale

    bool AVOID_FLUID_DAMPING = false;
    bool LOAD_SOLID_DIRECTLY = false;

    KinematicVector xmax, load;

    //initialize from job
    virtual void init(Job* job);
    virtual void run(Job* job);
    virtual KinematicVector getSolidLoading(Job *job, const KinematicVector &x);

};

#endif //MPM_V3_FVM_DRIVERS_HPP

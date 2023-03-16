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

    bool USE_THERMAL_FLOOR = false;
    bool USE_ARTIFICIAL_VISCOSITY = false;
    double dt0; //baseline time increment
    double lambda; //element quality (1.0 = perfect scaling, 0.0 = terrible scaling)
    double eta; //fluid viscosity
    double T_floor; //thermal floor
    double cv_floor; //heat capacity for thermal floor

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

class FVMBolideImpactDriver : public FVMVariableStepDriver{
public:
    FVMBolideImpactDriver(){
        FVMVariableStepDriver::object_name = "FVMBolideImpactDriver";
        object_name = "FVMBolideImpactDriver"; //set object name here
    }

    //static variables
    static double H;// = 0;        //Height in Meters
    static double V;// = 11e3;     //Velocity in Meters/Second
    static double Y0;// = 0;       //Initial Leading Position of Projectile in Meters

    //constant variables
    static constexpr int bmax = 7;
    static constexpr double Hb[8] = {0, 11e3, 20e3, 32e3, 47e3, 51e3, 71e3, 85e3};                     //Heights of Model in Meters
    static constexpr double Lb[8] = {-6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3, 0.0};        //Gradients of Model in Kelvin per Meter
    static constexpr double Tb[8] = {288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.65};  //Thermal Layers of Model in Kelvin
    static constexpr double Pb[8] = {101325, 22632, 5474.9, 802.135, 102.488, 61.8581, 3.6561, 0.3358}; //Pressure Payers of Model in Pascal
    static constexpr double g0 = 9.80665;      //Gravitational Acceleration in Meters per Second per Second
    static constexpr double R  = 287.053;      //Atmospheric Gas Constant Joules per Kilogram per Kelvin

    //static functions
    static double getAmbientPressure(double h);     //Function Returns Ambient Pressure from Standard Atmosphere Model
    static double getAmbientTemperature(double h);  //Function Returns Ambient Temperature from Standard Atmosphere Model
    static double getAmbientDensity(double h);      //Function Returns Ambient Density from Standard Atmosphere Model
    static double getVelocity();                    //Function Returns Velocity

    //initialization and running functions
    virtual void init(Job* job);
    virtual void run(Job* job);
};

class FVMBolideRestartDriver : public FVMBolideImpactDriver{
public:
    FVMBolideRestartDriver(){
        FVMVariableStepDriver::object_name = "FVMBolideRestartDriver";
        FVMBolideImpactDriver::object_name = "FVMBolideRestartDriver";
        object_name = "FVMBolideRestartDriver"; //set object name here
    }

    //filenames to read data from
    std::string point_file, volume_file;

    //initialization and running functions
    virtual void init(Job* job);
    virtual void restart(Job* job);
    virtual void run(Job* job);
};

#endif //MPM_V3_FVM_DRIVERS_HPP

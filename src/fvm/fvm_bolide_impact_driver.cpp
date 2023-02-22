//
// Created by hemi-user on 2/21/23.
//

//
// Created by aaron on 12/12/19.
// finite_volume_driver.cpp
//

//
// Created by aaron on 5/10/18.
// default_driver.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <time.h>

#include "job.hpp"
#include "mpm_objects.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"

#include "fvm_objects.hpp"
#include "fvm_drivers.hpp"
#include "fvm_artificial_viscosity.hpp"
#include "registry.hpp"

// Initialize Static Members
double FVMBolideImpactDriver::H;        //Height in Meters
double FVMBolideImpactDriver::V;        //Velocity in Meters/Second
double FVMBolideImpactDriver::Y0;       //Initial Leading Position of Projectile in Meters

//constant variables
const int FVMBolideImpactDriver::bmax;
const double FVMBolideImpactDriver::Hb[8];       //Heights of Model in Meters
const double FVMBolideImpactDriver::Lb[8];       //Gradients of Model in Kelvin per Meter
const double FVMBolideImpactDriver::Tb[8];       //Thermal Layers of Model in Kelvin
const double FVMBolideImpactDriver::Pb[8];       //Pressure Payers of Model in Pascal
const double FVMBolideImpactDriver::g0;       //Gravitational Acceleration in Meters per Second per Second
const double FVMBolideImpactDriver::R;        //Atmospheric Gas Constant Joules per Kilogram per Kelvin

/*----------------------------------------------------------------------------*/

void FVMBolideImpactDriver::init(Job* job){

    //call finite volume driver initializer
    FiniteVolumeDriver::init(job);

    if (TYPE != THERMAL){
        std::cerr << "FVMBolideImpactDriver requires THERMAL simulation. Exiting." << std::endl;
        exit(0);
    } else if (fluid_grid->object_name.compare("FVMGmsh3D") != 0){
        std::cerr << fluid_grid->object_name << " : FVMBolideImpactDriver requires FVMGmsh3D grid. Exiting." << std::endl;
        exit(0);
    }

    int GRID_DIM = 0;
    //assign grid dimension from job type
    if (job->JOB_TYPE == job->JOB_1D){
        GRID_DIM = 1;
    } else if (job->JOB_TYPE == job->JOB_2D){
        GRID_DIM = 2;
    } else if (job->JOB_TYPE == job->JOB_3D){
        GRID_DIM = 3;
    } else if (job->JOB_TYPE == job->JOB_2D_OOP){
        GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
    } else if (job->JOB_TYPE == job->JOB_AXISYM){
        GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
    } else {
        std::cerr << "Job doesn't have defined type for input " << job->JOB_TYPE << "." << std::endl;
        exit(0);
    }

    //check size of properties passed to driver object
    if (fp64_props.size() < 1+GRID_DIM+2) {
        //use default properties
        dt0 = job->dt;
        H = 0;
        V = 11e3;
        lambda = 0.5;
        eta = 0;
    } else if (fp64_props.size() < 1+GRID_DIM+4) {
        //assign initial height and velocity
        dt0 = job->dt;
        H = fp64_props[1 + GRID_DIM];
        V = fp64_props[2 + GRID_DIM];
        lambda = 0.5;
        eta = 0;
    } else {
        //assign initial height, velocity, and time-increment parameters
        dt0 = job->dt;
        H = fp64_props[1 + GRID_DIM];
        V = fp64_props[2 + GRID_DIM];
        lambda = fp64_props[3 + GRID_DIM];
        eta = fp64_props[4 + GRID_DIM];
    }

    //loop over str-props and assign relevant flags
    std::vector<std::string> options = {"USE_ARTIFICIAL_VISCOSITY"};
    for (int i=0; i<str_props.size(); i++){
        switch (Parser::findStringID(options, str_props[i])){
            case 0:
                //USE_ARTIFICIAL_VISCOSITY
                USE_ARTIFICIAL_VISCOSITY = true;
                std::cout << "Using artificial viscosity." << std::endl;
                break;
            default:
                //do nothing
                break;
        }
    }

    // set up thermal floor
    USE_THERMAL_FLOOR = true;
    cv_floor = 717.63; T_floor = 100;

    std::cout << "Using variable time-stepping. Default: " << dt0 << "s, Mesh Quality: " << lambda << ", Viscosity: " << eta << std::endl;
    std::cout << "Bolide Impact Simulation:\n";
    std::cout << "    Inclination: 45 Degrees\n";
    std::cout << "    Initial Velocity: " << V << " m/s\n";
    std::cout << "    Initial Height: " << H << " m\n";
    std::cout << "    Initial Temperature: " << getAmbientTemperature(H) << " K\n";
    std::cout << "    Initial Pressure: " << getAmbientPressure(H) << " Pa\n";
    std::cout << "    Initial Density: " << getAmbientDensity(H) << " kg/m^3" << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

double FVMBolideImpactDriver::getAmbientPressure(double h) {
    // Output
    double T = Tb[0];
    double P = Pb[0];
    bool success = false;

    // Standard Atmosphere Model (1976)
    for (int i=0; i<bmax; i++){
        if (!success && h < Hb[i+1]){
            // Determine Temperature
            T = Tb[i] + Lb[i] * (h - Hb[i]);

            // Determine Correct Model For Pressure
            if (std::abs(Lb[i]) > 1e-10){
                P = Pb[i] * std::pow((Tb[i] / T), (g0 / (R * Lb[i])));
            } else {
                P = Pb[i] * std::exp(-g0 * (h - Hb[i]) / (R * Tb[i]));
            }
            success = true;
        }
    }

    if (!success){
        // Height Above 85km
        P = Pb[bmax];
    }

    // Return Temperature
    return P;
}

double FVMBolideImpactDriver::getAmbientTemperature(double h){
    // Output
    double T = Tb[0];
    bool success = false;

    // Standard Atmosphere Model (1976)
    for (int i=0; i<bmax; i++){
        if (!success && h < Hb[i+1]){
            T = Tb[i] + Lb[i] * (h - Hb[i]);
            success = true;
        }
    }

    if (!success){
        // Height Above 85km
        T = Tb[bmax];
    }

    // Return Temperature
    return T;
}

double FVMBolideImpactDriver::getAmbientDensity(double h){

    double T = getAmbientTemperature(h);
    double P = getAmbientPressure(h);

    // Compute Density
    return (P / (R*T));
}


double FVMBolideImpactDriver::getVelocity(){
    return V;
}

/*----------------------------------------------------------------------------*/

void FVMBolideImpactDriver::run(Job* job) {
    //initialize FVM objects
    serializer->init(job, this);
    fluid_grid->init(job, this);
    solver->init(job, this);
    fluid_material->init(job, this);
    fluid_body->init(job, this);

    //set counters to zero
    int stepCount = 0;
    int frameCount = 0;

    struct timespec timeStart, timeFrame, timeFinish, timeStep;
    clock_gettime(CLOCK_MONOTONIC, &timeStart);
    timeFrame = timeStart;
    timeStep = timeStart;
    //clock_t clockSim = clock();
    //clock_t clockFrame = clock();
    double tSim = 0;
    double tFrame = 0;
    double tStep = 0;

    //initialize gravity
    generateGravity(job);
    applyGravity(job);

    //bolide simulation variables
    double V0 = V;  //initial velocity
    double H0 = H;  //initial height
    double cosTHETA = std::sqrt(0.5); //cosine of inclination angle
    double YMin = Y0;  //minimum bolide point y-position
    double YMax = Y0;  //maximum bolide point y-position (for radius calculation)
    double YDot = 0;   //y-velocity of leading bolide point
    double YTmp = 0;
    int YMinID  = 0;   //ID of point with minimum y-position

    bool YMin_Set = false;
    for (int i=0; i<job->bodies[0]->points->x.size(); i++){
        // Only Consider Active Points
        if (job->bodies[0]->points->active(i)){
            if (!YMin_Set){
                //Initialize YMin
                YMin = job->bodies[0]->points->x(i,1);
                YMax = job->bodies[0]->points->x(i,1);
                YMin_Set = true;
            } else {
                YTmp = job->bodies[0]->points->x(i,1);
                if (YTmp < YMin){
                    YMin = YTmp;
                } else if (YTmp > YMax){
                    YMax = YTmp;
                }
            }
        }
    }
    Y0 = YMin;

    // Compute Radius and Centroid of Bolide
    double R0 = (YMax - YMin) / 2.0;
    double YMed = (YMax + YMin) / 2.0;
    KinematicVector X0 = KinematicVector(job->JOB_TYPE);
    KinematicVector XTmp = KinematicVector(job->JOB_TYPE);
    X0.setZero(); X0(1) = YMed;

    // Set Initial Conditions
    if (true){
        for (int e=0; e<fluid_grid->element_count; e++){
            XTmp = fluid_grid->getElementCentroid(job, e);
            // Check That Element is Outside of Bolide
            if ((XTmp - X0).dot(XTmp - X0) > 2.25*R0*R0) {
                // Subtract Kinetic Energy From Total Energy
                fluid_body->rhoE(e) -= 0.5 * fluid_body->p(e).dot(fluid_body->p(e)) / fluid_body->rho(e);
                // Update Momentum
                fluid_body->p(e, 1) = fluid_body->rho(e) * V;
                // Add Kinetic Energy To Total Energy
                fluid_body->rhoE(e) += 0.5 * fluid_body->p(e).dot(fluid_body->p(e)) / fluid_body->rho(e);
            }
        }
    }


    //create new output file
    //open and clear file
    std::ostringstream oss;
    oss << "Bolide_History_" << (V/1e3) << "kms_" << (H/1e3) << "km.csv";
    std::string output_filename = oss.str();
    std::ofstream file (output_filename,std::ios::trunc);
    if (file.is_open()){
        //success!
        //write file header
        file << "Time [s], Time Increment [s], Velocity [km/s], Altitude [km], Pressure [Pa], Temperature [K], Density [kg/m^3], Leading Edge Position [m]\n";
        file << job->t << ", ";
        file << job->dt << ", ";
        file << (V/1e3) << ", ";
        file << (H/1e3) << ", ";
        file << getAmbientPressure(H) << ", ";
        file << getAmbientTemperature(H) << ", ";
        file << getAmbientDensity(H) << ", ";
        file << YMin << "\n";
        file.close();
    } else {
        std::cerr << "ERROR! Cannot open " << output_filename << "! Exiting." << std::endl;
        exit(0);
    }

    //element length scales
    Eigen::VectorXd l(fluid_grid->element_count);
    l.setZero();
    std::vector<int> e_faces;
    double l_min, l_e;
    for (int e=0; e<fluid_grid->element_count; e++){
        e_faces = fluid_grid->getElementFaces(e);
        for (int f=0; f<e_faces.size(); f++){
            l_e = (fluid_grid->getElementCentroid(job, e) - fluid_grid->getFaceCentroid(job, e_faces[f])).norm();
            if (f==0 || l_e < l_min){
                l_min = l_e;
            }
        }
        l(e) = l_min;
    }

    //temporary velocity magnitude and stable time step
    double c, v, dts;

    //artificial viscosity
    std::vector<ArtificialViscosityCalculator::fluxVector> fluxVectors;
    if (USE_ARTIFICIAL_VISCOSITY){
        fluxVectors.resize(fluid_grid->face_count);
    }

    //reality check
    bool failed = false;

    //run simulation until stop_time
    while (job->t <= stop_time){

        //FIRST: Identify Leading Material Point
        bool YMin_Set = false;
        YMinID = 0;
        for (int i=0; i<job->bodies[0]->points->x.size(); i++){
            // Only Consider Active Points
            if (job->bodies[0]->points->active(i)){
                if (!YMin_Set){
                    //Initialize YMin
                    YMin = job->bodies[0]->points->x(i,1);
                    YMinID = i;
                    YMin_Set = true;
                } else {
                    YTmp = job->bodies[0]->points->x(i,1);
                    if (YTmp < YMin){
                        YMin = YTmp;
                        YMinID = i;
                    }
                }
            }
        }

        //SECOND: Determine Velocity of Leading Material Point
        if (YMin_Set) {
            YDot = job->bodies[0]->points->x_t(YMinID, 0);
        }

        //THIRD: Adjust Reference Frame of Solution Fields to New Velocity
        if (YMin_Set){
            // (1) Adjust Material Point Velocities and Momenta
            for (int i=0; i<job->bodies[0]->points->x.size(); i++){
                if (job->bodies[0]->points->active(i)){
                    job->bodies[0]->points->x_t(i,1) -= YDot;
                    job->bodies[0]->points->mx_t(i) = job->bodies[0]->points->m(i) * job->bodies[0]->points->x_t(i);
                }
            }

            // (2) Adjust Fluid Momenta and Energies
            for (int e=0; e<fluid_grid->element_count; e++){
                // Subtract Kinetic Energy From Total Energy
                fluid_body->rhoE(e) -= 0.5 * fluid_body->p(e).dot(fluid_body->p(e)) / fluid_body->rho(e);
                // Update Momentum
                fluid_body->p(e,1) -= fluid_body->rho(e) * YDot;
                // Add Kinetic Energy To Total Energy
                fluid_body->rhoE(e) += 0.5 * fluid_body->p(e).dot(fluid_body->p(e)) / fluid_body->rho(e);
            }

            // (3) Adjust Velocity BC
            V -= YDot;
        }


        //FOURTH: check time increment for stability
        job->dt = dt0;
        for (int e=0; e<fluid_grid->element_count; e++){
            //check advection
            v = fluid_body->p(e).norm()/fluid_body->rho(e);
            c = fluid_material->getSpeedOfSound(job,
                                                this,
                                                fluid_body->rho(e),
                                                fluid_body->p(e),
                                                fluid_body->rhoE(e),
                                                fluid_body->n_e(e));
            dts = lambda * l(e)/(c+v);
            if (dts < job->dt){
                job->dt = dts;
            }

            //check viscosity
            if (eta > 0){
                dts = (lambda * lambda * l(e) * l(e) * fluid_body->rho(e)) / eta;
                if (dts < job->dt){
                    job->dt = dts;
                }
            }

            //check thermal floor
            if (TYPE == THERMAL && USE_THERMAL_FLOOR && fluid_material->getTemperature(job,
                                                                                       this,
                                                                                       fluid_body->rho(e),
                                                                                       fluid_body->p(e),
                                                                                       fluid_body->rhoE(e),
                                                                                       fluid_body->n_e(e)) <= T_floor){

                //Highlight Error
                std::cout << std::endl;
                std::cout << "    Thermal Floor: [" << e << "]: " << fluid_material->getTemperature(job,
                                                                                                    this,
                                                                                                    fluid_body->rho(e),
                                                                                                    fluid_body->p(e),
                                                                                                    fluid_body->rhoE(e),
                                                                                                    fluid_body->n_e(e));
                std::cout << " --> " << T_floor << std::endl;

                //Assign New Temperature
                fluid_body->rhoE(e) = fluid_body->rho(e) * cv_floor * T_floor
                                      + 0.5*fluid_body->p[e].dot(fluid_body->p[e])/fluid_body->rho(e);

            }

            if (!failed && (fluid_material->getTemperature(job,
                                                           this,
                                                           fluid_body->rho(e),
                                                           fluid_body->p(e),
                                                           fluid_body->rhoE(e),
                                                           fluid_body->n_e(e)) <= 0
                            ||
                            std::isnan(fluid_material->getTemperature(job,
                                                                      this,
                                                                      fluid_body->rho(e),
                                                                      fluid_body->p(e),
                                                                      fluid_body->rhoE(e),
                                                                      fluid_body->n_e(e))))){
                failed = true;
                std::cout << std::endl << std::endl;
                std::cout << "    First failue: [" << e << "]" << std::endl;
                std::cout << "    Density: " << fluid_body->rho(e) << std::endl;
                std::cout << "    Momentum: " << fluid_body->p(e).norm() << std::endl;
                std::cout << "    Energy: " << fluid_body->rhoE(e)<< std::endl;
                std::cout << "    Porosity: " << fluid_body->n_e(e)<< std::endl;
                std::cout << "    Position: " << EIGEN_MAP_OF_KINEMATIC_VECTOR(fluid_grid->getElementCentroid(job,e)).transpose();
                std::cout << std::endl;
                std::cout << "Exiting." << std::endl;
                exit(0);
            }
        }
        //report step length to console
        if (job->dt < dt0){
            std::cout << "dt: " << job->dt << "    \r" << std::flush;
        }


        // FIFTH: run solver
        solver->step(job,this);

        // SIXTH: Write Data to Files
        if (job->serializer->writeFrame(job) == 1) {
            //call fvm serializer to write frame as well:
            serializer->writeFrame(job,this);

            //write bolide simulation data
            file.open(output_filename,std::ios::app);
            if (file.is_open()){
                //success!
                file << job->t << ", ";
                file << job->dt << ", ";
                file << (V/1e3) << ", ";
                file << (H/1e3) << ", ";
                file << getAmbientPressure(H) << ", ";
                file << getAmbientTemperature(H) << ", ";
                file << getAmbientDensity(H) << ", ";
                file << YMin << "\n";
                file.close();
            } else {
                std::cerr << "ERROR! Cannot open " << output_filename << "! Exiting." << std::endl;
                exit(0);
            }

            //successful frame written
            //tFrame = (double)(clock() - clockFrame)/CLOCKS_PER_SEC;
            //tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
            clock_gettime(CLOCK_MONOTONIC,&timeFinish);
            tFrame = (timeFinish.tv_sec - timeFrame.tv_sec) + (timeFinish.tv_nsec - timeFrame.tv_nsec)/1000000000.0;
            tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
            timeFrame = timeFinish;
            printf("\33[2K");
            std::cout << "Frame Written [" << ++frameCount << "]. Time/Frame [" << tFrame << " s]. Elapsed Time [" << tSim << " s]." << std::flush;
        }
        std::cout << "\r";

        //SEVENTH: apply artificial viscosity
        if (USE_ARTIFICIAL_VISCOSITY){
            applyArtificialViscosityFluxes(job, fluxVectors);
        }

        //EIGHTH: update altitude
        H -= cosTHETA * V * job->dt;

        //NINTH: update time
        job->t += job->dt;
    }
    //tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
    clock_gettime(CLOCK_MONOTONIC,&timeFinish);
    tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
    std::cout << std::endl << std::endl << "Simulation Complete. Elapsed Time [" << tSim << "s]." << std::endl;
    return;
}

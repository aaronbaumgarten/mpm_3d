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

/*----------------------------------------------------------------------------*/

void FVMVariableStepDriver::init(Job* job){

    //call parent initializer
    FiniteVolumeDriver::init(job);

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
    if (fp64_props.size() <= 1+GRID_DIM) {
        //use default properties
        dt0 = job->dt;
        lambda = 0.5;
        eta = 0;
        T_floor = 0;
        cv_floor = 0;
    } else if (fp64_props.size() == 2+GRID_DIM) {
        dt0 = job->dt;
        lambda = fp64_props[1+GRID_DIM];
        eta = 0;
        T_floor = 0;
        cv_floor = 0;
    } else if (fp64_props.size() <= 4+GRID_DIM) {
        dt0 = job->dt;
        lambda = fp64_props[1+GRID_DIM];
        eta = fp64_props[2+GRID_DIM];
        T_floor = 0;
        cv_floor = 0;
    } else if (fp64_props.size() > 4+GRID_DIM) {
        dt0 = job->dt;
        lambda = fp64_props[1+GRID_DIM];
        eta = fp64_props[2+GRID_DIM];
        T_floor = fp64_props[3+GRID_DIM];
        cv_floor = fp64_props[4+GRID_DIM];
    }

    //loop over str-props and assign relevant flags
    std::vector<std::string> options = {"USE_ARTIFICIAL_VISCOSITY", "USE_THERMAL_FLOOR"};
    for (int i=0; i<str_props.size(); i++){
        switch (Parser::findStringID(options, str_props[i])){
            case 0:
                //USE_ARTIFICIAL_VISCOSITY
                USE_ARTIFICIAL_VISCOSITY = true;
                std::cout << "Using artificial viscosity." << std::endl;
                break;
            case 1:
                //USE_THERMAL_FLOOR
                USE_THERMAL_FLOOR = true;
                std::cout << "Using thermal floor. WARNING: energy WILL be inserted into simulation. T = " << T_floor;
                std::cout << " [K]. cv = " << cv_floor << " J/kg*K." << std::endl;
                break;
            default:
                //do nothing
                break;
        }
    }

    std::cout << "Using variable time-stepping. Default: " << dt0 << "s, Mesh Quality: " << lambda << ", Viscosity: " << eta << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void FVMVariableStepDriver::run(Job* job) {
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
        //check time increment
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

        //run solver
        solver->step(job,this);

        //clock_gettime(CLOCK_MONOTONIC, &timeFinish);
        //tStep = (timeFinish.tv_sec - timeStep.tv_sec) + (timeFinish.tv_nsec - timeStep.tv_nsec)/1000000000.0;
        //timeStep = timeFinish;

        //std::cout << "Step Completed [" << ++stepCount << "]. Step Time [" << tStep << " s]." << std::flush;
        if (job->serializer->writeFrame(job) == 1) {
            //call fvm serializer to write frame as well:
            serializer->writeFrame(job,this);

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

        //apply artificial viscosity
        if (USE_ARTIFICIAL_VISCOSITY){
            applyArtificialViscosityFluxes(job, fluxVectors);
        }

        job->t += job->dt;
    }
    //tSim = (double)(clock() - clockSim)/CLOCKS_PER_SEC;
    clock_gettime(CLOCK_MONOTONIC,&timeFinish);
    tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
    std::cout << std::endl << std::endl << "Simulation Complete. Elapsed Time [" << tSim << "s]." << std::endl;
    return;
}

//parallel functions
void FVMVariableStepDriver::applyArtificialViscosityFluxes(Job* job,
                                                           std::vector<ArtificialViscosityCalculator::fluxVector> &fluxVectors){
    if (job->thread_count > 1){
        //determine number of threads for element gradient calculation
        int thread_count = 0;
        if (fluid_grid->face_count >= job->thread_count){
            thread_count = job->thread_count;
        } else {
            thread_count = fluid_grid->face_count;
        }

        //boolean of completion status
        //volatile bool firstTaskComplete[thread_count] = {false};
        volatile bool *firstTaskComplete = new volatile bool[thread_count];
        for (int t=0; t<thread_count; t++){
            firstTaskComplete[t] = false;
        }

        //choose interval size
        int k_max = fluid_grid->face_count - 1;
        int k_interval = (fluid_grid->face_count/thread_count) + 1;
        int k_begin, k_end;

        for (int t=0; t<thread_count; t++) {
            //set interval
            k_begin = t * k_interval;
            k_end = k_begin + k_interval - 1;
            if (k_end > k_max){
                k_end = k_max;
            }

            //send job to threadpool
            job->threadPool.doJob(std::bind(calculateSubsetOfAVFluxes,
                                           job,
                                           this,
                                           std::ref(fluxVectors),
                                           k_begin,
                                           k_end,
                                           std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }

        delete[] firstTaskComplete;

    } else {
        //call construction function
        volatile bool taskDone = false;
        calculateSubsetOfAVFluxes(job, this, fluxVectors, 0, fluid_grid->face_count-1, taskDone);
    }

    //add vectors
    std::array<int,2> f_elem = std::array<int,2>();
    for (int f=0; f<fluid_grid->face_count; f++){
        f_elem = fluid_grid->getOrientedElementsByFace(f);
        if (f_elem[0] > -1 && f_elem[1] > -1) {
            fluid_body->rho(f_elem[0]) -= job->dt * fluxVectors[f].rho / fluid_grid->getElementVolume(f_elem[0]);
            fluid_body->rho(f_elem[1]) += job->dt * fluxVectors[f].rho / fluid_grid->getElementVolume(f_elem[1]);
            fluid_body->p(f_elem[0]) -= job->dt * fluxVectors[f].p / fluid_grid->getElementVolume(f_elem[0]);
            fluid_body->p(f_elem[1]) += job->dt * fluxVectors[f].p / fluid_grid->getElementVolume(f_elem[1]);
            fluid_body->rhoE(f_elem[0]) -= job->dt * fluxVectors[f].rhoE / fluid_grid->getElementVolume(f_elem[0]);
            fluid_body->rhoE(f_elem[1]) += job->dt * fluxVectors[f].rhoE / fluid_grid->getElementVolume(f_elem[1]);
        }
    }

    return;
}

void FVMVariableStepDriver::calculateSubsetOfAVFluxes(Job* job,
                                                       FiniteVolumeDriver* driver,
                                                       std::vector<ArtificialViscosityCalculator::fluxVector> &fluxVectors,
                                                       int f_begin, int f_end,
                                                       volatile bool &done){
    std::array<int,2> f_elem = std::array<int,2>();
    for (int f=f_begin; f<=f_end; f++){
        f_elem = driver->fluid_grid->getOrientedElementsByFace(f);
        if (f_elem[0] > -1 && f_elem[1] > -1) {
            fluxVectors[f] = ArtificialViscosityCalculator::getArtificialViscosityFlux(job, driver, f);
        }
    }

    done = true;

    return;
}

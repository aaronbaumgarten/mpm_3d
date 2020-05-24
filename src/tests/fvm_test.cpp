//
// Created by aaron on 12/23/19.
// fvm_test.cpp.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <fstream>
#include <job.hpp>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "mpm_objects.hpp"
#include "fvm_objects.hpp"

#include "registry.hpp"
#include "job.hpp"

#include "fvm/fvm_grids.hpp"
#include "fvm/fvm_bodies.hpp"
#include "fvm/fvm_materials.hpp"
#include "fvm/fvm_solvers.hpp"
#include "fvm/fvm_serializers.hpp"
#include "objects/grids/grids.hpp"

/*
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
};
 */

namespace FVM_TEST{
    class FVMTestDriver : public FiniteVolumeDriver {
    public:
        FVMTestDriver(){
            object_name = "FVMTestDriver"; //set object name here
        }

        int GRID_DIM;

        //functions which must be implemented by every driver
        void init(Job* job){
            //force job properties
            job->t = 0;
            job->JOB_TYPE = job->JOB_3D;

            job->grid = std::unique_ptr<Grid>(new CartesianLinear);
            job->grid->node_count = 1;

            GRID_DIM = 3;

            //set simulation ORDER
            ORDER = 1;

            //set simulation TYPE
            TYPE = FiniteVolumeDriver::THERMAL;

            //set gravity
            gravity = KinematicVector(job->JOB_TYPE);
            gravity.setZero();

            //assign grid
            fluid_grid = std::unique_ptr<FiniteVolumeGrid>(new FVMCartesian);
            fluid_grid->fp64_props = {2, 1, 1, 0, 0, 0, 0, 0, 0}; //Lx
            fluid_grid->int_props = {2, 1, 1, 8, 8, 8, 8, 8, 8};  //Nx

            //assign body
            fluid_body = std::unique_ptr<FiniteVolumeBody>(new FVMDefaultBody);
            fluid_body->fp64_props = {1.177, 298}; //rho and theta

            //assign material
            fluid_material = std::unique_ptr<FiniteVolumeMaterial>(new FVMIdealGas);
            fluid_material->fp64_props = {1.4, 287.1, 0, 0}; //gamma, R, eta, k

            //assign solver
            solver = std::unique_ptr<FiniteVolumeSolver>(new FVMDefaultSolver);

            //set time step for stability
            job->dt = 0.001; //< dx*sqrt(rho/kappa)
            stop_time = 1.0; //seconds

            return;
        }

        void run(Job* job){
            //initialize FVM objects
            fluid_grid->init(job, this);
            fluid_body->init(job, this);
            fluid_material->init(job, this);
            solver->init(job, this);

            //set up specific problems
            KinematicVector p_0 = KinematicVector(job->JOB_TYPE);
            KinematicVector p_1 = KinematicVector(job->JOB_TYPE);
            double rho_0, rho_1, rhoE_0, rhoE_1;

            //first check number of elements
            std::cout << "Number of Elements: " << fluid_grid->element_count << " ?= 2" << std::endl;

            //check element volumes
            std::cout << "Element Volumes: " << fluid_grid->getElementVolume(0) << " ?= 1" << std::endl;

            /*---------------------------------*/
            //Flux Test
            /*---------------------------------*/

            p_0[0] = std::rand()/(RAND_MAX+1.0);
            p_1[0] = std::rand()/(RAND_MAX+1.0);

            p_0[1] = std::rand()/(RAND_MAX+1.0);
            p_1[1] = std::rand()/(RAND_MAX+1.0);

            p_0[2] = std::rand()/(RAND_MAX+1.0);
            p_1[2] = std::rand()/(RAND_MAX+1.0);

            rho_0 = 0.5 + std::rand()/(RAND_MAX+1.0);
            rho_1 = 0.5 + std::rand()/(RAND_MAX+1.0);

            rhoE_0 = 0.5 + std::rand()/(RAND_MAX+1.0);
            rhoE_1 = 0.5 + std::rand()/(RAND_MAX+1.0);

            rhoE_0 *= 100000;
            rhoE_1 *= 100000;

            fluid_body->rho(0) = rho_0;
            fluid_body->rho(1) = rho_1;

            fluid_body->p[0] = p_0;
            fluid_body->p[1] = p_1;

            fluid_body->rhoE(0) = rhoE_0;
            fluid_body->rhoE(1) = rhoE_1;

            //call grid to reconstruct conserved fields
            fluid_grid->constructDensityField(job, this);
            fluid_grid->constructMomentumField(job, this);
            fluid_grid->constructEnergyField(job, this);

            //fluid fluxes associated with each volume (to compare against exact solution)
            Eigen::VectorXd density_fluxes = fluid_grid->calculateElementMassFluxes(job, this);
            KinematicVectorArray momentum_fluxes = fluid_grid->calculateElementMomentumFluxes(job, this);
            Eigen::VectorXd energy_fluxes = fluid_grid->calculateElementEnergyFluxes(job, this);

            //flux into cell 0 should be zero
            std::cout << "2 Cell Check #1:" << std::endl;
            std::cout << "rho:  " << rho_0 << " , " << rho_1 << std::endl;
            std::cout << "p:    " << p_0[0] << " , " << p_1[0] << std::endl;
            std::cout << "rhoE: " << rhoE_0 << " , " << rhoE_1 << std::endl;

            //analytical solution
            double P_0 = (rhoE_0 - 0.5*p_0.dot(p_0)/rho_0)*0.4;
            double P_1 = (rhoE_1 - 0.5*p_1.dot(p_1)/rho_1)*0.4;

            KinematicVector u_bar = (p_0/std::sqrt(rho_0) + p_1/std::sqrt(rho_1))
                                    /(std::sqrt(rho_0) + std::sqrt(rho_1));
            double H_bar = ((rhoE_0 + P_0)/std::sqrt(rho_0) + (rhoE_1 + P_1)/std::sqrt(rho_1))
                           /(std::sqrt(rho_0) + std::sqrt(rho_1));

            double a = std::sqrt(0.4*(H_bar - 0.5*u_bar.dot(u_bar)));

            double lambda_1 = std::abs(u_bar[0] - a);
            double lambda_2 = std::abs(u_bar[0]);
            double lambda_3 = std::abs(u_bar[0]);
            double lambda_4 = std::abs(u_bar[0] + a);

            std::cout << "lambdas: " << lambda_1 << ", " << lambda_4 << " ?= " << a << std::endl;

            KinematicVector s_bar = u_bar;
            s_bar[0] = 0;
            KinematicVector s_0 = p_0/rho_0;
            KinematicVector s_1 = p_1/rho_1;
            s_0[0] = 0;
            s_1[0] = 0;

            KinematicVector sa_2 = ((p_1 - p_0) - u_bar*(rho_1 - rho_0));
            sa_2[0] = 0;

            double a_3 = ((H_bar - u_bar.dot(u_bar))*(rho_1 - rho_0)
                          + u_bar[0]*(p_1[0] - p_0[0])
                          + s_bar.dot(rho_1 * s_1 - rho_0*s_0)
                          - (rhoE_1 - rhoE_0)) / (H_bar - 0.5*u_bar.dot(u_bar));

            double a_4 = 0.5/a * (p_1[0] - p_0[0] - (rho_1 - rho_0)*(u_bar[0] - a) - a*a_3);
            double a_1 = (rho_1 - rho_0) - a_3 - a_4;

            //solution
            double density_flux = 0.5*(p_0[0] + p_1[0]) - 0.5*(a_1*lambda_1 + a_3*lambda_3 + a_4*lambda_4);

            double momentum_flux = 0.5*(p_0[0]*p_0[0]/rho_0 + p_1[0]*p_1[0]/rho_1)
                                     + 0.5*(P_0 + P_1) - (rhoE_1 - 0.5*rho_1*s_1.dot(s_1))*0.4
                                     - 0.5*(a_1*lambda_1*(u_bar[0] - a)
                                            + a_3*lambda_3*u_bar[0]
                                            + a_4*lambda_4*(u_bar[0] + a));

            double energy_flux = 0.5*((rhoE_0 + P_0)*p_0[0]/rho_0 + (rhoE_1 + P_1)*p_1[0]/rho_1)
                                 - 0.5*(a_1*lambda_1*(H_bar - u_bar[0]*a)
                                        + lambda_2*sa_2.dot(s_bar)
                                        + a_3*lambda_3*0.5*u_bar.dot(u_bar)
                                        + a_4*lambda_4*(H_bar + u_bar[0]*a));

            std::cout << "density flux: " << density_fluxes(1) << " ?= " << density_flux;
            std::cout << " : " << (density_fluxes(1) - density_flux)/density_flux << std::endl;

            std::cout << "momentum flux: " << momentum_fluxes(1,0) << " ?= " << momentum_flux;
            std::cout << " : " << (momentum_fluxes(1,0) - momentum_flux)/momentum_flux << std::endl;

            std::cout << "energy flux: " << energy_fluxes(1) << " ?= " << energy_flux;
            std::cout << " : " << (energy_fluxes(1) - energy_flux)/energy_flux << std::endl;

            return;
        }
    };
}

void fvm_test(Job* job){
    std::unique_ptr<FiniteVolumeDriver> driver = std::unique_ptr<FiniteVolumeDriver>(new FVM_TEST::FVMTestDriver);
    driver->init(job);
    driver->run(job);
    return;
}
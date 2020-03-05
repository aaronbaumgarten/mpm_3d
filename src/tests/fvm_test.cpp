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
            ORDER = 2;

            //set simulation TYPE
            TYPE = FiniteVolumeDriver::ISOTHERMAL;

            //set gravity
            gravity = KinematicVector(job->JOB_TYPE);
            gravity.setZero();

            //assign grid
            fluid_grid = std::unique_ptr<FiniteVolumeGrid>(new FVMCartesian);
            fluid_grid->fp64_props = {3, 1, 1}; //Lx
            fluid_grid->int_props = {3, 1, 1};  //Nx

            //assign body
            fluid_body = std::unique_ptr<FiniteVolumeBody>(new FVMDefaultBody);
            fluid_body->fp64_props = {1, 1}; //rho and theta

            //assign material
            fluid_material = std::unique_ptr<FiniteVolumeMaterial>(new FVMBarotropicViscousFluid);
            fluid_material->fp64_props = {100, 0, 1}; //kappa, eta, rho_0

            //assign solver
            solver = std::unique_ptr<FiniteVolumeSolver>(new FVMDefaultSolver);

            //set time step for stability
            job->dt = 0.01; //< dx*sqrt(rho/kappa)
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
            KinematicVector p_2 = KinematicVector(job->JOB_TYPE);
            double rho_0, rho_1, rho_2;

            //first check number of elements
            std::cout << "Number of Elements: " << fluid_grid->element_count << " ?= 3" << std::endl;

            //check element volumes
            std::cout << "Element Volumes: " << fluid_grid->getElementVolume(0) << " ?= 1" << std::endl;

            /*---------------------------------*/
            //Problem 1
            /*---------------------------------*/

            //let p_0 = p_1 = RAND =/= p_2
            p_0[0] = std::rand()/(RAND_MAX+1.0);
            p_1[0] = p_0[0];
            p_2[0] = -std::rand()/(RAND_MAX+1.0);

            //let rho_0 = rho_1 = RAND =/= rho_2
            rho_0 = 0.5 + std::rand()/(RAND_MAX+1.0);
            rho_1 = rho_0;
            rho_2 = 0.5 + std::rand()/(RAND_MAX+1.0);

            fluid_body->rho(0) = rho_0;
            fluid_body->rho(1) = rho_1;
            fluid_body->rho(2) = rho_2;

            fluid_body->p[0] = p_0;
            fluid_body->p[1] = p_1;
            fluid_body->p[2] = p_2;

            //call grid to reconstruct conserved fields
            fluid_grid->constructDensityField(job, this);
            fluid_grid->constructMomentumField(job, this);

            //fluid fluxes associated with each volume (to compare against exact solution)
            Eigen::VectorXd density_fluxes = fluid_grid->calculateElementMassFluxes(job, this);
            KinematicVectorArray momentum_fluxes = fluid_grid->calculateElementMomentumFluxes(job, this);

            //flux into cell 0 should be zero
            std::cout << "3 Cell Check #1:" << std::endl;
            std::cout << "rho: " << rho_0 << " , " << rho_1 << " , " << rho_2 << std::endl;
            std::cout << "p:   " << p_0[0] << " , " << p_1[0] << " , " << p_2[0] << std::endl;

            std::cout << "mass_flux[0]     = " << density_fluxes(0) << " ?= " << -p_0[0] << std::endl;
            std::cout << "momentum_flux[0] = <" << momentum_fluxes(0,0) << ", " << momentum_fluxes(0,1) << "> ?= <";
            std::cout << -p_0[0]*p_0[0]/rho_0 << ", 0>" << std::endl;

            //exact solution for second boundary
            double rho_bar = std::sqrt(rho_1*rho_2);
            double u_bar = (p_1[0]/std::sqrt(rho_1) + p_2[0]/std::sqrt(rho_2))/(std::sqrt(rho_1) + std::sqrt(rho_2));
            double a_1 = 0.5*(rho_2 - rho_1) - 0.5*std::sqrt(rho_bar)/10.0 * (p_2[0] - p_1[0] - (rho_2 - rho_1)*u_bar);
            double a_2 = rho_2 - rho_1 - a_1;
            double dm_dt = 0.5*(p_1[0] + p_2[0] - a_1*std::abs(u_bar - 10.0/std::sqrt(rho_bar))
                                                - a_2*std::abs(u_bar + 10.0/std::sqrt(rho_bar))) - p_0[0];
            /*
            double dp_dt = 0.5*(p_1[0]*p_1[0]/rho_1 + p_2[0]*p_2[0]/rho_2)
                           + 100.0*std::log(rho_bar/1.0) + 0.5*(100.0/rho_bar*(rho_1 + rho_2 - 2*rho_bar))
                            - 0.5*(a_1*std::abs(u_bar - 10.0/std::sqrt(rho_bar))*(u_bar - 10.0/std::sqrt(rho_bar))
                                    + a_2*std::abs(u_bar + 10.0/std::sqrt(rho_bar))*(u_bar + 10.0/std::sqrt(rho_bar)));
                                    */
            double dp_dt = 0.5*(p_1[0]*p_1[0]/rho_1 + p_2[0]*p_2[0]/rho_2)
                           + 50.0*std::log(rho_1/1.0) + 50.0*std::log(rho_2/1.0)
                           - 0.5*(a_1*std::abs(u_bar - 10.0/std::sqrt(rho_bar))*(u_bar - 10.0/std::sqrt(rho_bar))
                                  + a_2*std::abs(u_bar + 10.0/std::sqrt(rho_bar))*(u_bar + 10.0/std::sqrt(rho_bar)));

            //std::cout << "alpha_ext: " << a_1 << " " << a_2 << std::endl;

            //account for pressure term from lhs of element 1
            dp_dt -= 100.0*std::log(rho_0/1.0) + p_0[0]*p_0[0]/rho_0;

            std::cout << "mass_flux[1]       = " << density_fluxes(1) << " ?= " << -dm_dt << " : " << density_fluxes(1) + dm_dt << std::endl;
            std::cout << "momentum_flux[1,0] = " << momentum_fluxes(1,0) << " ?= " << -dp_dt << " : " << momentum_fluxes(1,0) + dp_dt << std::endl;
            std::cout << "momentum_flux[1,1] = " << momentum_fluxes(1,1) << " ?= 0" << std::endl;

            std::cout << "SUM(mass_flux)     = " << density_fluxes.sum() << " ?= 0" << std::endl;


            /*---------------------------------*/
            //Problem 2
            /*---------------------------------*/

            //let p_0 =/= p_1 =/= p_2
            p_0[0] = 2.0;
            p_1[0] = 1.0+std::rand()/(RAND_MAX+1.0);;
            p_2[0] = -std::rand()/(RAND_MAX+1.0);

            //let rho_0 =/= rho_1 =/= rho_2
            rho_0 = 2.0;
            rho_1 = 1.0 + 0.5*std::rand()/(RAND_MAX+1.0);
            rho_2 = 0.5 + 0.5*std::rand()/(RAND_MAX+1.0);

            fluid_body->rho(0) = rho_0;
            fluid_body->rho(1) = rho_1;
            fluid_body->rho(2) = rho_2;

            fluid_body->p[0] = p_0;
            fluid_body->p[1] = p_1;
            fluid_body->p[2] = p_2;

            //call grid to reconstruct conserved fields
            fluid_grid->constructDensityField(job, this);
            fluid_grid->constructMomentumField(job, this);

            //fluid fluxes associated with each volume (to compare against exact solution)
            density_fluxes = fluid_grid->calculateElementMassFluxes(job, this);
            momentum_fluxes = fluid_grid->calculateElementMomentumFluxes(job, this);

            //flux into cell 0 should be zero
            std::cout << "3 Cell Check #2:" << std::endl;
            std::cout << "rho: " << rho_0 << " , " << rho_1 << " , " << rho_2 << std::endl;
            std::cout << "p:   " << p_0[0] << " , " << p_1[0] << " , " << p_2[0] << std::endl;

            double rho_x = (rho_2 - rho_0)/2.0;
            if (std::abs(2*(rho_1 - rho_0)) < std::abs(rho_x)){
                rho_x = 2*(rho_1 - rho_0);
            }
            if (std::abs(2*(rho_2 - rho_1)) < std::abs(rho_x)){
                rho_x = 2*(rho_2 - rho_1);
            }

            double p_x = (p_2[0] - p_0[0])/2.0;
            if (std::abs(2*(p_1[0] - p_0[0])) < std::abs(p_x)){
                p_x = 2*(p_1[0] - p_0[0]);
            }
            if (std::abs(2*(p_2[0] - p_1[0])) < std::abs(p_x)){
                p_x = 2*(p_2[0] - p_1[0]);
            }

            double rho_minus = rho_0;
            double rho_plus = rho_1 - rho_x/2.0;
            double p_minus = p_0[0];
            double p_plus = p_1[0] - p_x/2.0;
            rho_bar = std::sqrt(rho_plus*rho_minus);
            u_bar = (p_plus/std::sqrt(rho_plus) + p_minus/std::sqrt(rho_minus))/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
            double c_bar = 10.0/std::sqrt(rho_bar);
            a_1 = 0.5*(rho_plus - rho_minus) - 1/(2.0*c_bar)*(p_plus - p_minus - (rho_plus - rho_minus)*u_bar);
            a_2 = rho_plus - rho_minus - a_1;

            dm_dt = 0.5*(p_plus + p_minus - a_1*std::abs(u_bar - c_bar) - a_2*std::abs(u_bar+c_bar));
            /*
            dp_dt = 0.5*(p_plus*p_plus/rho_plus + p_minus*p_minus/rho_minus) + 100.0*std::log(rho_bar)
                    + 0.5*(100.0/rho_bar*(rho_plus + rho_minus - 2*rho_bar))
                    - 0.5*(a_1*std::abs(u_bar - c_bar)*(u_bar - c_bar) + a_2*std::abs(u_bar+c_bar)*(u_bar+c_bar))
                    - 100.0*std::log(rho_0);
                    */
            dp_dt = 0.5*(p_plus*p_plus/rho_plus + p_minus*p_minus/rho_minus) + 50.0*std::log(rho_plus) + 50.0*std::log(rho_minus)
                    - 0.5*(a_1*std::abs(u_bar - c_bar)*(u_bar - c_bar) + a_2*std::abs(u_bar+c_bar)*(u_bar+c_bar))
                    - 100.0*std::log(rho_0);

            std::cout << "mass_flux[0]     = " << density_fluxes(0) << " ?= " << -dm_dt << " : " << density_fluxes(0) + dm_dt << std::endl;
            std::cout << "momentum_flux[0] = <" << momentum_fluxes(0,0) << ", " << momentum_fluxes(0,1) << "> ?= <";
            std::cout << -dp_dt << ", 0> : " << momentum_fluxes(0,0)+dp_dt << std::endl;

            std::cout << "SUM(mass_flux)     = " << density_fluxes.sum() << " ?= 0" << std::endl;

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
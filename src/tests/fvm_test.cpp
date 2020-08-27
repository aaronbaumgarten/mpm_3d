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
#include "objects/bodies/bodies.hpp"
#include "objects/nodes/nodes.hpp"
#include "objects/points/points.hpp"
#include "objects/serializers/serializers.hpp"
#include "objects/boundaries/boundaries.hpp"
#include "objects/materials/materials.hpp"

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

    class FVMMPMDragTestDriver : public FiniteVolumeDriver {
    public:
        FVMMPMDragTestDriver(){
            object_name = "FVMMPMDragTestDriver"; //set object name here
        }

        int GRID_DIM;
        int N_i = 101;
        int N_e = 100;
        int N_q = 100000;
        double rho_s = 1000;

        //functions which must be implemented by every driver
        void init(Job* job){
            //force job properties
            job->t = 0;
            job->JOB_TYPE = job->JOB_1D;

            GRID_DIM = 1;

            //set simulation ORDER
            ORDER = 2;

            //set simulation TYPE
            TYPE = FiniteVolumeDriver::THERMAL;

            //set gravity
            gravity = KinematicVector(job->JOB_TYPE);
            gravity.setZero();

            //assign body
            fluid_body = std::unique_ptr<FiniteVolumeBody>(new FVMDefaultBody);
            fluid_body->fp64_props = {1.177, 298}; //rho and theta

            //assign material
            fluid_material = std::unique_ptr<FiniteVolumeMaterial>(new FVMSlurryGasPhase);
            fluid_material->fp64_props = {1.4, 287.1, 1, 0, rho_s, 1}; //gamma, R, eta, k, solid_rho, grain_diam
            fluid_material->int_props = {0}; //solid body id (hopefully doesn't cause issue)

            //assign solver
            solver = std::unique_ptr<FiniteVolumeSolver>(new FVMDefaultSolver);

            //set time step for stability
            job->dt = 0.001; //< dx*sqrt(rho/kappa)
            stop_time = 1.0; //seconds

            return;
        }

        double get_solid_velocity(double x){
            return 1.0 + 0.5*(x - 0.75)*(x - 0.75);
        }

        double get_fluid_velocity(double x){
            return 0.9 - 0.5*(x - 0.25)*(x - 0.25);
        }

        double get_porosity(double x){
            return 0.9 - 0.5*x*x;
        }

        double get_fluid_density(double x){
            return 1.0 + 10.0*(x - 0.5)*(x - 0.5);
        }

        double get_drag(double x){
            double v_s = get_solid_velocity(x);
            double v_f = get_fluid_velocity(x);
            double rho_f = get_fluid_density(x);
            double n = get_porosity(x);
            //let eta and d be 1.0
            double Re = (v_s - v_f)*rho_f;
            return 18*n*(1-n)*(10 * (1-n)/(n*n)
                               + n*n*(1 + 1.5*std::sqrt(1-n))
                               + 0.413*Re/(24*n*n)
                                 * (1/n + 3*n*(1-n) + 8.4*std::pow(Re,-0.343))
                                 /(1 + std::pow(10,3*(1-n))*std::pow(Re,-0.5*(1+4*(1-n)))))*(v_s-v_f);
        }

        void run(Job* job){
            //do nothing
            std::vector<double> L1_err = get_L1_errors(job, 11, 10);
            return;
        }

        std::vector<double> get_L1_errors(Job* job, int n_i, int n_e, int n_q = 100000, std::string str_type = "CartesianLinear"){

            //assign grid dimensions
            N_i = n_i;
            N_e = n_e;
            N_q = n_q;

            //set up mpm grid
            if (str_type == "CartesianCubic"){
                job->grid = std::unique_ptr<Grid>(new CartesianCubic);
            } else if (str_type == "Linear1DNonUniform"){
                job->grid = std::unique_ptr<Grid>(new Linear1DNonUniform);
            } else {
                job->grid = std::unique_ptr<Grid>(new CartesianLinear);
            }
            job->grid->node_count = N_i;
            job->grid->GRID_DIM = 1;
            job->grid->fp64_props = {1.0};
            job->grid->int_props = {N_i-1};

            //set up mpm body
            job->bodies.push_back(std::unique_ptr<Body>(new DefaultBody));
            job->bodies[0]->points = std::unique_ptr<Points>(new DefaultPoints);
            job->bodies[0]->nodes = std::unique_ptr<Nodes>(new DefaultNodes);
            job->bodies[0]->nodes->m = Eigen::VectorXd(N_i);
            job->bodies[0]->nodes->x_t = KinematicVectorArray(N_i, job->JOB_1D);

            //assign grid
            if (str_type == "Linear1DNonUniform"){
                fluid_grid = std::unique_ptr<FiniteVolumeGrid>(new FVMLinear1DNonUniform);
            } else {
                fluid_grid = std::unique_ptr<FiniteVolumeGrid>(new FVMCartesian);
            }
            fluid_grid->fp64_props = {1.0}; //Lx
            fluid_grid->int_props = {N_e};  //Nx
            fluid_grid->str_props = {"USE_LOCAL_POROSITY_CORRECTION"};

            if (true){ //N_e < N_i){
                fluid_grid->fp64_props.push_back(0.1/N_q);
                fluid_grid->str_props.push_back("USE_ENHANCED_M_MATRIX");
            }

            //initialize FVM objects
            job->grid->init(job);
            fluid_grid->init(job, this);
            fluid_material->init(job, this);
            fluid_body->init(job, this);
            solver->init(job, this);

            //first check number of elements
            std::cout << "Number of Elements: " << fluid_grid->element_count << std::endl;

            //check number of nodes
            std::cout << "Number of Nodes: " << job->grid->node_count << std::endl;

            //assign fluid fields
            for (int e=0; e<fluid_grid->element_count;e++){
                //zero out density and energy
                fluid_body->rho(e) = 0;
                fluid_body->p[e].setZero();

                //integrate average field values
                for (int q=0; q<fluid_grid->qpe; q++){
                    fluid_body->rho(e) += fluid_grid->getQuadratureWeight(e*fluid_grid->qpe + q)
                                          *get_fluid_density(fluid_grid->getQuadraturePosition(job, e*fluid_grid->qpe + q)[0]);
                    fluid_body->p(e,0) += fluid_grid->getQuadratureWeight(e*fluid_grid->qpe + q)
                                          *get_fluid_density(fluid_grid->getQuadraturePosition(job, e*fluid_grid->qpe + q)[0])
                                          *get_fluid_velocity(fluid_grid->getQuadraturePosition(job, e*fluid_grid->qpe + q)[0]);
                }

                //average per volume
                fluid_body->rho(e) /= fluid_grid->getElementVolume(e);
                fluid_body->p(e,0) /= fluid_grid->getElementVolume(e);
            }

            //assign solid fields
            for (int i=0; i<job->grid->node_count; i++){
                job->bodies[0]->nodes->m(i) = job->grid->nodeVolume(job, i) * rho_s * (1.0 - get_porosity(job->grid->nodeIDToPosition(job,i)[0]));
                job->bodies[0]->nodes->x_t(i,0) = get_solid_velocity(job->grid->nodeIDToPosition(job,i)[0]);
            }

            //call grid to reconstruct conserved fields
            fluid_grid->mapMixturePropertiesToQuadraturePoints(job, this);
            fluid_grid->constructDensityField(job, this);
            fluid_grid->constructMomentumField(job, this);
            fluid_grid->constructEnergyField(job, this);
            fluid_grid->constructPorosityField(job, this);

            //initialize drag output
            KinematicVectorArray f_d = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
            KinematicVectorArray f_d_e = KinematicVectorArray(fluid_grid->element_count, job->JOB_TYPE);

            //get drag coefficients (should be ok with dt=0)
            job->dt = 0;
            Eigen::VectorXd K_n = fluid_grid->getCorrectedDragCoefficients(job, this);

            //get drag estimate for f_i^drag and F_e^drag
            fluid_grid->calculateSplitIntegralCorrectedDragForces(job, this, f_d, f_d_e, K_n);

            double f_d_nodes = 0;
            double f_d_elements = 0;
            double c = 0; //for safe summation
            double y = 0; //for safe summation
            double t = 0; //for safe summation

            //f_d * V_i = -f_i^drag
            for (int i=0; i<job->grid->node_count; i++) {
                f_d_nodes -= f_d(i,0) * job->grid->nodeVolume(job,i);
                //y = -f_d(i,0) * job->grid->nodeVolume(job,i) - c;
                //t = f_d_nodes + y;
                //c = (t - f_d_nodes) - y;
                //f_d_nodes = t;
            }

            c = 0;
            //adjust F_e^drag to integral form
            for (int e=0; e<fluid_grid->element_count; e++){
                f_d_elements -= f_d_e(e,0) * fluid_grid->getElementVolume(e);
                //y = -f_d_e(e,0) * fluid_grid->getElementVolume(e) - c;
                //t = f_d_elements + y;
                //c = (t - f_d_elements) - y;
                //f_d_elements = t;
            }

            //estimate exact drag with quadrature
            double W_q = 1.0/N_q;
            double X_q = W_q/2.0;
            double f_d_exact = 0;
            c = 0;
            for (int q=0; q<N_q; q++){
                f_d_exact += W_q * get_drag(X_q);
                //y = W_q * get_drag(X_q) - c;
                //t = f_d_exact + y;
                //c = (t - f_d_exact) - y;
                //f_d_exact = t;
                X_q += W_q;
            }

            //save guess
            double f_d_exact_0 = f_d_exact;

            //recalculate
            N_q *= 10;
            W_q = 1.0/N_q;
            X_q = W_q/2.0;
            f_d_exact = 0;
            c = 0;
            for (int q=0; q<N_q; q++){
                f_d_exact += W_q * get_drag(X_q);
                //y = W_q * get_drag(X_q) - c;
                //t = f_d_exact + y;
                //c = (t - f_d_exact) - y;
                //f_d_exact = t;
                X_q += W_q;
            }

            //check answer against exact value
            //int_0^1(f_d) = 9.488652696271 +/- 8e-12
            //std::cout << "Exact Drag Integral: 9.488652696271" << std::endl;
            std::cout << "Exact Drag Integral: " << f_d_exact << " : " << (f_d_exact - f_d_exact_0) << std::endl;
            std::cout << "Nodal Drag Integral: " << f_d_nodes << " : " << (f_d_nodes - f_d_exact) << std::endl;
            std::cout << "Elem. Drag Integral: " << f_d_elements << " : " << (f_d_elements - f_d_exact) << std::endl;

            //calculate mapping for mpm grid and fvm elements
            Eigen::VectorXd f_i_mpm = Eigen::VectorXd(job->grid->node_count);
            f_i_mpm.setZero();

            Eigen::VectorXd f_e_fvm = Eigen::VectorXd(fluid_grid->element_count);
            f_e_fvm.setZero();

            std::vector<int> nvec(0);
            std::vector<double> valvec(0);
            KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
            W_q = 1.0/N_q;
            X_q = W_q/2.0;
            for (int q=0; q<N_q; q++) {
                nvec.resize(0);
                valvec.resize(0);
                tmpVec(0) = X_q;
                job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
                for (int j = 0; j < nvec.size(); j++) {
                    f_i_mpm(nvec[j]) -= W_q*valvec[j]*get_drag(X_q);
                }

                //get element id
                int e = fluid_grid->whichElement(job, tmpVec);
                if (e > -1) {
                    f_e_fvm(e) -= W_q * get_drag(X_q);
                }

                X_q += W_q;
            }

            //calculate f_i L1 norm
            double f_i_L1 = 0;
            for (int i=0; i<job->grid->node_count;i++){
                f_i_L1 += std::abs(f_i_mpm(i) - f_d(i,0) * job->grid->nodeVolume(job,i));
            }

            std::cout << "Nodal Drag Error (L1): " << f_i_L1 << std::endl;


            //calculate f_e L1 norm
            double f_e_L1 = 0;
            for (int e=0; e<fluid_grid->element_count; e++){
                f_e_L1 += std::abs(f_e_fvm(e) - f_d_e(e,0) * fluid_grid->getElementVolume(e));
            }

            std::cout << "Elem. Drag Error (L1): " << f_e_L1 << std::endl;

            std::vector<double> result = {f_i_L1, f_e_L1};

            return result;
        }
    };

    class FVMMPMBuoyancyTestDriver : public FiniteVolumeDriver {
    public:
        FVMMPMBuoyancyTestDriver(){
            object_name = "FVMMPMBuoyancyTestDriver"; //set object name here
        }

        int GRID_DIM;
        int N_i = 101;
        int N_e = 100;
        int N_q = 100000;
        double rho_s = 1000;

        //functions which must be implemented by every driver
        void init(Job* job){
            //force job properties
            job->t = 0;
            job->JOB_TYPE = job->JOB_1D;

            GRID_DIM = 1;

            //set simulation ORDER
            ORDER = 2;

            //set simulation TYPE
            TYPE = FiniteVolumeDriver::THERMAL;

            //set gravity
            gravity = KinematicVector(job->JOB_TYPE);
            gravity.setZero();

            //assign body
            fluid_body = std::unique_ptr<FiniteVolumeBody>(new FVMDefaultBody);
            fluid_body->fp64_props = {1.177, 298}; //rho and theta

            //assign material
            fluid_material = std::unique_ptr<FiniteVolumeMaterial>(new FVMSlurryGasPhase);
            fluid_material->fp64_props = {1.4, 287.1, 1, 0, rho_s, 1}; //gamma, R, eta, k, solid_rho, grain_diam
            fluid_material->int_props = {0}; //solid body id (hopefully doesn't cause issue)

            //assign solver
            solver = std::unique_ptr<FiniteVolumeSolver>(new FVMDefaultSolver);

            //set time step for stability
            job->dt = 0.001; //< dx*sqrt(rho/kappa)
            stop_time = 1.0; //seconds

            return;
        }

        double get_solid_velocity(double x){
            return 1.0 + 0.5*(x - 0.75)*(x - 0.75);
        }

        double get_fluid_velocity(double x){
            return 0.9 - 0.5*(x - 0.25)*(x - 0.25);
        }

        double get_porosity(double x){
            return 0.9 - 0.5*x*x;
        }

        double get_porosity_gradient(double x){
            return -x;
        }

        double get_fluid_density(double x){
            return 1.0 + 10.0*(x - 0.5)*(x - 0.5);
        }

        double get_pressure(double x){
            return 100000.0 + 10000.0*(x - 0.6)*(x - 0.6);
        }

        double get_energy(double x){
            return get_porosity(x)*get_pressure(x)/(0.4) + 0.5*get_fluid_density(x)*get_fluid_velocity(x)*get_fluid_velocity(x);
        }

        double get_buoyancy(double x){
            return get_pressure(x)*get_porosity_gradient(x);
        }

        void run(Job* job){
            //do nothing
            std::vector<double> L1_err = get_L1_errors(job, 10241, 10);
            return;
        }

        std::vector<double> get_L1_errors(Job* job, int n_i, int n_e, int n_q = 100000, std::string str_type = "CartesianLinear"){

            //assign grid dimensions
            N_i = n_i;
            N_e = n_e;
            N_q = n_q;

            //set up mpm grid
            if (str_type == "CartesianCubic"){
                job->grid = std::unique_ptr<Grid>(new CartesianCubic);
            } else if (str_type == "Linear1DNonUniform"){
                job->grid = std::unique_ptr<Grid>(new Linear1DNonUniform);
            } else {
                job->grid = std::unique_ptr<Grid>(new CartesianLinear);
            }
            job->grid->node_count = N_i;
            job->grid->GRID_DIM = 1;
            job->grid->fp64_props = {1.0};
            job->grid->int_props = {N_i-1};

            //set up mpm body
            job->bodies.push_back(std::unique_ptr<Body>(new DefaultBody));
            job->bodies[0]->points = std::unique_ptr<Points>(new DefaultPoints);
            job->bodies[0]->nodes = std::unique_ptr<Nodes>(new DefaultNodes);
            job->bodies[0]->nodes->m = Eigen::VectorXd(N_i);
            job->bodies[0]->nodes->x_t = KinematicVectorArray(N_i, job->JOB_1D);

            //assign grid
            if (str_type == "Linear1DNonUniform"){
                fluid_grid = std::unique_ptr<FiniteVolumeGrid>(new FVMLinear1DNonUniform);
            } else {
                fluid_grid = std::unique_ptr<FiniteVolumeGrid>(new FVMCartesian);
            }
            fluid_grid->fp64_props = {1.0}; //Lx
            fluid_grid->int_props = {N_e};  //Nx
            fluid_grid->str_props = {"USE_LOCAL_POROSITY_CORRECTION"};

            if (true){ //N_e < N_i){
                fluid_grid->fp64_props.push_back(0.1/N_q);
                fluid_grid->str_props.emplace_back("USE_ENHANCED_M_MATRIX");
                //fluid_grid->str_props.emplace_back("USE_GAUSS_GREEN");
            }

            //initialize FVM objects
            job->grid->init(job);
            fluid_grid->init(job, this);
            fluid_material->init(job, this);
            fluid_body->init(job, this);
            solver->init(job, this);

            //first check number of elements
            std::cout << "Number of Elements: " << fluid_grid->element_count << std::endl;

            //check number of nodes
            std::cout << "Number of Nodes: " << job->grid->node_count << std::endl;

            //assign fluid fields
            for (int e=0; e<fluid_grid->element_count;e++){
                //zero out density and energy
                fluid_body->rho(e) = 0;
                fluid_body->p[e].setZero();
                fluid_body->rhoE(e) = 0;

                //integrate average field values
                for (int q=0; q<fluid_grid->qpe; q++){
                    fluid_body->rho(e) += fluid_grid->getQuadratureWeight(e*fluid_grid->qpe + q)
                                          *get_fluid_density(fluid_grid->getQuadraturePosition(job, e*fluid_grid->qpe + q)[0]);
                    fluid_body->p(e,0) += fluid_grid->getQuadratureWeight(e*fluid_grid->qpe + q)
                                          *get_fluid_density(fluid_grid->getQuadraturePosition(job, e*fluid_grid->qpe + q)[0])
                                          *get_fluid_velocity(fluid_grid->getQuadraturePosition(job, e*fluid_grid->qpe + q)[0]);
                    fluid_body->rhoE(e) += fluid_grid->getQuadratureWeight(e*fluid_grid->qpe + q)
                                             *get_energy(fluid_grid->getQuadraturePosition(job, e*fluid_grid->qpe + q)[0]);
                }

                //average per volume
                fluid_body->rho(e) /= fluid_grid->getElementVolume(e);
                fluid_body->p(e,0) /= fluid_grid->getElementVolume(e);
                fluid_body->rhoE(e) /= fluid_grid->getElementVolume(e);
            }

            //assign solid fields
            for (int i=0; i<job->grid->node_count; i++){
                job->bodies[0]->nodes->m(i) = job->grid->nodeVolume(job, i) * rho_s * (1.0 - get_porosity(job->grid->nodeIDToPosition(job,i)[0]));
                job->bodies[0]->nodes->x_t(i,0) = get_solid_velocity(job->grid->nodeIDToPosition(job,i)[0]);
            }

            //call grid to reconstruct conserved fields
            fluid_grid->mapMixturePropertiesToQuadraturePoints(job, this);
            fluid_grid->constructDensityField(job, this);
            fluid_grid->constructMomentumField(job, this);
            fluid_grid->constructEnergyField(job, this);
            fluid_grid->constructPorosityField(job, this);

            //find and report max pressure
            double max_pressure = 0;
            for (int e=0; e<fluid_grid->element_count; e++){
                if (fluid_material->getPressure(job,
                                                this,
                                                fluid_body->rho(e),
                                                fluid_body->p[e],
                                                fluid_body->rhoE(e),
                                                get_porosity(fluid_grid->getElementCentroid(job, e)[0])) > max_pressure){
                    max_pressure = fluid_material->getPressure(job,
                                                               this,
                                                               fluid_body->rho(e),
                                                               fluid_body->p[e],
                                                               fluid_body->rhoE(e),
                                                               get_porosity(fluid_grid->getElementCentroid(job, e)[0]));
                }
            }
            std::cout << "Max Pressure: " << max_pressure << " ?= 103600" << std::endl;

            //find and report max porosity gradient
            KinematicVectorArray grad_n = KinematicVectorArray(fluid_grid->int_quad_count + fluid_grid->ext_quad_count, job->JOB_TYPE);
            grad_n = fluid_grid->gradQ.operate(fluid_body->n,MPMScalarSparseMatrix::TRANSPOSED);
            double max_abs = 0;
            double max_val = 0;
            for (int q=0; q<fluid_grid->int_quad_count; q++){
                if (std::abs(grad_n(q,0)) > max_abs){
                    max_abs = std::abs(grad_n(q,0));
                    max_val = grad_n(q,0);
                }
            }
            std::cout << "Max Porosity Gradient: " << max_val << " ?= -1.0" << std::endl;

            //initialize buoyancy output
            KinematicVectorArray f_b = KinematicVectorArray(job->grid->node_count, job->JOB_TYPE);
            KinematicVectorArray f_b_e = KinematicVectorArray(fluid_grid->element_count, job->JOB_TYPE);

            //get drag estimate for f_i^drag and F_e^drag
            fluid_grid->calculateSplitIntegralBuoyantForces(job, this, f_b, f_b_e);

            //remove edge contributions
            double rho, rhoE, n, P;
            KinematicVector p = KinematicVector(job->JOB_TYPE);
            //approximate local density, momentum, and energy
            rho = fluid_body->rho(0) + fluid_body->rho_x(0,0)*(0.0 - fluid_grid->getElementCentroid(job,0)[0]);
            rhoE = fluid_body->rhoE(0) + fluid_body->rhoE_x(0,0)*(0.0 - fluid_grid->getElementCentroid(job,0)[0]);
            p = fluid_body->p[0];
            p(0) += fluid_body->p_x(0,0,0)*(0.0 - fluid_grid->getElementCentroid(job,0)[0]);
            n = fluid_body->n_e(0);

            //fill in pressure
            P = fluid_material->getPressure(job,
                                            this,
                                            rho,
                                            p,
                                            rhoE,
                                            n); //n_q(e*qpe + q));

            f_b(0,0) -= (1.0 - get_porosity(0.0))*P / job->grid->nodeVolume(job,0);

            //approximate local density, momentum, and energy
            rho = fluid_body->rho(N_e-1) + fluid_body->rho_x(N_e-1,0)*(1.0 - fluid_grid->getElementCentroid(job,N_e-1)[0]);
            rhoE = fluid_body->rhoE(N_e-1) + fluid_body->rhoE_x(N_e-1,0)*(1.0 - fluid_grid->getElementCentroid(job,N_e-1)[0]);
            p = fluid_body->p[N_e-1];
            p(0) += fluid_body->p_x(N_e-1,0,0)*(1.0 - fluid_grid->getElementCentroid(job,N_e-1)[0]);
            n = fluid_body->n_e(N_e-1);

            //fill in pressure
            P = fluid_material->getPressure(job,
                                            this,
                                            rho,
                                            p,
                                            rhoE,
                                            n); //n_q(e*qpe + q));

            f_b(f_b.size()-1,0) += (1.0 - get_porosity(1.0))*P / job->grid->nodeVolume(job,f_b.size()-1);

            double f_b_nodes = 0;
            double f_b_elements = 0;
            double c = 0; //for safe summation
            double y = 0; //for safe summation
            double t = 0; //for safe summation

            //f_d * V_i = -f_i^drag
            for (int i=0; i<job->grid->node_count; i++) {
                f_b_nodes -= f_b(i,0) * job->grid->nodeVolume(job,i);
                //y = -f_d(i,0) * job->grid->nodeVolume(job,i) - c;
                //t = f_d_nodes + y;
                //c = (t - f_d_nodes) - y;
                //f_d_nodes = t;
            }

            c = 0;
            //adjust F_e^drag to integral form
            for (int e=0; e<fluid_grid->element_count; e++){
                f_b_elements -= f_b_e(e,0) * fluid_grid->getElementVolume(e);
                //y = -f_d_e(e,0) * fluid_grid->getElementVolume(e) - c;
                //t = f_d_elements + y;
                //c = (t - f_d_elements) - y;
                //f_d_elements = t;
            }

            //estimate exact drag with quadrature
            double W_q = 1.0/N_q;
            double X_q = W_q/2.0;
            double f_b_exact = 0;
            c = 0;
            for (int q=0; q<N_q; q++){
                f_b_exact += W_q * get_buoyancy(X_q);
                //y = W_q * get_drag(X_q) - c;
                //t = f_d_exact + y;
                //c = (t - f_d_exact) - y;
                //f_d_exact = t;
                X_q += W_q;
            }

            //save guess
            double f_b_exact_0 = f_b_exact;

            //recalculate
            N_q *= 10;
            W_q = 1.0/N_q;
            X_q = W_q/2.0;
            f_b_exact = 0;
            c = 0;
            for (int q=0; q<N_q; q++){
                f_b_exact += W_q * get_buoyancy(X_q);
                //y = W_q * get_drag(X_q) - c;
                //t = f_d_exact + y;
                //c = (t - f_d_exact) - y;
                //f_d_exact = t;
                X_q += W_q;
            }

            //check answer against exact value
            //int_0^1(f_d) = 9.488652696271 +/- 8e-12
            //std::cout << "Exact Drag Integral: 9.488652696271" << std::endl;
            std::cout << "Exact Buoy Integral: " << f_b_exact << " : " << (f_b_exact - f_b_exact_0) << std::endl;
            std::cout << "Nodal Buoy Integral: " << f_b_nodes << " : " << (f_b_nodes - f_b_exact) << std::endl;
            std::cout << "Elem. Buoy Integral: " << f_b_elements << " : " << (f_b_elements - f_b_exact) << std::endl;

            //calculate mapping for mpm grid and fvm elements
            Eigen::VectorXd f_i_mpm = Eigen::VectorXd(job->grid->node_count);
            f_i_mpm.setZero();

            Eigen::VectorXd f_e_fvm = Eigen::VectorXd(fluid_grid->element_count);
            f_e_fvm.setZero();

            std::vector<int> nvec(0);
            KinematicVectorArray gradvec(0,job->JOB_TYPE);
            KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
            double p_err = 0;
            double rho_err = 0;
            W_q = 1.0/N_q;
            X_q = W_q/2.0;
            for (int q=0; q<N_q; q++) {
                nvec.resize(0);
                gradvec.resize(0);
                tmpVec(0) = X_q;
                job->grid->evaluateBasisFnGradient(job, tmpVec, nvec, gradvec);
                for (int j = 0; j < nvec.size(); j++) {
                    f_i_mpm(nvec[j]) += W_q*gradvec(j,0)*(1.0 - fluid_body->n(nvec[j]))*get_pressure(X_q);
                }

                //get element id
                int e = fluid_grid->whichElement(job, tmpVec);
                if (e > -1) {
                    f_e_fvm(e) -= W_q * get_buoyancy(X_q);


                    //approximate local density, momentum, and energy
                    rho = fluid_body->rho(e) + fluid_body->rho_x(e,0)*(X_q - fluid_grid->getElementCentroid(job,e)[0]);
                    rhoE = fluid_body->rhoE(e) + fluid_body->rhoE_x(e,0)*(X_q - fluid_grid->getElementCentroid(job,e)[0]);
                    p = fluid_body->p[e];
                    p(0) += fluid_body->p_x(e,0,0)*(X_q - fluid_grid->getElementCentroid(job,e)[0]);
                    n = fluid_body->n_e(e);

                    //fill in pressure
                    P = fluid_material->getPressure(job,
                                                    this,
                                                    rho,
                                                    p,
                                                    rhoE,
                                                    n); //n_q(e*qpe + q));

                    if (std::abs(get_pressure(X_q) - P) > p_err){
                        p_err = std::abs(get_pressure(X_q) - P);
                    }

                    if (std::abs(get_fluid_density(X_q)/get_porosity(X_q) - rho/n) > rho_err){
                        rho_err = std::abs(get_fluid_density(X_q)/get_porosity(X_q) - rho/n);
                    }
                }

                X_q += W_q;
            }
            std::cout << "Pressure Error: " << p_err << std::endl;
            std::cout << "Density Error: " << rho_err << std::endl;

            //calculate f_i L1 norm (without edges)
            double f_i_L1 = 0;
            for (int i=0; i<job->grid->node_count;i++){
                f_i_L1 += std::abs(f_i_mpm(i) - f_b(i,0) * job->grid->nodeVolume(job,i));
                //std::cout << f_i_mpm(i) << " ?= " << f_b(i,0) * job->grid->nodeVolume(job,i) << std::endl;
            }

            std::cout << "Nodal Drag Error (L1): " << f_i_L1 << std::endl;


            //calculate f_e L1 norm
            double f_e_L1 = 0;
            for (int e=0; e<fluid_grid->element_count; e++){
                f_e_L1 += std::abs(f_e_fvm(e) - f_b_e(e,0) * fluid_grid->getElementVolume(e));
            }

            std::cout << "Elem. Drag Error (L1): " << f_e_L1 << std::endl;

            std::vector<double> result = {f_i_L1, f_e_L1};

            return result;
        }
    };

    class FVMMPMPorosityTestDriver : public FiniteVolumeDriver {
    public:
        FVMMPMPorosityTestDriver(){
            object_name = "FVMMPMPorosityTestDriver"; //set object name here
        }

        int GRID_DIM;
        int N_i = 101;
        int N_e = 100;
        int N_p = 100;
        int N_q = 100000;
        double rho_s = 1000;

        //functions which must be implemented by every driver
        void init(Job* job){
            //force job properties
            job->t = 0;
            job->JOB_TYPE = job->JOB_1D;

            GRID_DIM = 1;

            //set simulation ORDER
            ORDER = 2;

            //set simulation TYPE
            TYPE = FiniteVolumeDriver::THERMAL;

            //set gravity
            gravity = KinematicVector(job->JOB_TYPE);
            gravity.setZero();

            //assign body
            fluid_body = std::unique_ptr<FiniteVolumeBody>(new FVMDefaultBody);
            fluid_body->fp64_props = {1.177, 298}; //rho and theta

            //assign material
            fluid_material = std::unique_ptr<FiniteVolumeMaterial>(new FVMSlurryGasPhase);
            fluid_material->fp64_props = {1.4, 287.1, 1, 0, rho_s, 1}; //gamma, R, eta, k, solid_rho, grain_diam
            fluid_material->int_props = {0}; //solid body id (hopefully doesn't cause issue)

            //assign solver
            solver = std::unique_ptr<FiniteVolumeSolver>(new FVMDefaultSolver);

            //set time step for stability
            job->dt = 0.001; //< dx*sqrt(rho/kappa)
            stop_time = 1.0; //seconds

            return;
        }

        double get_porosity(double x){
            return 0.9 - 0.5*x*x;
        }

        double get_displacement(double x){
            return 0.75*(x*x - x);
        }

        double get_displaced_position(double x){
            return 0.25*x + 0.75*x*x;
        }

        double get_original_position(double y){
            double a = 0.75;
            double b = 0.25;
            double c = -y;
            return (-b + std::sqrt(b*b - 4.0*a*c))/(2.0*a);
        }

        int which_point(double x, int N_p){
            //convert position into undeformed space
            double y = get_original_position(x);
            return (int)std::floor(y * N_p);
        }

        void run(Job* job){
            //do nothing
            double L1_err = get_L1_error(job, 10241, 10);
            return;
        }

        double get_L1_error(Job* job, int n_i, int n_p, int n_q = 100000, std::string str_type = "CartesianLinear"){

            //assign grid dimensions
            N_i = n_i;
            N_e = n_p; //for simplicity, collocate point centroids with FVM cell centers
            N_p = n_p;
            N_q = n_q;

            //set up mpm grid
            if (str_type == "CartesianCubic"){
                job->grid = std::unique_ptr<Grid>(new CartesianCubic);
            } else if (str_type == "Linear1DNonUniform"){
                job->grid = std::unique_ptr<Grid>(new Linear1DNonUniform);
            } else {
                job->grid = std::unique_ptr<Grid>(new CartesianLinear);
            }
            job->grid->node_count = N_i;
            job->grid->GRID_DIM = 1;
            job->grid->fp64_props = {1.0};
            job->grid->int_props = {N_i-1};

            //set up mpm body
            job->bodies.push_back(std::unique_ptr<Body>(new DefaultBody));
            job->bodies[0]->points = std::unique_ptr<Points>(new DefaultPoints);
            job->bodies[0]->nodes = std::unique_ptr<Nodes>(new DefaultNodes);
            job->bodies[0]->nodes->m = Eigen::VectorXd(N_i);
            job->bodies[0]->nodes->x = KinematicVectorArray(N_i, job->JOB_1D);
            job->bodies[0]->nodes->x_t = KinematicVectorArray(N_i, job->JOB_1D);

            //need to define x, m, v on points
            job->bodies[0]->points->m = Eigen::VectorXd(N_p);
            job->bodies[0]->points->v = Eigen::VectorXd(N_p);
            job->bodies[0]->points->x = KinematicVectorArray(N_p, job->JOB_1D);
            job->bodies[0]->points->active = Eigen::VectorXi(N_p);

            //assign grid
            fluid_grid = std::unique_ptr<FiniteVolumeGrid>(new FVMCartesian);//(new FVMLinear1DNonUniform);
            fluid_grid->fp64_props = {1.0}; //Lx
            fluid_grid->int_props = {N_e};  //Nx

            //initialize FVM objects
            job->grid->init(job);
            fluid_grid->init(job, this);
            fluid_material->init(job, this);
            fluid_body->init(job, this);
            solver->init(job, this);

            //first check number of elements
            std::cout << "Number of Elements: " << fluid_grid->element_count << std::endl;

            //check number of nodes
            std::cout << "Number of Nodes: " << job->grid->node_count << std::endl;

            //assign fluid fields
            for (int e=0; e<fluid_grid->element_count;e++){
                //zero out density and energy
                fluid_body->rho(e) = 0;
                fluid_body->p[e].setZero();
                fluid_body->rhoE(e) = 0;
            }

            //assign solid fields
            for (int i=0; i<job->grid->node_count; i++){
                job->bodies[0]->nodes->m(i) = job->grid->nodeVolume(job, i) * rho_s * (1.0 - get_porosity(job->grid->nodeIDToPosition(job,i)[0]));
                job->bodies[0]->nodes->x_t(i,0) = 0.0;
            }
            double x_tmp, x0_tmp, x1_tmp;
            for (int p=0; p<job->bodies[0]->points->x.size(); p++){
                job->bodies[0]->points->x[p] = fluid_grid->getElementCentroid(job, p); //equivalence b/w points and fluid grid
                job->bodies[0]->points->v(p) = fluid_grid->getElementVolume(p);
                job->bodies[0]->points->m(p) = rho_s
                                               *job->bodies[0]->points->v(p)
                                               *(1.0 - get_porosity(job->bodies[0]->points->x(p,0)));
                job->bodies[0]->points->active(p) = 1;
            }

            //initialize MPM objects
            job->bodies[0]->init(job);

            //reset point values to displaced state
            for (int p=0; p<job->bodies[0]->points->x.size(); p++){
                x_tmp = fluid_grid->getElementCentroid(job, p)[0];
                x0_tmp = fluid_grid->getFaceCentroid(job, p)[0];
                x1_tmp = fluid_grid->getFaceCentroid(job, p+1)[0];

                job->bodies[0]->points->x(p,0) = get_displaced_position(x_tmp);
                job->bodies[0]->points->v(p) = get_displaced_position(x1_tmp) - get_displaced_position(x0_tmp);
                job->bodies[0]->points->m(p) = rho_s
                                               *job->bodies[0]->points->v(p)
                                               *(1.0 - get_porosity(job->bodies[0]->points->x(p,0)));
                job->bodies[0]->points->active(p) = 1;
            }

            //create mapping matrix \tilde{S}_ip
            if (str_type == "CartesianCubic"){
                job->bodies[0]->generateMap(job, DefaultBody::CPDI_OFF);
            } else {
                job->bodies[0]->generateMap(job, DefaultBody::CPDI_ON);
            }

            //call fluid to construct porosity field
            job->bodies[0]->nodes->m = job->bodies[0]->S*job->bodies[0]->points->m;
            fluid_material->calculatePorosity(job, this);

            //recalculate N_q (for consistency)
            N_q *= 10;

            //calculate exact S_ip and A_ij
            MPMScalarSparseMatrix S = MPMScalarSparseMatrix(N_i,N_p);
            MPMScalarSparseMatrix A = MPMScalarSparseMatrix(N_i,N_i);

            std::vector<int> nvec(0);
            std::vector<double> valvec(0);
            KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
            double W_q = 1.0/N_q;
            double X_q = W_q/2.0;
            for (int q=0; q<N_q; q++) {
                nvec.resize(0);
                valvec.resize(0);
                tmpVec(0) = X_q;
                job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);

                for (int i = 0; i < nvec.size(); i++) {
                    for (int j = 0; j < nvec.size(); j++) {
                        A.push_back(nvec[i],nvec[j],W_q*valvec[i]*valvec[j]);
                    }
                }

                //get point id
                int e = which_point(X_q, N_p); //fluid_grid->whichElement(job, tmpVec);
                if (e > -1) {
                    for (int i=0; i<nvec.size(); i++){
                        S.push_back(nvec[i],e,W_q*valvec[i]/job->bodies[0]->points->v(e));
                    }
                }

                X_q += W_q;
            }

            //calculate error
            Eigen::VectorXd m_j_exact = S*job->bodies[0]->points->m;

            Eigen::VectorXd rho_bar_s = Eigen::VectorXd(job->grid->node_count);
            for (int i=0; i<job->grid->node_count; i++){
                rho_bar_s(i) = rho_s * (1.0 - fluid_body->n(i));
            }

            Eigen::VectorXd m_j_approx = A*rho_bar_s;

            double uh = 0;
            for (int i=0; i<job->grid->node_count; i++){
                uh += rho_bar_s(i) * job->grid->nodeVolume(job,i);
            }

            //calculate f_i L1 norm (without edges)
            double f_i_L1 = 0;
            for (int i=0; i<job->grid->node_count;i++){
                f_i_L1 += std::abs(m_j_exact(i) - m_j_approx(i));
                //std::cout << m_j_exact(i) << " ?= " << m_j_approx(i) << std::endl;
            }

            std::cout << "Total Error: " <<  m_j_approx.sum() << " ... " << m_j_approx.sum() - m_j_exact.sum() << std::endl;
            std::cout << "Nodal Drag Error (L1): " << f_i_L1 << std::endl;

            /*
            Eigen::VectorXd ones = Eigen::VectorXd::Ones(job->bodies[0]->points->x.size());
            Eigen::VectorXd exact = S*ones;
            Eigen::VectorXd approx = S*ones;
            for (int i=0; i<job->grid->node_count; i++){

            }
            std::cout << "8th node error:" << std::endl;
             */

            return f_i_L1;
        }
    };

    class FVMMPMMoMSDriver : public FiniteVolumeDriver {
    public:
        FVMMPMMoMSDriver(){
            object_name = "FVMMPMMoMSDriver"; //set object name here
        }

        int GRID_DIM;
        int N_i = 101;
        int N_e = 100;
        int N_p = 100;
        int N_q = 1;
        double rho_s = 1000;
        double rho_f = 117.7; //1.177;
        double eta_0 = 1.0;
        double d = 0.1;
        double fps = 100;
        double R = 2.871e-2;
        double E = 1e4;
        double nu = 0.3;
        double G = E / (2.0 * (1.0 + nu));

        //functions which must be implemented by every driver
        void init(Job* job){
            //force job properties for 2D simulation
            job->t = 0;
            job->JOB_TYPE = job->JOB_2D;
            job->DIM = 2;

            GRID_DIM = 2;

            //set simulation ORDER
            ORDER = 2;

            //set simulation TYPE
            TYPE = FiniteVolumeDriver::THERMAL;

            //set gravity
            gravity = KinematicVector(job->JOB_TYPE);
            gravity.setZero();

            //assign body
            fluid_body = std::unique_ptr<FiniteVolumeBody>(new FVMDefaultBody);
            fluid_body->fp64_props = {rho_f, 298}; //rho and theta

            //assign material
            fluid_material = std::unique_ptr<FiniteVolumeMaterial>(new FVMSlurryGasPhase);
            fluid_material->fp64_props = {1.4, R, eta_0, 0, rho_s, d}; //gamma, R, eta, k, solid_rho, grain_diam
            fluid_material->int_props = {0}; //solid body id (hopefully doesn't cause issue)

            //assign solver
            solver = std::unique_ptr<FiniteVolumeSolver>(new FVMMixtureSolverRK4);
            solver->str_props = {"sand"};
            solver->int_props = {1, 1};

            //assign serializer
            serializer = std::unique_ptr<FiniteVolumeSerializer>(new FVMDefaultVTK);
            serializer->fp64_props = {fps};
            serializer->str_props = {"output", "test"};

            //set time step for stability
            job->dt = 3e-5; //< dx*sqrt(rho/kappa)
            stop_time = 1.0; //seconds

            //setup serializer
            job->serializer = std::unique_ptr<Serializer>(new DefaultVTK);
            job->serializer->fp64_props = {fps};
            job->serializer->str_props = {"output","output","test"};

            return;
        }

        double get_porosity(const KinematicVector &x){
            return 0.4;

            KinematicVector tmpX = x;
            tmpX[0] = 2.0*x[0] - 1.0;
            tmpX[1] = 2.0*x[1] - 1.0;
            return 0.3 + 0.2*(tmpX[0]*tmpX[0] + tmpX[1]*tmpX[1]);
        }

        KinematicVector get_dn_dx(Job* job, const KinematicVector &x){
            return KinematicVector(job->JOB_TYPE);

            KinematicVector tmpX = x;
            tmpX[0] = 2.0*x[0] - 1.0;
            tmpX[1] = 2.0*x[1] - 1.0;
            KinematicVector tmp = KinematicVector(job->JOB_TYPE);
            tmp[0] = 0.8*tmpX[0];
            tmp[1] = 0.8*tmpX[1];
            return tmp;
        }

        double get_A(double t){
            return -t*t;
        }

        double get_dA_dt(double t){
            return -2*t;
        }

        double get_B(double t){
            return 4.0*t*t;
        }

        double get_dB_dt(double t){
            return 8.0*t;
        }

        double get_int_B(double t){
            return 4.0*t*t*t/3.0;
        }

        KinematicVector get_solid_velocity(Job* job, const KinematicVector &x, double t){
            KinematicVector tmpX = x;
            tmpX[0] = 2.0*x[0] - 1.0;
            tmpX[1] = 2.0*x[1] - 1.0;
            KinematicVector tmp = KinematicVector(job->JOB_TYPE);
            double r2 = tmpX[0]*tmpX[0] + tmpX[1]*tmpX[1];
            if (r2 <= 1.0) {
                tmp[0] = -tmpX[1] * (r2 - 1) * (r2 - 1);
                tmp[1] = tmpX[0] * (r2 - 1) * (r2 - 1);
            }
            return get_B(t)*tmp;
        }

        KinematicVector get_fluid_velocity(Job* job, const KinematicVector &x, double t){
            KinematicVector tmpX = x;
            tmpX[0] = 2.0*x[0] - 1.0;
            tmpX[1] = 2.0*x[1] - 1.0;
            KinematicVector tmp = KinematicVector(job->JOB_TYPE);
            tmp[0] = 4 * tmpX[1] * (tmpX[0]*tmpX[0] - 1) * (tmpX[0] * tmpX[0] - 1) * (tmpX[1]*tmpX[1] - 1);
            tmp[1] = -4 * tmpX[0] * (tmpX[0]*tmpX[0] - 1) * (tmpX[1] * tmpX[1] - 1) * (tmpX[1]*tmpX[1] - 1);
            return get_A(t)*tmp;
        }

        KinematicVector get_drag(Job* job, const KinematicVector &x, double t){
            KinematicVector v_s = get_solid_velocity(job, x, t);
            KinematicVector v_f = get_fluid_velocity(job, x, t);
            double bar_rho_f = get_porosity(x)*rho_f;
            double n = get_porosity(x);
            double Re = d*(v_s - v_f).norm()*bar_rho_f/eta_0;
            if (Re < 1e-10) {
                return eta_0 / (d*d) * 18 * n * (1 - n) * (10 * (1 - n) / (n * n)
                                                           + n * n * (1 + 1.5 * std::sqrt(1 - n)))
                       * (v_s - v_f);
            } else {
                return eta_0 / (d*d) * 18 * n * (1 - n) * (10 * (1 - n) / (n * n)
                                                           + n * n * (1 + 1.5 * std::sqrt(1 - n))
                                                           + 0.413 * Re / (24 * n * n)
                                                             * (1 / n + 3 * n * (1 - n) + 8.4 * std::pow(Re, -0.343))
                                                             / (1 + std::pow(10, 3 * (1 - n)) * std::pow(Re, -0.5 * (1 + 4 * (1 - n)))))
                       * (v_s - v_f);
            }
        }

        KinematicVector get_solid_body_force(Job* job, const KinematicVector &x, double t){
            KinematicVector tmpX = x;
            tmpX[0] = 2.0*x[0] - 1.0;
            tmpX[1] = 2.0*x[1] - 1.0;

            KinematicVector f_d = get_drag(job, x, t);
            double n = get_porosity(x);

            KinematicVector tmp = KinematicVector(job->JOB_TYPE);

            //rho dv/dt
            double r2 = tmpX[0]*tmpX[0] + tmpX[1]*tmpX[1];
            if (r2 <= 1.0) {
                tmp[0] = -(1 - n) * rho_s * get_dB_dt(t) * tmpX[1] * (r2 - 1) * (r2 - 1);
                tmp[1] = (1 - n) * rho_s * get_dB_dt(t) * tmpX[0] * (r2 - 1) * (r2 - 1);
            }

            //f_d
            tmp += f_d;

            //-div(T)
            double D = get_int_B(t);
            if (r2 <= 1.0){
                tmp[0] -= 4.0 * G*(-8.0*D*(3.0*r2 - 2.0)*tmpX[1] - 16.0/3.0 * D * D * (11.0*r2*(r2 - 1)*(r2 - 1.0) + 4.0*r2*(r2 - 1.0))*tmpX[0]);
                tmp[1] -= 4.0 * G*( 8.0*D*(3.0*r2 - 2.0)*tmpX[0] - 16.0/3.0 * D * D * (11.0*r2*(r2 - 1)*(r2 - 1.0) + 4.0*r2*(r2 - 1.0))*tmpX[1]);
            }

            return tmp;
            //pg 99, nb 7
        }

        KinematicVector get_fluid_body_force(Job* job, const KinematicVector &x, double t){
            KinematicVector tmpX = x;
            tmpX[0] = 2.0*x[0] - 1.0;
            tmpX[1] = 2.0*x[1] - 1.0;

            KinematicVector f_d = get_drag(job, x, t);
            double n = get_porosity(x);

            KinematicVector tmp = KinematicVector(job->JOB_TYPE);

            //rho dv/dt
            tmp[0] = n * rho_f * get_dA_dt(t) * 4.0 * (tmpX[0]*tmpX[0] - 1) * (tmpX[0]*tmpX[0] - 1) * tmpX[1] * (tmpX[1]*tmpX[1] - 1);
            tmp[1] = -n * rho_f * get_dA_dt(t) * 4.0 * tmpX[0] * (tmpX[0]*tmpX[0] - 1) * (tmpX[1]*tmpX[1] - 1) * (tmpX[1]*tmpX[1] - 1);

            //f_d
            tmp -= f_d;

            //dn/dx
            KinematicVector gradn = get_dn_dx(job, x);

            //-div(T)
            tmp[0] -= 4.0 * get_A(t) * eta_0 * (1 + 2.5*(1.0-n)) * (16.0 *(tmpX[1]*tmpX[1] - 1)*tmpX[1]*(3.0*tmpX[0]*tmpX[0] - 1) + 24.0*tmpX[1]*(tmpX[0]*tmpX[0] - 1)*(tmpX[0]*tmpX[0] - 1));
            tmp[1] -= 4.0 * get_A(t) * eta_0 * (1 + 2.5*(1.0-n)) * (-16.0 *(tmpX[0]*tmpX[0] - 1)*tmpX[0]*(3.0*tmpX[1]*tmpX[1] - 1) - 24.0*tmpX[0]*(tmpX[1]*tmpX[1] - 1)*(tmpX[1]*tmpX[1] - 1));

            KinematicTensor tmpMat = KinematicTensor(job->JOB_TYPE);
            tmpMat(0,0) = 2.0* 16.0*tmpX[0]*(tmpX[0]*tmpX[0] - 1)*tmpX[1]*(tmpX[1]*tmpX[1]-1);
            tmpMat(0,1) = 2.0* 2.0 * ((tmpX[0]*tmpX[0] - 1)*(tmpX[0]*tmpX[0] - 1)*(3.0*tmpX[1]*tmpX[1] - 1) - (tmpX[1]*tmpX[1] - 1)*(tmpX[1]*tmpX[1] - 1)*(3.0*tmpX[0]*tmpX[0] - 1));
            tmpMat(1,0) = tmpMat(0,1);
            tmpMat(1,1) = -tmpMat(0,0);

            tmp += 5.0/2.0*get_A(t)*2.0*eta_0*tmpMat*gradn;

            //+div(\rho v \otimes v)
            tmp[0] += n * rho_f * 2.0 * 16.0 * get_A(t) * get_A(t)
                      * tmpX[0] * (tmpX[0]*tmpX[0] - 1) * (tmpX[0]*tmpX[0] - 1) * (tmpX[0]*tmpX[0] - 1)
                      * (tmpX[1]*tmpX[1]*tmpX[1]*tmpX[1] - 1) * (tmpX[1]*tmpX[1] - 1);

            tmp[1] += n * rho_f * 2.0 * 16.0 * get_A(t) * get_A(t)
                      * tmpX[1] * (tmpX[1]*tmpX[1] - 1) * (tmpX[1]*tmpX[1] - 1) * (tmpX[1]*tmpX[1] - 1)
                      * (tmpX[0]*tmpX[0]*tmpX[0]*tmpX[0] - 1) * (tmpX[0]*tmpX[0] - 1);

            tmpMat(0,0) = tmpX[1]*tmpX[1]*(tmpX[0]*tmpX[0] - 1)*(tmpX[0]*tmpX[0] - 1);
            tmpMat(0,1) = -tmpX[0]*tmpX[1]*(tmpX[0]*tmpX[0] - 1)*(tmpX[1]*tmpX[1] - 1);
            tmpMat(1,0) = tmpMat(0,1);
            tmpMat(1,1) = tmpX[0]*tmpX[0]*(tmpX[1]*tmpX[1] - 1)*(tmpX[1]*tmpX[1] - 1);
            tmpMat *= 16.0 * get_A(t) * get_A(t) * (tmpX[0]*tmpX[0] - 1) * (tmpX[0]*tmpX[0] - 1) * (tmpX[1]*tmpX[1] - 1) * (tmpX[1]*tmpX[1] - 1);

            tmp += rho_f * tmpMat * gradn;

            return tmp;
            //pg 99, nb 7
        }

        double get_fluid_energy_adjustment(Job* job, const KinematicVector &x, double t){

            KinematicVector tmpX = x;
            tmpX[0] = 2.0*x[0] - 1.0;
            tmpX[1] = 2.0*x[1] - 1.0;

            KinematicVector v_s = get_solid_velocity(job, x, t);
            KinematicVector v_f = get_fluid_velocity(job, x, t);
            double bar_rho_f = get_porosity(x)*rho_f;
            double n = get_porosity(x);
            KinematicVector gradn = get_dn_dx(job, x);
            KinematicVector f_d = get_drag(job, x, t);

            //+0.5 * v_f * div(\rho v \otimes v)
            KinematicTensor tmpMat = KinematicTensor(job->JOB_TYPE);
            KinematicVector tmp = KinematicVector(job->JOB_TYPE);
            tmp[0] = n * rho_f * 2.0 * 16.0 * get_A(t) * get_A(t)
                      * tmpX[0] * (tmpX[0]*tmpX[0] - 1) * (tmpX[0]*tmpX[0] - 1) * (tmpX[0]*tmpX[0] - 1)
                      * (tmpX[1]*tmpX[1]*tmpX[1]*tmpX[1] - 1) * (tmpX[1]*tmpX[1] - 1);

            tmp[1] = n * rho_f * 2.0 * 16.0 * get_A(t) * get_A(t)
                      * tmpX[1] * (tmpX[1]*tmpX[1] - 1) * (tmpX[1]*tmpX[1] - 1) * (tmpX[1]*tmpX[1] - 1)
                      * (tmpX[0]*tmpX[0]*tmpX[0]*tmpX[0] - 1) * (tmpX[0]*tmpX[0] - 1);

            tmpMat(0,0) = tmpX[1]*tmpX[1]*(tmpX[0]*tmpX[0] - 1)*(tmpX[0]*tmpX[0] - 1);
            tmpMat(0,1) = -tmpX[0]*tmpX[1]*(tmpX[0]*tmpX[0] - 1)*(tmpX[1]*tmpX[1] - 1);
            tmpMat(1,0) = tmpMat(0,1);
            tmpMat(1,1) = tmpX[0]*tmpX[0]*(tmpX[1]*tmpX[1] - 1)*(tmpX[1]*tmpX[1] - 1);
            tmpMat *= 16.0 * get_A(t) * get_A(t) * (tmpX[0]*tmpX[0] - 1) * (tmpX[0]*tmpX[0] - 1) * (tmpX[1]*tmpX[1] - 1) * (tmpX[1]*tmpX[1] - 1);

            tmp += rho_f * tmpMat * gradn;

            double p_0 = rho_f*R*298.0;
            return -p_0*gradn.dot(v_s - v_f) + p_0*gradn.dot(v_f)/(0.4) + 0.5*v_f.dot(tmp) - f_d.dot(v_s);
        }

        double get_fluid_mass_adjustment(Job* job, const KinematicVector &x, double t){
            KinematicVector v_s = get_solid_velocity(job, x, t);
            KinematicVector v_f = get_fluid_velocity(job, x, t);
            double bar_rho_f = get_porosity(x)*rho_f;
            double n = get_porosity(x);
            KinematicVector gradn = get_dn_dx(job, x);

            return rho_f*v_f.dot(gradn);
        }

        KinematicVector getFluidLoading(Job *job, const KinematicVector &x){
            return get_fluid_body_force(job, x, job->t);
        }

        KinematicVector getSolidLoading(Job *job, const KinematicVector &x){
            return get_solid_body_force(job, x, job->t);
        }

        void run(Job* job){
            //do nothing
            std::vector<double> L1_err = get_L1_error(job, 11, 10, 10, 1, 1e-5, "test");
            return;
        }

        std::vector<double> get_L1_error(Job* job, int n_i, int n_p, int n_e, int n_q, double dt, std::string sim_name, int num_threads = 20){
            //assign job dt
            job->dt = dt;

            //reset job time
            job->t = 0;

            //assign grid dimensions
            N_i = n_i;
            N_e = n_e;
            N_p = n_p;
            N_q = n_q;

            //set up mpm grid
            job->grid = std::unique_ptr<Grid>(new CartesianLinear);
            job->grid->node_count = N_i*N_i;
            job->grid->GRID_DIM = 2;
            job->grid->fp64_props = {1.0, 1.0};
            job->grid->int_props = {(N_i-1), (N_i-1)};

            //set up mpm body
            job->bodies.clear();
            job->bodies.push_back(std::unique_ptr<Body>(new DefaultBody));
            job->activeBodies = {1};

            //assign simulation name
            job->bodies[0]->name = sim_name; //"sand";
            job->bodies[0]->id = 0;
            solver->str_props = {sim_name};
            serializer->str_props[1] = sim_name;

            job->bodies[0]->points = std::unique_ptr<Points>(new ThreadPoolPoints);

            int len = N_p*N_p;
            job->bodies[0]->points->x = KinematicVectorArray(len, job->JOB_TYPE);
            job->bodies[0]->points->u = KinematicVectorArray(len, job->JOB_TYPE);
            job->bodies[0]->points->x_t = KinematicVectorArray(len, job->JOB_TYPE);
            job->bodies[0]->points->mx_t = KinematicVectorArray(len, job->JOB_TYPE);
            job->bodies[0]->points->b = KinematicVectorArray(len, job->JOB_TYPE);
            job->bodies[0]->points->m.resize(len);
            job->bodies[0]->points->v.resize(len);
            job->bodies[0]->points->v0.resize(len);
            job->bodies[0]->points->active.resize(len);
            job->bodies[0]->points->T = MaterialTensorArray(len);
            job->bodies[0]->points->L = KinematicTensorArray(len, job->JOB_TYPE);

            //zero out all entries to start
            job->bodies[0]->points->x.setZero();
            job->bodies[0]->points->u.setZero();
            job->bodies[0]->points->x_t.setZero();
            job->bodies[0]->points->m.setZero();
            job->bodies[0]->points->v.setZero();
            job->bodies[0]->points->v0.setZero();
            job->bodies[0]->points->mx_t.setZero();
            job->bodies[0]->points->b.setZero();
            job->bodies[0]->points->T.setZero();
            job->bodies[0]->points->L.setZero();
            job->bodies[0]->points->active.setOnes();

            //initial point positions, volumes/masses, pressure
            for (int i=0; i<N_p; i++){
                for (int j=0; j<N_p; j++){
                    job->bodies[0]->points->x(i*N_p + j,0) = (i + 0.5)/N_p;
                    job->bodies[0]->points->x(i*N_p + j,1) = (j + 0.5)/N_p;
                    job->bodies[0]->points->v(i*N_p + j) = 1.0/(N_p*N_p);
                    job->bodies[0]->points->m(i*N_p + j) = rho_s/(N_p*N_p) * (1.0 - get_porosity(job->bodies[0]->points->x(i*N_p + j)));
                }
            }

            job->bodies[0]->nodes = std::unique_ptr<Nodes>(new DefaultNodes);

            job->bodies[0]->boundary = std::unique_ptr<Boundary>(new CartesianBox);
            job->bodies[0]->activeBoundary = 1;

            job->bodies[0]->material = std::unique_ptr<Material>(new CompressibleNeohookeanElasticity);
            job->bodies[0]->material->fp64_props = {E, nu};
            job->bodies[0]->activeMaterial = 1;

            //assign grid
            fluid_grid = std::unique_ptr<FiniteVolumeGrid>(new FVMCartesian);
            fluid_grid->fp64_props = {1.0, 1.0, 0, 0, 0, 0}; //Lx
            fluid_grid->int_props = {N_e, N_e, 6, 6, 6, 6, N_q};  //Nx
            fluid_grid->str_props = {"USE_ENHANCED_QUADRATURE","USE_LOCAL_POROSITY_CORRECTION"};

            //initialize MPM objects
            std::cout << "Job properties (JOB_TYPE = " << job->JOB_TYPE << ", t = " << job->t << ", dt = " << job->dt << ")." << std::endl;
            job->thread_count = num_threads;
            if (job->thread_count > 1){
                std::cout << "Using " << job->thread_count << " threads." << std::endl;
            }
            std::cout << "Job Initialized." << std::endl;

            //do not start threadpool if thread_count == 1
            if (job->thread_count > 1) {
                job->threadPool = ThreadPool(job->thread_count);
            }

            job->serializer->init(job);
            job->grid->init(job);
            job->bodies[0]->init(job);

            //initialize FVM objects
            serializer->init(job, this);
            fluid_grid->init(job, this);
            solver->init(job, this);
            fluid_material->init(job, this);
            fluid_body->init(job, this);

            //first check number of elements
            std::cout << "Number of Elements: " << fluid_grid->element_count << std::endl;

            //check number of nodes
            std::cout << "Number of Nodes: " << job->grid->node_count << std::endl;

            //check number of points
            std::cout << "Number of Points: " << job->bodies[0]->points->x.size() << std::endl;

            //set counters to zero
            int stepCount = 0;
            int frameCount = 0;

            struct timespec timeStart, timeFrame, timeFinish;
            clock_gettime(CLOCK_MONOTONIC, &timeStart);
            timeFrame = timeStart;
            double tSim = 0;
            double tFrame = 0;
            double eSolid;
            double eFluid;
            double etmp;
            bool converged = false;
            bool exploded = false;
            std::vector<double> result = std::vector<double>();

            std::cout << "\n" << std::flush;
            //run simulation until stop_time
            while (job->t <= stop_time && !converged){

                //ensure simulation hasn't nan'd out
                for (int e=0; e<fluid_grid->element_count; e++){
                    if (!std::isfinite(fluid_body->p(e).norm())){
                        exploded = true;
                        break;
                    }
                }

                //run solver
                if (!exploded) {
                    solver->step(job, this);
                }

                //add correction to fluid density and energy
                for (int e=0; e<fluid_grid->element_count; e++) {
                    fluid_body->rho(e) += job->dt * get_fluid_mass_adjustment(job, fluid_grid->getElementCentroid(job, e), job->t);
                    fluid_body->rhoE(e) += job->dt * get_fluid_energy_adjustment(job, fluid_grid->getElementCentroid(job, e), job->t);
                }

                //increment time
                job->t += job->dt;

                if (job->serializer->writeFrame(job) == 1) {
                    //call fvm serializer to write frame as well:
                    serializer->writeFrame(job,this);

                    //successful frame written
                    clock_gettime(CLOCK_MONOTONIC,&timeFinish);
                    tFrame = (timeFinish.tv_sec - timeFrame.tv_sec) + (timeFinish.tv_nsec - timeFrame.tv_nsec)/1000000000.0;
                    tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
                    timeFrame = timeFinish;

                    //calculate solid and fluid phase error
                    eSolid = 0;
                    for (int i=0; i<job->grid->node_count; i++){
                        if (job->bodies[0]->nodes->m(i) > 1e-10){
                            etmp = (job->bodies[0]->nodes->x_t(i) - get_solid_velocity(job, job->bodies[0]->nodes->x(i), job->t)).norm();
                            eSolid += etmp*etmp * job->grid->nodeVolume(job,i);
                        }
                    }
                    eFluid = 0;
                    for (int e=0; e<fluid_grid->element_count; e++){
                        etmp = (fluid_body->p(e)/fluid_body->rho(e) - get_fluid_velocity(job, fluid_grid->getElementCentroid(job,e), job->t)).norm();
                        eFluid += etmp * etmp * fluid_grid->getElementVolume(e);
                    }

                    //sqrt for L2 norm
                    eSolid = std::sqrt(eSolid);
                    eFluid = std::sqrt(eFluid);

                    //add chain of errors to output
                    result.push_back(eFluid);
                    result.push_back(eSolid);

                    std::cout << "\33[2K" << "\x1b[A" << "\33[2K" << "\r";
                    //printf("\33[2K");
                    std::cout << "Frame Written [" << ++frameCount << "]. Time/Frame [" << tFrame << " s]. Elapsed Time [" << tSim << " s]." << std::endl;
                    std::cout << "L2 Velocity Error : Solid [" << eSolid << "], Fluid [" << eFluid << "]" << std::flush;
                }
            }
            clock_gettime(CLOCK_MONOTONIC,&timeFinish);
            tSim = (timeFinish.tv_sec - timeStart.tv_sec) + (timeFinish.tv_nsec - timeStart.tv_nsec)/1000000000.0;
            std::cout << std::endl << std::endl << "Simulation Complete. Elapsed Time [" << tSim << "s]." << std::endl;

            return result;
        }
    };
}

void fvm_test(Job* job){
    std::unique_ptr<FiniteVolumeDriver> driver = std::unique_ptr<FiniteVolumeDriver>(new FVM_TEST::FVMTestDriver);
    driver->init(job);
    driver->run(job);
    return;
}

void fvm_mpm_drag_test(Job* job){
    FVM_TEST::FVMMPMDragTestDriver driver = FVM_TEST::FVMMPMDragTestDriver();
    driver.init(job);
    //driver->run(job);
    //driver.get_L1_errors(job, 100001, 100000);

    int n_i_max = 10241;
    int n_e_max = 10240;

    //set up study
    std::vector<int> n_i = {11,21,41,81,161,321,641,1281,2561,5121,10241};
    std::vector<int> n_e = {10,20,40,80,160,320,640,1280,2560,5120,10240};

    /*
    //set up data containers
    std::vector<double> f_ii_L1 = std::vector<double>(n_i.size()); //L1 error for f_i vs. n_i
    std::vector<double> f_ei_L1 = std::vector<double>(n_i.size()); //L1 error for f_e vs. n_i
    std::vector<double> f_ie_L1 = std::vector<double>(n_e.size()); //L1 error for f_i vs. n_e
    std::vector<double> f_ee_L1 = std::vector<double>(n_e.size()); //L1 error for f_e vs. n_e

    //run study
    std::vector<double> result;
    for (int i=0; i<n_i.size(); i++){
        result = driver.get_L1_errors(job, n_i[i], n_e_max, 102400, "Linear1DNonUniform");
        f_ii_L1[i] = result[0];
        f_ei_L1[i] = result[1];
    }

    for (int i=0; i<n_e.size(); i++){
        result = driver.get_L1_errors(job, n_i_max, n_e[i], 102400, "Linear1DNonUniform");
        f_ie_L1[i] = result[0];
        f_ee_L1[i] = result[1];
    }

    //print results
    std::cout << std::endl;
    std::cout << "N_i Convergence Study" << std::endl;
    for (int i=0; i<n_i.size(); i++){
        std::cout << n_i[i] << ", " << n_e_max << ": " << f_ii_L1[i] << " : " << f_ei_L1[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "N_e Convergence Study" << std::endl;
    for (int i=0; i<n_e.size(); i++){
        std::cout << n_i_max << ", " << n_e[i] << ": " << f_ie_L1[i] << " : " << f_ee_L1[i] << std::endl;
    }
     */

    //second study
    //set up data containers
    std::vector<int> n_i_result = std::vector<int>();
    std::vector<int> n_e_result = std::vector<int>();
    std::vector<double> f_i_L1 = std::vector<double>(); //L1 error for f_i
    std::vector<double> f_e_L1 = std::vector<double>(); //L1 error for f_e

    //run study
    std::vector<double> result;
    for (int i=0; i<n_i.size(); i++){
        for (int j=0; j<n_e.size(); j++) {
            result = driver.get_L1_errors(job, n_i[i], n_e[j], 102400, "Linear1DNonUniform");
            n_i_result.push_back(n_i[i]);
            n_e_result.push_back(n_e[j]);
            f_i_L1.push_back(result[0]);
            f_e_L1.push_back(result[1]);
        }
    }

    //print results
    std::cout << std::endl;
    std::cout << "Convergence Study" << std::endl;
    for (int i=0; i<n_i_result.size(); i++){
        std::cout << n_i_result[i] << ", " << n_e_result[i] << ": " << f_i_L1[i] << " : " << f_e_L1[i] << std::endl;
    }

    return;
}

void fvm_mpm_buoyancy_test(Job* job){
    FVM_TEST::FVMMPMBuoyancyTestDriver driver = FVM_TEST::FVMMPMBuoyancyTestDriver();
    driver.init(job);
    //driver->run(job);
    //driver.get_L1_errors(job, 100001, 100000);

    int n_i_max = 10241;
    int n_e_max = 10240;

    //set up study
    std::vector<int> n_i = {11,21,41,81,161,321,641,1281,2561,5121,10241};
    std::vector<int> n_e = {10,20,40,80,160,320,640,1280,2560,5120,10240};

    //set up data containers
    std::vector<double> f_ie_L1 = std::vector<double>(n_e.size()); //L1 error for f_i vs. n_e
    std::vector<double> f_ee_L1 = std::vector<double>(n_e.size()); //L1 error for f_e vs. n_e

    //run study
    std::vector<double> result;
    for (int i=0; i<n_e.size(); i++){
        result = driver.get_L1_errors(job, n_i[i], n_e[i], 102400);//, "Linear1DNonUniform");
        f_ie_L1[i] = result[0];
        f_ee_L1[i] = result[1];
    }

    //print results
    std::cout << std::endl;
    std::cout << "N_e Convergence Study" << std::endl;
    for (int i=0; i<n_e.size(); i++){
        std::cout << n_i[i] << ", " << n_e[i] << ": " << f_ie_L1[i] << " : " << f_ee_L1[i] << std::endl;
    }

    return;
}

void fvm_mpm_porosity_test(Job* job){
    FVM_TEST::FVMMPMPorosityTestDriver driver = FVM_TEST::FVMMPMPorosityTestDriver();
    driver.init(job);
    //driver->run(job);
    //driver.get_L1_errors(job, 100001, 100000);

    int n_i_max = 10241;
    int n_p_max = 10240;

    //set up study
    std::vector<int> n_i = {11,21,41,81,161,321,641,1281,2561,5121,10241};
    std::vector<int> n_p = {10,20,40,80,160,320,640,1280,2560,5120,10240};

    //set up data containers
    std::vector<double> f_i_L1_cst_hp1 = std::vector<double>(n_i.size()); //L1 error (10240 p)
    std::vector<double> f_i_L1_cst_hp2 = std::vector<double>(n_i.size()); //L1 error (5120 p)
    std::vector<double> f_i_L1_cst_ppc1 = std::vector<double>(n_i.size()-4); //L1 error (16 ppc)
    std::vector<double> f_i_L1_cst_ppc2 = std::vector<double>(n_i.size()-4); //L1 error (8 ppc)
    std::vector<double> f_i_L1_cst_hi1 = std::vector<double>(n_p.size()-6); //L1 error (641 n)
    std::vector<double> f_i_L1_cst_hi2 = std::vector<double>(n_p.size()-6); //L1 error (321 n)


    //run studies
    for (int i=0; i<n_i.size(); i++){
        f_i_L1_cst_hp1[i] = driver.get_L1_error(job, n_i[i], 10240, 102400, "Linear1DNonUniform");
        f_i_L1_cst_hp2[i] = driver.get_L1_error(job, n_i[i], 5120, 102400, "Linear1DNonUniform");
    }
    for (int i=0; i<n_i.size()-4; i++){
        f_i_L1_cst_ppc1[i] = driver.get_L1_error(job, n_i[i], n_p[i+4], 102400, "Linear1DNonUniform");
        f_i_L1_cst_ppc2[i] = driver.get_L1_error(job, n_i[i], n_p[i+3], 102400, "Linear1DNonUniform");
    }
    for (int i=0; i<n_p.size()-6; i++){
        f_i_L1_cst_hi1[i] = driver.get_L1_error(job, 641, n_p[i+6], 102400, "Linear1DNonUniform");
        f_i_L1_cst_hi2[i] = driver.get_L1_error(job, 321, n_p[i+5], 102400, "Linear1DNonUniform");
    }

    //print results
    std::cout << std::endl;
    std::cout << "N_i Convergence Study" << std::endl;
    for (int i=0; i<n_i.size(); i++){
        std::cout << n_i[i] << ", " << 10240 << ": " << f_i_L1_cst_hp1[i] << std::endl;
    }

    std::cout << std::endl;
    std::cout << "N_i Convergence Study" << std::endl;
    for (int i=0; i<n_i.size(); i++){
        std::cout << n_i[i] << ", " << 5120 << ": " << f_i_L1_cst_hp2[i] << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Lmpp Convergence Study" << std::endl;
    for (int i=0; i<n_i.size()-4; i++){
        std::cout << n_i[i] << ", " << n_p[i+4] << ": " << f_i_L1_cst_ppc1[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Lmpp Convergence Study" << std::endl;
    for (int i=0; i<n_i.size()-4; i++){
        std::cout << n_i[i] << ", " << n_p[i+3] << ": " << f_i_L1_cst_ppc2[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "N_p Convergence Study" << std::endl;
    for (int i=0; i<n_p.size()-6; i++){
        std::cout << 641 << ", " << n_p[i+6] << ": " << f_i_L1_cst_hi1[i] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "N_p Convergence Study" << std::endl;
    for (int i=0; i<n_p.size()-6; i++){
        std::cout << 321 << ", " << n_p[i+5] << ": " << f_i_L1_cst_hi1[i] << std::endl;
    }


    return;
}

void fvm_mpm_moms_test(Job* job) {
    FVM_TEST::FVMMPMMoMSDriver driver = FVM_TEST::FVMMPMMoMSDriver();
    driver.init(job);
    //driver->run(job);
    //driver.get_L1_error(job, 41, 40, 40, 1, 1e-4, "414040v2", 4);

    //return;

    /*
    double dt = 1e-5;
    std::vector<double> results = driver.get_L1_error(job, 41, 160, 40, 1, dt, "4016040v2", 10);

    //save values to file
    std::ofstream file2 = std::ofstream();
    file2.open("output/data_v2.txt");

    //print results to console
    file2 << 1e-5 << ", " << 41 << ", " << 160 << ", " << 40;
    for (int ii=0; ii<results.size(); ii++){
        file2 << ", " << results[ii];
    }
    file2 << std::endl;

    file2.close();

    return;
     */

    //run tests for h_i convergence
    double dt = 5e-4; //0.05 / 40.0 / 40.0;
    std::vector<double> Dt = {dt,  dt,   dt,   dt,   dt,   dt,  dt,  dt,  dt,  dt, dt,  dt,  2.0*dt, dt/2.0, dt/4.0, dt/8.0};
    std::vector<int> N_i =   {5+1, 10+1, 20+1, 30+1, 60+1, 61,  61,  61,  61,  61, 61,  61,  61,     61,     61,     61};
    std::vector<int> N_p =   {20,  40,   80,   120,  240,  240, 240, 240, 240, 60, 120, 180, 240,    240,    240,    240};
    std::vector<int> N_e =   {60,  60,   60,   60,   60,   5,   10,  20,  30,  60, 60,  60,  60,     60,     60,     60};
    std::vector<int> N_q =   { 1,   1,    1,    1,    2,   12,   6,   3,   2,   2,  2,   2,   2,      2,      2,      2};

    std::vector<double> results = {};
    std::vector<double> tmp;
    int len = 0;
    std::string sim_name;
    for (int i=0; i<N_i.size(); i++){
        sim_name = std::to_string(N_i[i]) + std::to_string(N_p[i]) + std::to_string(N_e[i]) + std::to_string(Dt[i]);
        tmp = driver.get_L1_error(job, N_i[i], N_p[i], N_e[i], N_q[i], Dt[i], sim_name, 20);
        len = tmp.size();
        for (int ii=0; ii<len; ii++){
            results.push_back(tmp[ii]);
        }
    }

    //print results to console
    std::cout << "Convergence Study" << std::endl;
    std::cout << "dt, N_i, N_p, N_e, E_s(0), E_f(0), ..., E_s(1.0), E_f(1.0)" << std::endl;
    for (int i=0; i<N_i.size(); i++){
        std::cout << Dt[i] << ", " << N_i[i] << ", " << N_p[i] << ", " << N_e[i];
        for (int ii=0; ii<len; ii++){
            std::cout << ", " << results[i*len + ii];
        }
        std::cout << std::endl;
    }

    //save values to file
    std::ofstream file = std::ofstream();
    file.open("output/data.txt");

    file << "Convergence Study" << std::endl;
    file << "dt, N_i, N_p, N_e, E_s(0), E_f(0), ..., E_s(1.0), E_f(1.0)" << std::endl;
    for (int i=0; i<N_i.size(); i++){
        file << Dt[i] << ", " << N_i[i] << ", " << N_p[i] << ", " << N_e[i];
        for (int ii=0; ii<len; ii++){
            file << ", " << results[i*len + ii];
        }
        file << std::endl;
    }

    file.close();


    return;
}

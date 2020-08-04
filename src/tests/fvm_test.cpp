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
//
// Created by aaron on 5/10/18.
// drivers.hpp
//

#ifndef MPM_V3_DRIVERS_HPP
#define MPM_V3_DRIVERS_HPP

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "job.hpp"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <time.h>
#include <random>

/*
 * IN THIS FILE, DEFINE DRIVER OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
class Driver : public MPMObject{
public:
    //functions which must be implemented by every driver
    virtual void init(Job*) = 0;                                        //initialize from Job
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save to file (in given directory) and return filename
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load from file
    virtual void run(Job*) = 0;                                         //run mpm according to problem
    virtual void generateGravity(Job*) = 0;                             //generate gravity
    virtual void applyGravity(Job*) = 0;                                //apply gravity
};
 */

class DefaultDriver : public Driver{
public:
    DefaultDriver(){
        object_name = "DefaultDriver"; //set object name here
    }

    double stop_time;
    KinematicVector gravity;

    virtual void init(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    virtual void run(Job* job);
    virtual void generateGravity(Job* job);
    virtual void applyGravity(Job* job);
};

/*----------------------------------------------------------------------------*/

class ColumnCollapseDriver : public Driver{
public:
    ColumnCollapseDriver(){
        object_name = "ColumnCollapseDriver";
    }

    double stop_time;
    KinematicVector gravity;

    void setPressure(Job* job);

    void init(Job* job);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void run(Job* job);
    void generateGravity(Job* job);
    void applyGravity(Job* job);
};

/*----------------------------------------------------------------------------*/

class UserDefinedGravityDriver : public Driver{
public:
    UserDefinedGravityDriver(){
        object_name = "UserDefinedGravityDriver";
    }

    double stop_time;
    KinematicVector gravity;

    void init(Job* job);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void run(Job* job);
    void generateGravity(Job* job);
    void applyGravity(Job* job);
};

/*----------------------------------------------------------------------------*/

class CavityFlowDriver : public Driver{
public:
    CavityFlowDriver(){
        object_name = "CavityFlowDriver";
    }

    double stop_time;
    double Ly, hy;
    double v_set;

    void init(Job* job);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void run(Job* job);
    void generateGravity(Job* job);
    void applyGravity(Job* job);
};


/*----------------------------------------------------------------------------*/
//driver for ballistic impact problems
//will launch ballistic object at specified velocity at t=0

class BallisticDriver : public DefaultDriver{
public:
    BallisticDriver(){
        object_name = "BallisticDriver"; //set object name here
    }

    int ballistic_id;
    double stop_time;
    KinematicVector gravity, velocity;

    void init(Job* job);

    void run(Job* job);
};

/*----------------------------------------------------------------------------*/
//Fish Random Sampler (6.883 project, spring 2020)
class FishRandomSampleDriver : public Driver{
public:
    FishRandomSampleDriver(){
        object_name = "FishRandomSampleDriver"; //set object name here
    }

    int num_samples;
    std::string output_file;
    double alpha_min, alpha_max;
    double beta_min, beta_max;
    double gamma_min, gamma_max;

    int fish_body;

    double stop_time;
    KinematicVector gravity;

    virtual void init(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    virtual void run(Job* job);
    virtual void generateGravity(Job* job);
    virtual void applyGravity(Job* job);
};

class FishMixtureModelRandomSampleDriver : public FishRandomSampleDriver{
public:
    FishMixtureModelRandomSampleDriver(){
        object_name = "FishMixtureModelRandomSampleDriver"; //set object name here
    }

    //from https://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
    struct MultivariateNormalRandomVariable
    {
        Eigen::VectorXd mean;
        Eigen::MatrixXd transform;

        void setMean(Eigen::VectorXd const& meanIn){
            mean = meanIn;
        }

        void setCovar(Eigen::MatrixXd const& covarIn){
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covarIn);
            transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
        }

        Eigen::VectorXd operator()() const
        {
            static std::mt19937 gen{ std::random_device{}() };
            static std::normal_distribution<> dist;

            return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](double x) { return dist(gen); });
        }
    };

    std::string input_file;
    std::vector<MultivariateNormalRandomVariable> normal_distributions;
    std::vector<double> probability_distribution;

    virtual void init(Job* job);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    virtual void run(Job* job);
    virtual void generateGravity(Job* job);
    virtual void applyGravity(Job* job);
};

/*----------------------------------------------------------------------------*/
//driver for chute flow with Sachith's sand model
class ChuteFlowDriver : public DefaultDriver{
public:
    ChuteFlowDriver(){
        object_name = "ChuteFlowDriver"; //set object name here
    }

    double stop_time;
    KinematicVector gravity;

    double g, theta, mu_1, mu_2, b, d, rho_s, h, phi;
    Eigen::VectorXd V_0;

    //need to create and write an output file
    std::string output_filename;

    void init(Job* job);

    void run(Job* job);
    void generateGravity(Job* job);
    void applyGravity(Job* job);

    KinematicVector getVelocity(Job* job, KinematicVector const &x);
    void writeErrorInfo(Job* job);
};

#endif //MPM_V3_DRIVERS_HPP

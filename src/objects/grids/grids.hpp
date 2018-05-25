//
// Created by aaron on 5/15/18.
// grids.hpp
//

#ifndef MPM_V3_GRIDS_HPP
#define MPM_V3_GRIDS_HPP

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

/*
 * IN THIS FILE, DEFINE GRID OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
 * class Grid : public MPMObject{
public:
    size_t node_count;      //number of nodes which define grid
    size_t element_count;   //number of elements in grid

    virtual void init(Job*) = 0; //initialize from Job

    virtual void writeFrame(Job*, Serializer*) = 0;                     //send frame data to Serializer
    virtual std::string saveState(Job*, Serializer*, std::string) = 0;  //save to file
    virtual int loadState(Job*, Serializer*, std::string) = 0;          //load from file

    virtual void writeHeader(Job*, Body*, Serializer*, std::ofstream&, int) = 0; //write cell types
    virtual int whichElement(Job*, KinematicVector&) = 0;                                                       //return element given position
    virtual bool inDomain(Job*, KinematicVector&) = 0;                                                          //check if location is in domain
    virtual KinematicVector nodeIDToPosition(Job*, int) = 0;                                                    //return position of given node
    virtual void evaluateBasisFnValue(Job*, KinematicVector&, std::vector<int>&, std::vector<double>&) = 0;     //fill id and value for nodal weights of position
    virtual void evaluateBasisFnGradient(Job*, KinematicVector&, std::vector<int>&, KinematicVectorArray&) = 0; //fill id and gradient for nodal gradient of position
    virtual double nodeVolume(Job*, int) = 0;                                                                   //return node volume of id
    virtual double elementVolume(Job*, int) = 0;                                                                //return element volume of id
    virtual int nodeTag(Job*, int) = 0;                                                                         //return 'tag' of node id
};
 */

class CartesianLinear : public Grid{
public:
    CartesianLinear(){
        object_name = "CartesianLinear";
    }

    KinematicVector Lx, hx;
    Eigen::VectorXi Nx;
    KinematicVectorArray x_n;
    Eigen::MatrixXi nodeIDs; //element to node map
    Eigen::MatrixXi A; //0,1 for directions
    size_t npe; //nodes per element
    Eigen::VectorXd v_n; //nodal volume
    double v_e; //element volume

    static int cartesianWhichElement(Job* job, KinematicVector& xIN, KinematicVector& LxIN, KinematicVector& hxIN, Eigen::VectorXi& NxIN);

    virtual void init(Job* job);
    virtual void hiddenInit(Job* job);

    virtual void writeFrame(Job* job, Serializer* serializer);
    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    virtual void writeHeader(Job* job, Body* body, Serializer* serializer, std::ofstream& nfile, int SPEC); //write cell types

    virtual int whichElement(Job* job, KinematicVector& xIN);
    virtual bool inDomain(Job* job, KinematicVector& xIN);
    virtual KinematicVector nodeIDToPosition(Job* job, int idIN);

    virtual void evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL);
    virtual void evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD);
    virtual double nodeVolume(Job* job, int idIN);
    virtual double elementVolume(Job* job, int idIN);
    virtual int nodeTag(Job* job, int idIN);
};

/*----------------------------------------------------------------------------*/

class CartesianPeriodic : public CartesianLinear{
public:
    CartesianPeriodic(){
        object_name = "CartesianPeriodic";
    }

    Eigen::VectorXi nntoni; //node number to node id

    void hiddenInit(Job* job);

    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void fixPosition(Job* job, KinematicVector& xIN);
    int whichElement(Job* job, KinematicVector& xIN);
    bool inDomain(Job* job, KinematicVector& xIN);

    void evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL);
    void evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD);
};

/*----------------------------------------------------------------------------*/

class CartesianCubic : public Grid{
public:
    CartesianCubic(){
        object_name = "CartesianCubic";
    }

    KinematicVector Lx, hx;
    Eigen::VectorXi Nx;
    KinematicVectorArray x_n, edge_n;
    Eigen::MatrixXi nodeIDs; //element to node map
    Eigen::MatrixXi A; //0,1 for directions
    size_t npe; //nodes per element
    Eigen::VectorXd v_n; //nodal volume
    double v_e; //element volume

    double s(double, double);
    double g(double, double);

    void init(Job* job);
    void hiddenInit(Job* job);

    void writeFrame(Job* job, Serializer* serializer);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    void writeHeader(Job* job, Body* body, Serializer* serializer, std::ofstream& nfile, int SPEC); //write cell types

    int whichElement(Job* job, KinematicVector& xIN);
    bool inDomain(Job* job, KinematicVector& xIN);
    KinematicVector nodeIDToPosition(Job* job, int idIN);

    void evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL);
    void evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD);
    double nodeVolume(Job* job, int idIN);
    double elementVolume(Job* job, int idIN);
    int nodeTag(Job* job, int idIN);
};

#endif //MPM_V3_GRIDS_HPP

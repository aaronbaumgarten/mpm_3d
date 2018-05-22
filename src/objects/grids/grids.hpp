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

    void init(Job* job);
    void hiddenInit(Job* job);

    void writeFrame(Job* job, Serializer* serializer);
    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    int whichElement(Job* job, KinematicVector& xIN);
    static int cartesianWhichElement(Job* job, KinematicVector& xIN, KinematicVector& LxIN, KinematicVector& hxIN, Eigen::VectorXi& NxIN);
    bool inDomain(Job* job, KinematicVector& xIN);
    KinematicVector nodeIDToPosition(Job* job, int idIN);

    void evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL);
    void evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD);
    double nodeVolume(Job* job, int idIN);
    double elementVolume(Job* job, int idIN);
    int nodeTag(Job* job, int idIN);
};

#endif //MPM_V3_GRIDS_HPP

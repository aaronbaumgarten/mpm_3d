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
    int node_count;      //number of nodes which define grid
    int element_count;   //number of elements in grid

    int GRID_DIM;

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
    virtual double nodeSurfaceArea(Job*, int) = 0;                                                          //return node surface volume
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
    int npe; //nodes per element
    Eigen::VectorXd v_n; //nodal volume
    double v_e; //element volume
    Eigen::VectorXd s_n; //surface integral

    static int cartesianWhichElement(Job* job, KinematicVector& xIN, KinematicVector& LxIN, KinematicVector& hxIN, Eigen::VectorXi& NxIN, int GRID_DIM_IN = -1);

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
    virtual double nodeSurfaceArea(Job* job, int idIN);
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

class CartesianCustom : public CartesianLinear{
public:
    CartesianCustom(){
        object_name = "CartesianCustom";
    }

    static const int WALL = 0; //for periodic props definition
    static const int PERIODIC = 1; //for periodic definition

    Eigen::VectorXi nntoni; //node number to node id
    Eigen::VectorXi periodic_props; //boundary props

    void init(Job* job);
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
    int npe; //nodes per element
    Eigen::VectorXd v_n; //nodal volume
    double v_e; //element volume
    Eigen::VectorXd s_n; //surface integral

    MPMScalarSparseMatrix S_grid = MPMScalarSparseMatrix(0,0); //node value to function value map

    static double s(double, double);
    static double g(double, double);

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
    virtual double nodeSurfaceArea(Job* job, int idIN);
};

/*----------------------------------------------------------------------------*/

class CartesianCubicCustom : public CartesianCubic{
public:
    CartesianCubicCustom(){
        object_name = "CartesianCubicCustom";
    }

    static const int WALL = 0; //for periodic props definition
    static const int PERIODIC = 1; //for periodic definition

    Eigen::VectorXi periodic_props; //boundary props
    Eigen::VectorXi nntoni; //node number to node id

    virtual void init(Job* job);
    virtual void hiddenInit(Job* job);

    virtual std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Serializer* serializer, std::string fullpath);

    virtual void fixPosition(Job* job, KinematicVector& xIN);
    virtual int whichElement(Job* job, KinematicVector& xIN);
    virtual bool inDomain(Job* job, KinematicVector& xIN);

    virtual void evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL);
    virtual void evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD);
};

/*----------------------------------------------------------------------------*/

class CartesianCubic_Offset : public CartesianCubic{
public:
    CartesianCubic_Offset(){
        object_name = "CartesianCubic_Offset";
    }

    double x_offset; //this grid really should only be used for axisymmetric simulations

    void init(Job* job);
    void hiddenInit(Job* job);

    std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    int loadState(Job* job, Serializer* serializer, std::string fullpath);

    int whichElement(Job* job, KinematicVector& xIN);
    bool inDomain(Job* job, KinematicVector& xIN);

    KinematicVector nodeIDToPosition(Job* job, int idIN);
};

/*----------------------------------------------------------------------------*/

class TriangularGridLinear : public Grid{
public:
    TriangularGridLinear(){
        object_name = "TriangularGridLinear";
    }

    double lc, msh_version;
    std::string msh_filename;
    KinematicVector Lx, hx, x_min;
    Eigen::VectorXi Nx;
    Eigen::MatrixXi nodeIDs; //element to node map
    Eigen::VectorXi nodeTags; //node tags

    int npe = 3;
    KinematicVectorArray x_n;
    KinematicTensorArray A, Ainv; //map xi to x and x to xi

    Eigen::VectorXd v_n, v_e; //nodal/element volumes
    Eigen::VectorXd s_n; //surface integral

    std::vector<int> search_cells; //search grid (cell to element map)
    std::vector<int> search_offsets; //search offsets

    virtual int ijk_to_n(Job* job, const Eigen::VectorXi& ijk);
    virtual Eigen::VectorXi n_to_ijk(Job* job, int n);

    virtual bool line_segment_intersect(Eigen::VectorXd s0_p0, Eigen::VectorXd s0_p1, Eigen::VectorXd s1_p0, Eigen::VectorXd s1_p1);
    virtual int whichSearchCell(const KinematicVector& xIN);
    virtual bool inElement(Job* job, const KinematicVector& xIN, int idIN);

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
    virtual double nodeSurfaceArea(Job* job, int idIN);
};

/*----------------------------------------------------------------------------*/

class TetrahedralGridLinear : public TriangularGridLinear{
public:
    TetrahedralGridLinear(){
        object_name = "TetrahedralGridLinear";
        npe = 4;
    }

    //cutoff for number of search cells to switch to large domain version of code.
    const static int LARGE_DOMAIN_CUTOFF = 200000;

    //search grid for very large matrices
    bool large_domain = false;
    std::vector<std::vector<int>> search_cells_large_domain;

    //search method for very large, irregular grids
    bool use_n_to_e = false;
    std::vector<std::vector<int>> node_to_element_map;
    std::vector<std::vector<int>> element_neighbors;

    //stride is 6 for 3D
    std::vector<int> element_min_max; //list of minimum and maximum ijk of nodes

    virtual int whichSearchCell(const KinematicVector& xIN);
    virtual bool inElement(Job* job, const KinematicVector& xIN, int idIN);

    virtual void init(Job* job);
    virtual void hiddenInit(Job* job);

    //void writeFrame(Job* job, Serializer* serializer);
    //std::string saveState(Job* job, Serializer* serializer, std::string filepath);
    //int loadState(Job* job, Serializer* serializer, std::string fullpath);

    virtual void writeHeader(Job* job, Body* body, Serializer* serializer, std::ofstream& nfile, int SPEC); //write cell types

    virtual int whichElement(Job* job, KinematicVector& xIN);
    virtual bool inDomain(Job* job, KinematicVector& xIN);
    virtual KinematicVector nodeIDToPosition(Job* job, int idIN);

    virtual void evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL);
    virtual void evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD);
    virtual double nodeVolume(Job* job, int idIN);
    virtual double elementVolume(Job* job, int idIN);
    virtual int nodeTag(Job* job, int idIN);
    virtual double nodeSurfaceArea(Job* job, int idIN);

    virtual int whichElement(Job* job, KinematicVector& xIN, int& elem_guess); //return element given position
    virtual void evaluateBasisFnValue(Job* job, KinematicVector& xIN, std::vector<int>& nID, std::vector<double>& nVAL, int& elem_guess);     //fill id and value for nodal weights of position
    virtual void evaluateBasisFnGradient(Job* job, KinematicVector& xIN, std::vector<int>& nID, KinematicVectorArray& nGRAD, int& elem_guess); //fill id and gradient for nodal gradient of position

};

class Linear1DNonUniform : public Grid{
public:
    Linear1DNonUniform(){
        object_name = "Linear1DNonUniform";
    }

    KinematicVector Lx, hx;
    Eigen::VectorXi Nx;
    KinematicVectorArray x_n;
    Eigen::MatrixXi nodeIDs; //element to node map
    Eigen::MatrixXi A; //0,1 for directions
    int npe; //nodes per element
    Eigen::VectorXd v_n; //nodal volume
    Eigen::VectorXd v_e; //element volume
    Eigen::VectorXd s_n; //surface integral

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
    virtual double nodeSurfaceArea(Job* job, int idIN);
};

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

#endif //MPM_V3_GRIDS_HPP

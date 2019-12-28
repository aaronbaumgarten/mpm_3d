//
// Created by aaron on 12/12/19.
// fvm_grids.hpp
//

#ifndef MPM_V3_FVM_GRIDS_HPP
#define MPM_V3_FVM_GRIDS_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <fstream>

#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"
#include "mpm_objects.hpp"
#include "fvm_objects.hpp"

/*
//finite volume grid class
//responsible for tracking element, face definitions and cross integrals with mpm grid
class FiniteVolumeGrid: public MPMObject{
public:
    //boundary condition types
    static const int DIRICHLET = 0;
    static const int NEUMANN = 1;
    static const int PERIODIC = 2;

    int face_count;      //number of faces which define grid
    int node_count;      //number of nodes which define grid
    int element_count;   //number of elements in grid

    int GRID_DIM = -1; //dimension of grid (might be different than simulation dimension)

    MPMScalarSparseMatrix M = MPMScalarSparseMatrix(0,0); //M_ji maps jth element to ith mpm node
    KinematicVectorSparseMatrix gradM = KinematicVectorSparseMatrix(0,0); //gradM_ji maps jth element to ith node gradient

    virtual void init(Job*, FiniteVolumeDriver*) = 0; //initialize from mpm Job

    virtual void writeHeader(std::ofstream&, int) = 0; //write cell types

    virtual double getElementVolume(int) = 0;      //return element volume of id
    virtual int getElementTag(int) = 0;            //return 'tag' of id
    virtual double getFaceArea(int) = 0;           //return face area of id
    virtual int getFaceTag(int) = 0;               //return 'tag of face id
    virtual KinematicVector getFaceNormal(Job*, int) = 0; //return normal vector associated with given id

    virtual std::vector<int> getElementFaces(int) = 0; //return faces associated with given element id
    virtual std::array<int,2> getOrientedElementsByFace(int) = 0; //return elements associated with given face id A -|-> B
    virtual std::vector<int> getElementNeighbors(int) = 0; //return list of elements in neighborhood
    virtual KinematicVector getElementCentroid(int) = 0; //return element centroid
    virtual KinematicVector getFaceCentroid(Job*, int) = 0;

    virtual void generateMappings(Job*, FiniteVolumeDriver*) = 0; //generate M_ji and gradM_ji maps
    virtual void constructMomentumField(Job*, FiniteVolumeDriver*) = 0; //construct momentum from FV body
    virtual void constructDensityField(Job*, FiniteVolumeDriver*) = 0; //construct momentum from FV body
    virtual KinematicTensorArray getVelocityGradients(Job*, FiniteVolumeDriver*) = 0; //return velocity gradient in each element

    //functions to compute element-wise fluxes of field variables using reconstructed velocity field
    virtual Eigen::VectorXd calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, Eigen::VectorXd&) = 0;
    //virtual KinematicVectorArray calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, KinematicVectorArray&, int) = 0;
    //virtual KinematicTensorArray calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, KinematicTensorArray&, int) = 0;
    //virtual MaterialVectorArray calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, MaterialVectorArray&, int) = 0;
    //virtual MaterialTensorArray calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, MaterialTensorArray&, int) = 0;

    //functions to compute element mass flux
    virtual Eigen::VectorXd calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver) = 0;

    //functions to compute element momentum fluxes
    virtual KinematicVectorArray calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver) = 0;
};
*/

class FVMCartesian : public FiniteVolumeGrid{
public:
    FVMCartesian(){
        object_name = "FVMCartesian";
    }

    //grid dimensions
    std::vector<int> Nx;
    std::vector<double> Lx, hx;
    double element_volumes;
    std::vector<double> face_areas;
    std::vector<int> face_normals; //x,y,z only
    int num_neighbors;

    //boundary conditions
    //-x,+x,-y,+y,-z,+z
    std::vector<int> bc_tags;
    std::vector<KinematicVector> bc_values;
    std::vector<int> face_bcs;                      //vector of boundary definition (-x -> 0, +y -> 3, etc)

    //grid definitions
    std::vector<std::array<int,2>> face_elements;   //oriented face list (A -|-> B)
    std::vector<std::vector<int>> element_faces;    //unordered list of faces
    std::vector<std::vector<int>> element_neighbors; //list of nearest neighbors

    //gradient reconstruction vectors for solving least squares problem
    std::vector<Eigen::MatrixXd> A_e;
    std::vector<Eigen::MatrixXd> A_inv; //psuedo inverse
    std::vector<Eigen::VectorXd> b_e;

    //mapping functions
    std::vector<int> e_to_ijk(int e); //convert element id to ijk element definition
    int ijk_to_e(std::vector<int> ijk); //convert ijk element definition to e
    std::vector<int> n_to_ijk(int n); //convert node id to ijk
    int ijk_to_n(std::vector<int> ijk); //convert ijk to node id

    virtual void init(Job* job, FiniteVolumeDriver* driver);

    virtual void writeHeader(std::ofstream& file, int TYPE);

    virtual double getElementVolume(int e);
    virtual int getElementTag(int e);
    virtual double getFaceArea(int f);
    virtual int getFaceTag(int f);
    virtual KinematicVector getFaceNormal(Job* job, int f);

    virtual std::vector<int> getElementFaces(int e);
    virtual std::array<int,2> getOrientedElementsByFace(int e);
    virtual std::vector<int> getElementNeighbors(int e);
    virtual KinematicVector getElementCentroid(Job* job, int e); //return element centroid
    virtual KinematicVector getFaceCentroid(Job* job, int f);

    virtual void generateMappings(Job* job, FiniteVolumeDriver* driver);
    virtual void constructMomentumField(Job* job, FiniteVolumeDriver* driver);
    virtual void constructDensityField(Job* job, FiniteVolumeDriver* driver);
    virtual KinematicTensorArray getVelocityGradients(Job* job, FiniteVolumeDriver* driver);

    //functions to compute element-wise fluxes of field variables using reconstructed velocity field
    virtual Eigen::VectorXd calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, Eigen::VectorXd& values);

    //functions to compute element mass flux
    virtual Eigen::VectorXd calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver);

    //functions to compute element momentum fluxes
    virtual KinematicVectorArray calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver);
};

/*----------------------------------------------------------------------------*/

class FVMGmsh2D : public FiniteVolumeGrid{
public:
    FVMGmsh2D(){
        object_name = "FVMGmsh2D";
    }

    //grid definitions
    std::string filename;
    Eigen::VectorXd face_areas;
    KinematicVectorArray face_normals; //x,y,z only

    Eigen::MatrixXi nodeIDs;        //element to node map
    int npe = 3;                    //nodes per element
    KinematicVectorArray x_n, x_e, x_f;  //node positions, element centroids, face_centroid
    Eigen::VectorXd v_e;            //element volumes

    //boundary conditions
    Eigen::VectorXi bc_tags;        //tags for each face
    KinematicVectorArray bc_values; //bc values at each face

    //grid definitions
    std::vector<std::array<int,2>> face_nodes;      //un-oriented list of nodes a---b
    std::vector<std::array<int,2>> face_elements;   //oriented face list (A -|-> B)
    std::vector<std::vector<int>> element_faces;    //unordered list of faces
    std::vector<std::vector<int>> element_neighbors; //list of nearest neighbors

    //gradient reconstruction vectors for solving least squares problem
    std::vector<Eigen::MatrixXd> A_e;
    std::vector<Eigen::MatrixXd> A_inv; //psuedo inverse
    std::vector<Eigen::VectorXd> b_e;

    virtual void init(Job* job, FiniteVolumeDriver* driver);

    virtual void writeHeader(std::ofstream& file, int TYPE);

    virtual double getElementVolume(int e);
    virtual int getElementTag(int e);
    virtual double getFaceArea(int f);
    virtual int getFaceTag(int f);
    virtual KinematicVector getFaceNormal(Job* job, int f);

    virtual std::vector<int> getElementFaces(int e);
    virtual std::array<int,2> getOrientedElementsByFace(int e);
    virtual std::vector<int> getElementNeighbors(int e);
    virtual KinematicVector getElementCentroid(Job* job, int e); //return element centroid
    virtual KinematicVector getFaceCentroid(Job* job, int f);

    virtual void generateMappings(Job* job, FiniteVolumeDriver* driver);
    virtual void constructMomentumField(Job* job, FiniteVolumeDriver* driver);
    virtual void constructDensityField(Job* job, FiniteVolumeDriver* driver);
    virtual KinematicTensorArray getVelocityGradients(Job* job, FiniteVolumeDriver* driver);

    //functions to compute element-wise fluxes of field variables using reconstructed velocity field
    virtual Eigen::VectorXd calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, Eigen::VectorXd& values);

    //functions to compute element mass flux
    virtual Eigen::VectorXd calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver);

    //functions to compute element momentum fluxes
    virtual KinematicVectorArray calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver);
};

#endif //MPM_V3_FVM_GRIDS_HPP

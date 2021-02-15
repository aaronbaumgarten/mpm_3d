//
// Created by aaron on 5/14/18.
// points.hpp
//

#ifndef MPM_V3_POINTS_HPP
#define MPM_V3_POINTS_HPP

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <regex>
#include <algorithm>
#include <sstream>
#include <Eigen/Core>
#include <ctime>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

/*
 * IN THIS FILE, DEFINE POINTS OBJECTS.
 * EACH OBJECT MUST BE ADDED TO THE REGISTRY IN src/registry
 * BEFORE USE.
 */

/*
 * class Points : public MPMObject{
public:
    std::string file;

    KinematicVectorArray x, u, x_t, mx_t, b;    //state vectors
    KinematicTensorArray L;                     //velocity gradient
    MaterialTensorArray T;                      //cauchy stress
    Eigen::VectorXd m, v, v0, extent;           //weight measures
    Eigen::VectorXi active;                     //active?

    virtual void init(Job*, Body*) = 0;                         //initialize from Job and Body
    virtual void readFromFile(Job*, Body*, std::string) = 0;    //construct points from given file
    virtual void generateMap(Job*, Body*, int) = 0;             //generate S and gradS
    virtual void updateIntegrators(Job*, Body*) = 0;            //update integrators (extent, etc.)

    virtual void writeHeader(Job*, Body*, Serializer*, std::ofstream&, int) = 0;
    virtual void writeFrame(Job*, Body*, Serializer*) = 0;                      //send frame data to Serializer
    virtual std::string saveState(Job*, Body*, Serializer*, std::string) = 0;   //save to file (in given directory)
    virtual int loadState(Job*, Body*, Serializer*, std::string) = 0;           //load data from full path

    virtual void generateLoads(Job*, Body*) = 0;   //arbitrary loading during simulation
    virtual void applyLoads(Job*, Body*) = 0;
};
 */

class DefaultPoints : public Points{
public:
    DefaultPoints(){
        object_name = "DefaultPoints";
    }

    Eigen::VectorXd extent;
    Eigen::MatrixXi A; //for mapping corners to position

    bool use_elem;          //determines whether to use element history for points
    Eigen::MatrixXi elem;   //for holding previous elemental positions (to limit cell searches to cell crossings)

    virtual void init(Job* job, Body* body);
    virtual void readFromFile(Job* job, Body* body, std::string fileIN);
    virtual void generateMap(Job* job, Body* body, int SPEC);             //generate S and gradS
    virtual void updateIntegrators(Job* job, Body* body);                 //update integrators (extent, etc.)

    virtual void writeHeader(Job* job, Body* body, Serializer* serializer, std::ofstream& pfile, int SPEC);
    virtual void writeFrame(Job* job, Body* body, Serializer* serializer);
    virtual std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    virtual int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);

    virtual void generateLoads(Job* job, Body* body);   //arbitrary loading during simulation
    virtual void applyLoads(Job* job, Body* body);
};

/*----------------------------------------------------------------------------*/

//general object for creating sets of points
class Part : public MPMObject{
public:
    virtual bool encompasses(KinematicVector&) = 0;
    virtual void init(Job*) = 0;
};

class Ball : public Part{
public:
    Ball(){
        object_name = "Ball";
    }

    double r;           //radius
    KinematicVector o;  //origin

    void init(Job* job);
    bool encompasses(KinematicVector& xIN);
};

class Box : public Part{
public:
    Box(){
        object_name = "Box";
    }

    KinematicVector x_min, x_max;  //bounds

    void init(Job* job);
    bool encompasses(KinematicVector& xIN);
};

class SineWave : public Part{
public:
    SineWave(){
        object_name = "SineWave";
    }

    double amplitude, wavelength, x_0, y_0;

    void init(Job* job);
    bool encompasses(KinematicVector& xIN);

};

class SandPile : public Part{
public:
    SandPile(){
        object_name = "SandPile";
    }

    double mu;
    int H_index;
    KinematicVector peak;

    void init(Job* job);
    bool encompasses(KinematicVector& xIN);

};

/*----------------------------------------------------------------------------*/

class GmshPoints : public DefaultPoints{
public:
    GmshPoints(){
        object_name = "GmshPoints";
    }

    double msh_version;

    //linear material points per cell and density of points
    double rho;
    int lmpp;

    //part list
    std::vector<std::unique_ptr<Part>> part_list;

    //msh file
    std::string msh_file, out_file;

    //void init(Job* job, Body* body);
    void readFromFile(Job* job, Body* body, std::string fileIN);
};

/*----------------------------------------------------------------------------*/

class CartesianPoints : public DefaultPoints{
public:
    CartesianPoints(){
        object_name = "CartesianPoints";
    }

    KinematicVector Lx;
    Eigen::VectorXi Nx;

    //linear material points per cell and density of points
    double rho;
    int lmpp;

    //part list
    std::vector<std::unique_ptr<Part>> part_list;

    //msh file
    std::string msh_file, out_file;

    //void init(Job* job, Body* body);
    void readFromFile(Job* job, Body* body, std::string fileIN);
};

/*---------------------------------------------------------------------------*/

class ThreadPoolPoints : public DefaultPoints{
public:
    ThreadPoolPoints(){
        object_name = "ThreadPoolPoints";
    }

    std::vector<MPMScalarSparseMatrix> S_vec = std::vector<MPMScalarSparseMatrix>(0);
    std::vector<KinematicVectorSparseMatrix> gradS_vec = std::vector<KinematicVectorSparseMatrix>(0);

    static void fillMap(Job* job, Body* body, ThreadPoolPoints* points,
                        int SPEC,
                        MPMScalarSparseMatrix& S,
                        KinematicVectorSparseMatrix& gradS,
                        int i_begin, int i_end, volatile bool &done);

    static void combineMaps(Job* job, Body* body, ThreadPoolPoints* points,
                            MPMScalarSparseMatrix& S,
                            KinematicVectorSparseMatrix& gradS,
                            int i_begin_S, int i_begin_gradS, volatile bool &done);

    //virtual void init(Job* job, Body* body);
    //virtual void readFromFile(Job* job, Body* body, std::string fileIN);
    virtual void generateMap(Job* job, Body* body, int SPEC);             //generate S and gradS
    //virtual void updateIntegrators(Job* job, Body* body);                 //update integrators (extent, etc.)

    //virtual void writeHeader(Job* job, Body* body, Serializer* serializer, std::ofstream& pfile, int SPEC);
    //virtual void writeFrame(Job* job, Body* body, Serializer* serializer);
    //virtual std::string saveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
    //virtual int loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);

    //virtual void generateLoads(Job* job, Body* body);   //arbitrary loading during simulation
    //virtual void applyLoads(Job* job, Body* body);
};


/*----------------------------------------------------------------------------*/

class ImprovedQuadraturePoints : public DefaultPoints{
public:
    ImprovedQuadraturePoints(){
        object_name = "ImprovedQuadraturePoints";
    }

    /* input to this file denotes the type of improved quadrature to use
     * 0 -- No quadrature improvement
     * 1 -- uGIMP
     * 2 -- cpGIMP
     * 3 -- CPDI
     * 4 -- CPDI2
     */

    /* second input to this file denotes the type of position correction to use
     * 1 -- Avoid-a-void
     * 2 -- SPH correction
     * 3 -- \delta-correction, strain
     * 4 -- \delta-correction, displacement
     */

    static const int STANDARD = 0;
    static const int UGIMP = 1;
    static const int CPGIMP = 2;
    static const int CPDI = 3;
    static const int CPDI2 = 4;

    static const int NO_CORRECTION = 0;
    static const int AVAV = 1;
    static const int SPH_LIKE = 2;
    static const int DELTA_STRAIN = 3;
    static const int DELTA_DISP = 4;
    static const int DELTA_NEWTON = 5;

    int QUADRULE = 0;
    int POSITIONRULE = 0;

    //uGIMP already implemented
    //cpGIMP requires adjusting F on each point
    //CPDI requires only F
    //CPDI2 requires corner positions be advected

    //error measure, true integral of node basis function, estimated integral of node basis function
    Eigen::VectorXd e, H, V_i, v_i, eonV;

    //gradient of error measure
    KinematicVectorArray grad_e, grad_H, grad_eonV;

    //deformation gradient (points)
    KinematicTensorArray F;

    //number of corners per material point
    int cpmp;

    //list of corner positions (length cpmp*points)
    KinematicVectorArray corner_positions;

    //grid variables for delta correction scheme
    double alpha, h, REL_TOL;
    KinematicVectorArray del_pos;
    int max_iter;
    bool use_velocity_limiter = false;
    bool use_strain_limiter = false;
    double limiter = 1.0;

    //variables for avoid-a-void scheme
    double alg_inf = 1e10; // a really big number
    double buffer_scale, r, merge_dist;
    int skip_counter, skip_value;
    Eigen::VectorXd d;                      //grid point distance f'n value
    std::vector<int> edge_list;             //list of indices associated with each edge
    std::vector<int> buffer_list;           //list of points in buffer (inactive)

    //avoid-a-void search cell definitions and functions
    int number_of_neighbors = 3;
    int number_of_second_neighbors = 5;
    KinematicVector x_min, x_max, Lx, hx;
    Eigen::VectorXi Nx;
    std::vector<int> node_to_cell_map;
    std::vector<int> point_to_cell_map;
    std::vector<std::vector<int>> cell_to_point_map;
    std::vector<std::vector<int>> cell_to_node_map;
    Eigen::VectorXd point_dist;

    int pos_to_cell(Job* job, KinematicVector x);  //function to identify which cell 'x' is in
    std::vector<int> cell_to_ijk(Job* job, int i); //function to switch between cell id's and ijk position
    int ijk_to_cell(Job* job, std::vector<int> ijk);   //function to switch between ijk positions and cell id

    //sph-like needs additional contasnt
    double sph_const = 0.1;
    double u_max;

    //internal variables for i/o
    int sampledFrames = 0;
    std::string outputFolder;

    //initialization function
    void init(Job* job, Body* body);

    //generate mapping matrix associated with point integration
    void generateMap(Job* job, Body* body, int SPEC);

    //update integrators (F, corner_positions)
    void updateIntegrators(Job* job, Body* body);

    //output corners to file (maybe vtk?)
    void writeFrame(Job* job, Body* body, Serializer* serializer);
};

/*----------------------------------------------------------------------------*/
/*
class EnhancedStrainPoints : public DefaultPoints{
public:
    EnhancedStrainPoints(){
        object_name = "EnhancedStrainPoints";
    }

    Eigen::VectorXd e, V_i, v_i;
    KinematicVectorArray grad_e;
    double alpha, h;

    void init(Job* job, Body* body);

    void generateLoads(Job* job, Body* body);
    void applyLoads(Job* job, Body* body);

    void writeFrame(Job* job, Body* body, Serializer* serializer);
};
 */

#endif //MPM_V3_POINTS_HPP

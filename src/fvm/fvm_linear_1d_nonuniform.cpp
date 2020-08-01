//
// Created by aaron on 12/12/19.
// fvm_cartesian.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <Eigen/Dense>
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
#include "fvm_grids.hpp"
#include "threadpool.hpp"
#include "../objects/grids/grids.hpp"

//helper functions
int FVMLinear1DNonUniform::whichElement(Job* job, KinematicVector xIN){
    //convert Lx, Nx, Hx to KinematicVector
    KinematicVector lx = KinematicVector(Lx.data(), job->JOB_TYPE);
    KinematicVector dx = KinematicVector(hx.data(), job->JOB_TYPE);
    Eigen::VectorXi nx = Eigen::VectorXi(GRID_DIM);
    for (int i=0; i<GRID_DIM; i++){
        nx(i) = Nx[i];
    }
    return CartesianLinear::cartesianWhichElement(job, xIN, lx, dx, nx);
} //return which element this position is within

std::vector<double> FVMLinear1DNonUniform::getBoundingBox(){
    std::vector<double> result = std::vector<double>(2*GRID_DIM);
    for (int i=0; i<GRID_DIM; i++){
        result[2*i] = 0.0;
        result[2*i+1] = Lx[i];
    }
    return result;
}  //return min/max of each coordinate


void FVMLinear1DNonUniform::init(Job* job, FiniteVolumeDriver* driver){
    if (job->JOB_TYPE != job->JOB_1D){
        std::cout << "ERROR! FVMLinear1DNonUniform grid designed for 1D simulations only! Exiting." << std::endl;
        exit(0);
    }

    GRID_DIM = 1;

    //call initializer for base class
    FVMGridBase::init(job, driver);

    //check size of properties passed to grid object
    if (fp64_props.size() < GRID_DIM || int_props.size() < GRID_DIM) {
        std::cout << fp64_props.size() << ", " << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least %i properties defined (Lx, Nx).\n",
                __FILE__, __func__, (2*GRID_DIM));
        exit(0);
    } else {
        //set side lengths and discretization
        Lx = std::vector<double>(GRID_DIM);
        Nx = std::vector<int>(GRID_DIM);
        hx = std::vector<double>(GRID_DIM);
        for (int pos=0; pos<GRID_DIM; pos++){
            Lx[pos] = fp64_props[pos];
            Nx[pos] = int_props[pos];
            hx[pos] = Lx[pos]/Nx[pos];
        }
    }

    //initialize tags and values vectors;
    std::vector<FiniteVolumeMethod::BCContainer> tmp_bc_info = std::vector<FiniteVolumeMethod::BCContainer>(2*GRID_DIM);

    if (int_props.size() == GRID_DIM){
        //no bc tags given
        //set to zero dirichlet
        if (driver->TYPE == FiniteVolumeDriver::ISOTHERMAL) {
            for (int i = 0; i < tmp_bc_info.size(); i++) {
                tmp_bc_info[i].tag = FiniteVolumeGrid::VELOCITY_INLET;
                tmp_bc_info[i].values = std::array<double, 2>();
                tmp_bc_info[i].vector = KinematicVector(job->JOB_TYPE);
                tmp_bc_info[i].vector.setZero();
            }
        } else if (driver->TYPE == FiniteVolumeDriver::THERMAL) {
            for (int i = 0; i < tmp_bc_info.size(); i++) {
                tmp_bc_info[i].tag = FiniteVolumeGrid::ADIABATIC_WALL;
                tmp_bc_info[i].values = std::array<double, 2>();
                tmp_bc_info[i].vector = KinematicVector(job->JOB_TYPE);
                tmp_bc_info[i].vector.setZero();
            }
        } else {
            std::cerr << "ERROR: FVMLinear1DNonUniform is not implemented for FiniteVolumeDriver TYPE " << driver->TYPE << "! Exiting.";
            exit(0);
        }
    } else {
        //counter to avoid overflow error
        int fp64_iterator = GRID_DIM; //already read in GRID_DIM properties

        //bc_tags given
        for (int i=0; i<tmp_bc_info.size(); i++){
            tmp_bc_info[i].tag = int_props[i+GRID_DIM];
            tmp_bc_info[i].values = std::array<double, 2>();
            tmp_bc_info[i].values[0] = 0; tmp_bc_info[i].values[1] = 0;
            tmp_bc_info[i].vector = KinematicVector(job->JOB_TYPE);
            //check that bc tag is in coded list
            if (tmp_bc_info[i].tag != VELOCITY_INLET &&
                    tmp_bc_info[i].tag != VELOCITY_TEMP_INLET &&
                    tmp_bc_info[i].tag != VELOCITY_DENSITY_INLET &&
                    tmp_bc_info[i].tag != PRESSURE_INLET &&
                    tmp_bc_info[i].tag != PRESSURE_OUTLET &&
                    tmp_bc_info[i].tag != DAMPED_OUTLET &&
                    tmp_bc_info[i].tag != ADIABATIC_WALL &&
                    tmp_bc_info[i].tag != THERMAL_WALL &&
                    tmp_bc_info[i].tag != SYMMETRIC_WALL &&
                    tmp_bc_info[i].tag != SUPERSONIC_INLET &&
                    tmp_bc_info[i].tag != SUPERSONIC_OUTLET &&
                    tmp_bc_info[i].tag != PERIODIC &&
                    tmp_bc_info[i].tag != DAMPED_WALL &&
                    tmp_bc_info[i].tag != STAGNATION_INLET){
                std::cerr << "ERROR: Boundary tag " << tmp_bc_info[i].tag << " not defined for FVMLinear1DNonUniform grid object! Exiting." << std::endl;
                exit(0);
            }

            //add required input parameters for boundary type
            switch (tmp_bc_info[i].tag){
                case VELOCITY_INLET:
                    //check that boundary condition is valid for this simulation type
                    if (driver->TYPE != FiniteVolumeDriver::ISOTHERMAL){
                        std::cerr << "ERROR: VELOCITY_INLET boundary condition requires ISOTHERMAL simulation type for stability. Exiting.";
                    } else {
                        //first GRID_DIM properties are the components of velocity at boundary
                        for (int pos=0; pos<GRID_DIM; pos++){
                            if (fp64_props.size() > fp64_iterator) {
                                tmp_bc_info[i].vector[pos] = fp64_props[fp64_iterator];
                                fp64_iterator++;
                            } else {
                                std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                                exit(0);
                            }
                        }
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : VELOCITY_INLET : u = ";
                    std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(tmp_bc_info[i].vector).transpose() << std::endl;
                    break;

                case VELOCITY_TEMP_INLET:
                    //check that boundary condition is valid for this simulation type
                    if (driver->TYPE != FiniteVolumeDriver::THERMAL){
                        std::cerr << "ERROR: VELOCITY_TEMP_INLET boundary condition requires THERMAL simulation type. Exiting.";
                    } else {
                        //first GRID_DIM properties are the components of velocity at boundary
                        for (int pos=0; pos<GRID_DIM; pos++){
                            if (fp64_props.size() > fp64_iterator) {
                                tmp_bc_info[i].vector[pos] = fp64_props[fp64_iterator];
                                fp64_iterator++;
                            } else {
                                std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                                exit(0);
                            }
                        }
                        //next property is temperature
                        if (fp64_props.size() > fp64_iterator) {
                            tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                            fp64_iterator++;
                        } else {
                            std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                            exit(0);
                        }
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : VELOCITY_TEMP_INLET : u = ";
                    std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(tmp_bc_info[i].vector).transpose();
                    std::cout << ", T = " << tmp_bc_info[i].values[0] << std::endl;
                    break;

                case VELOCITY_DENSITY_INLET:
                    //first GRID_DIM properties are the components of velocity at boundary
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (fp64_props.size() > fp64_iterator) {
                            tmp_bc_info[i].vector[pos] = fp64_props[fp64_iterator];
                            fp64_iterator++;
                        } else {
                            std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                            exit(0);
                        }
                    }
                    //next property is density
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : VELOCITY_DENSITY_INLET : u = ";
                    std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(tmp_bc_info[i].vector).transpose();
                    std::cout << ", rho = " << tmp_bc_info[i].values[0] << std::endl;
                    break;

                case PRESSURE_INLET:
                    //first property is density
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }
                    //second property is temperature
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[1] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : PRESSURE_INLET : P = " << tmp_bc_info[i].values[0];
                    std::cout << ", T = " << tmp_bc_info[i].values[1] << std::endl;
                    break;

                case PRESSURE_OUTLET:
                    //first property is density
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }
                    //second property is temperature
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[1] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : PRESSURE_OUTLET : P = " << tmp_bc_info[i].values[0];
                    std::cout << ", T* = " << tmp_bc_info[i].values[1] << std::endl;
                    break;

                case DAMPED_OUTLET:
                    //first property is density
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }
                    //second property is temperature
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[1] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : DAMPED_OUTLET : P = " << tmp_bc_info[i].values[0];
                    std::cout << ", T* = " << tmp_bc_info[i].values[1] << std::endl;
                    break;

                case ADIABATIC_WALL:
                    //for readability, take in one unused property
                    fp64_iterator++;

                    //print boundary condition info
                    std::cout << " - " << i << " : ADIABATIC_WALL" << std::endl;
                    break;

                case THERMAL_WALL:
                    //one property, temperature
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : THERMAL_WALL : T = " << tmp_bc_info[i].values[0] << std::endl;
                    break;

                case SYMMETRIC_WALL:
                    //for readability, take in one unused property
                    fp64_iterator++;

                    //print boundary condition info
                    std::cout << " - " << i << " : SYMMETRIC_WALL" << std::endl;
                    break;

                case SUPERSONIC_INLET:
                    //first GRID_DIM properties are the components of velocity at boundary
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (fp64_props.size() > fp64_iterator) {
                            tmp_bc_info[i].vector[pos] = fp64_props[fp64_iterator];
                            fp64_iterator++;
                        } else {
                            std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                            exit(0);
                        }
                    }
                    //next property is density
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }
                    //final property is temperature
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[1] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : SUPERSONIC_INLET : u = ";
                    std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(tmp_bc_info[i].vector).transpose();
                    std::cout << ", rho = " << tmp_bc_info[i].values[0];
                    std::cout << ", T = " << tmp_bc_info[i].values[1] << std::endl;
                    break;

                case SUPERSONIC_OUTLET:
                    //for readability, take in one unused property
                    fp64_iterator++;

                    //print boundary condition info
                    std::cout << " - " << i << " : SUPERSONIC_OUTLET" << std::endl;
                    break;

                case PERIODIC:
                    //for readability, take in one unused property
                    fp64_iterator++;

                    //print boundary condition info
                    std::cout << " - " << i << " : PERIODIC" << std::endl;
                    break;

                case DAMPED_WALL:
                    //one property, damping coefficient
                    if (fp64_props.size() > fp64_iterator) {
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : DAMPED_WALL : nu = " << tmp_bc_info[i].values[0] << std::endl;
                    break;

                case STAGNATION_INLET:
                    //two properties: P^t and T^t
                    if (fp64_props.size() > fp64_iterator+1){
                        tmp_bc_info[i].values[0] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                        tmp_bc_info[i].values[1] = fp64_props[fp64_iterator];
                        fp64_iterator++;
                    } else {
                        std::cerr << "ERROR: Not enough fp64 properties given. Exiting." << std::endl;
                        exit(0);
                    }

                    //print boundary condition info
                    std::cout << " - " << i << " : STAGNATION_INLET : P^t = " << tmp_bc_info[i].values[0];
                    std::cout << ", T^t = " << tmp_bc_info[i].values[1] << std::endl;
                    break;

                default:
                    //do nothing
                    break;
            }

        }

        //check bc_tags for consistency
        for (int i=0; i<GRID_DIM; i++){
            if ((tmp_bc_info[2*i].tag == FiniteVolumeGrid::PERIODIC
                 || tmp_bc_info[2*i+1].tag == FiniteVolumeGrid::PERIODIC)
                   && tmp_bc_info[2*i].tag != tmp_bc_info[2*i+1].tag){
                std::cerr << "ERROR: Boundary conditions defined on FVMLinear1DNonUniform grid do not match! Exiting." << std::endl;
                exit(0);
            }
        }
    }

    //initialize grid information
    if (GRID_DIM == 1){
        //number of grid elements
        element_count = Nx[0];
        node_count = Nx[0] + 1;
        face_count = Nx[0] + 1;

        //quadrature rule
        if (driver->ORDER == 1 || USE_REDUCED_QUADRATURE){
            qpe = 1; //quad points per element
            qpf = 1; //quad points per face
        } else if (driver->ORDER == 2){
            qpe = 2;
            qpf = 1;
        } else {
            std::cerr << "ERROR: FVMLinear1DNonUniform not defined for simulation ORDER > 2." << std::endl;
            qpe = 2;
            qpf = 1;
        }
        ext_quad_count = face_count*qpf;
        int_quad_count = element_count*qpe;

        //volumes and areas
        v_e = Eigen::VectorXd(element_count);
        v_e.setConstant(hx[0]);
        face_areas = Eigen::VectorXd(face_count);
        face_areas.setConstant(1.0);

        //face normals
        face_normals = KinematicVectorArray(face_count, job->JOB_TYPE);
        for (int f=0; f<face_count; f++){
            face_normals(f,0) = 1; //+x
        }

        //define default face to element definitions.
        face_elements = std::vector<std::array<int,2>>(face_count);
        face_elements[0][0] = -1;
        face_elements[0][1] = 0;
        face_elements[face_count-1][0] = element_count-1;
        face_elements[face_count-1][1] = -1;
        for (int i=1; i<face_count-1; i++){
            face_elements[i][0] = i-1;
            face_elements[i][1] = i;
        }

        //define which boundary is associated with each face
        bc_info = std::vector<FiniteVolumeMethod::BCContainer>(face_count);
        //initalize bc_info
        for (int f=0; f<face_count; f++){
            bc_info[f].tag = -1;
            bc_info[f].values = std::array<double,2>();
            bc_info[f].values[0] = 0;
            bc_info[f].values[1] = 0;
            bc_info[f].vector = KinematicVector(job->JOB_TYPE);
            bc_info[f].vector.setZero();
        }

        bc_info[0]            = tmp_bc_info[0]; //-x
        bc_info[face_count-1] = tmp_bc_info[1]; //+x
        for (int i=1; i<face_count-1; i++){
            bc_info[i].tag = -1;
        }

        //specify boundary behavior on periodic ends
        if (bc_info[0].tag == FiniteVolumeGrid::PERIODIC){
            //remove right face and replace left face
            face_elements[0][0] = element_count-1;
            face_elements[face_count-1][0] = -1;
        }

    } else {
        std::cout << "Uh?" << std::endl;
    }

    //define face, element, and quadrature centroids
    x_e = KinematicVectorArray(element_count, job->JOB_TYPE);
    x_f = KinematicVectorArray(face_count, job->JOB_TYPE);
    x_q = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE); //volume integrals first, then surface
    w_q = Eigen::VectorXd(int_quad_count + ext_quad_count);
    q_b = std::vector<bool>(int_quad_count + ext_quad_count);

    //initialize porosity and solid velocity fields
    n_q = Eigen::VectorXd(int_quad_count + ext_quad_count);
    n_q.setConstant(1.0);
    gradn_q = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE);
    gradn_q.setZero();
    v_sq = KinematicVectorArray(int_quad_count + ext_quad_count, job->JOB_TYPE);
    v_sq.setZero();

    //quadrature points are default NOT on boundary
    for (int i=0; i<int_quad_count+ext_quad_count; i++){
        q_b[i] = false;
    }

    x_f(0,0) = 0.0;
    x_f(face_count-1,0) = Lx[0];
    //face centroid
    for (int f=0; f<face_count; f++) {
        x_f(f,0) = hx[0] * (f - 0.25 + 0.5*std::rand()/RAND_MAX);
    }

    //element centroid
    std::vector<int> tmp;
    for (int e=0; e<element_count; e++){
        x_e(e,0) = (x_f(e,0) + x_f(e+1,0))/2.0;
        v_e(e) = (x_f(e,0) - x_f(e+1,0));
    }

    //quadrature points
    double offset_factor = 1.0/sqrt(3.0);

    //loop over elements for interior quadrature points
    for (int e=0; e<element_count; e++){
        if (qpe == 1){
            //quad point collocated with element centroid
            x_q[e*qpe] = x_e[e];
            w_q(e*qpe) = v_e[e];
        } else if (qpe == 2){
            //quad points at +/- sqrt(1/3) along x-axis
            x_q(e*qpe + 0,0) = x_e(e,0) - offset_factor*v_e[0]/2.0;
            x_q(e*qpe + 1,0) = x_e(e,0) + offset_factor*v_e[0]/2.0;
            w_q[e*qpe + 0] = v_e[e]/qpe;
            w_q[e*qpe + 1] = v_e[e]/qpe;
        }
    }

    //loop over faces for exterior quadrature points
    for (int f=0; f<face_count; f++){
        for (int q=0; q<qpf; q++) {
            x_q[int_quad_count + f * qpf + q] = x_f[f];
            w_q[int_quad_count + f * qpf + q] = face_areas(f)/qpf;
            if (bc_info[f].tag > -1 && bc_info[f].tag != PERIODIC) {
                //face is bounding face and quadrature point should be included in boundary integrals
                q_b[int_quad_count + f * qpf + q] = true;
            }
        }
    }

    //with face to element map defined, we can construct the element to face map robustly as follows:
    element_faces = std::vector<std::vector<int>>(element_count);
    for (int f=0; f<face_count; f++){
        //loop over each face
        if (face_elements[f][0] >= 0){
            //if A is an element, f is a face of A
            element_faces[face_elements[f][0]].push_back(f);
        }
        if (face_elements[f][1] >= 0){
            //if B is an element, f is a face of B
            element_faces[face_elements[f][1]].push_back(f);
        }
    }

    //consistency check!!!
    for (int e=0; e<element_count; e++){
        if (element_faces[e].size() != 2*GRID_DIM){
            std::cout << "ERROR: Element " << e << " has wrong number of faces! ";
            std::cout << element_faces[e].size() << " != " << 2*GRID_DIM << "." << std::endl;
        }
    }

    //generate map of element neighbors which share a face or node
    element_neighbors = std::vector<std::vector<int>>(element_count);
    std::vector<int> ijk = std::vector<int>(GRID_DIM);
    //std::vector<int> tmp = ijk;
    tmp = ijk;
    int tmp_e;
    for (int e=0; e<element_count; e++){
        ijk = e_to_ijk(e);
        tmp = ijk;
        //loop over i
        for (int i=-1; i<2; i+=2){
            tmp[0] = ijk[0] + i;
            //check for periodic BC
            for (int f=0; f<element_faces[e].size(); f++){
                //loop -x
                if (face_normals(f,0) == 1 && bc_info[f].tag == PERIODIC){
                    if (ijk[0] == 0 && i==-1) {
                        tmp[0] = Nx[0] - 1;
                    } else if (ijk[0] == Nx[0]-1 && i==1){
                        tmp[0] = 0;
                    }
                }
            }

            tmp_e = ijk_to_e(tmp);
            if (tmp_e >= 0 && tmp_e != e && tmp_e < element_count){
                element_neighbors[e].push_back(tmp_e);
            }
        }

        if (GRID_DIM > 1) {
            tmp = ijk;
            //loop over j
            for (int j = -1; j < 2; j += 2) {
                tmp[1] = ijk[1] + j;
                //check for periodic BC
                for (int f = 0; f < element_faces[e].size(); f++) {
                    //loop -y
                    if (face_normals(f, 1) == 1 && bc_info[f].tag == PERIODIC) {
                        if (ijk[1] == 0 && j == -1) {
                            tmp[1] = Nx[1] - 1;
                        } else if (ijk[1] == Nx[1] - 1 && j == 1) {
                            tmp[1] = 0;
                        }
                    }
                }

                tmp_e = ijk_to_e(tmp);
                if (tmp_e >= 0 && tmp_e != e && tmp_e < element_count) {
                    element_neighbors[e].push_back(tmp_e);
                }
            }
        }

        if (GRID_DIM > 2){
            tmp = ijk;
            for (int k=-1; k<2; k+=2){
                //loop over k
                tmp[2] = ijk[2] + k;
                //check for periodic BC
                for (int f=0; f<element_faces[e].size(); f++){
                    //loop -z
                    if (face_normals(f,2) == 1 && bc_info[f].tag == PERIODIC){
                        if (ijk[2] == 0 && k == -1) {
                            tmp[2] = Nx[2] - 1;
                        } else if (ijk[2] == Nx[2] - 1 && k == 1) {
                            tmp[2] = 0;
                        }
                    }
                }

                tmp_e = ijk_to_e(tmp);
                if (tmp_e >= 0 && tmp_e != e && tmp_e < element_count){
                    element_neighbors[e].push_back(tmp_e);
                }
            }
        }

        /*
        std::cout << "[" << e << "] : ";
        for (int ii=0; ii<element_neighbors[e].size(); ii++){
            std::cout << element_neighbors[e][ii] << ", ";
        }
        std::cout << std::endl;
         */
    }

    //consistency check!!!
    int num_neighbors = 3;
    for (int i = 1; i<GRID_DIM; i++){
        num_neighbors *= 3;
    }
    num_neighbors -= 1;

    for (int e=0; e<element_count; e++){
        if (element_neighbors[e].size() > num_neighbors){
            std::cout << "ERROR: Element " << e << " has too many neighbors! ";
            std::cout << element_neighbors[e].size() << " > " << num_neighbors - 1 << "!" << std::endl;
        }
    }

    //lastly, construct A and b for least squares calculations
    A_e = std::vector<Eigen::MatrixXd>(element_count);
    A_inv = std::vector<Eigen::MatrixXd>(element_count);
    b_e = std::vector<Eigen::VectorXd>(element_count);
    int length_of_A; //including element neighbors and dirichlet faces
    //initialize matrix A and vector B
    KinematicVector x_0, x;
    double tmp_dif;
    for (int e = 0; e<element_count; e++) {
        //size A and b for element based on number of neighboring cells and dirichlet faces
        length_of_A = element_neighbors[e].size();

        //get element centroid
        x_0 = getElementCentroid(job, e);
        A_e[e] = Eigen::MatrixXd(length_of_A,GRID_DIM);
        b_e[e] = Eigen::VectorXd(length_of_A);
        //create system of equations
        for (int ii = 0; ii < element_neighbors[e].size(); ii++) {
            x = getElementCentroid(job, element_neighbors[e][ii]);
            for (int pos = 0; pos < GRID_DIM; pos++) {
                tmp_dif = x[pos] - x_0[pos];
                //if periodic, then x - x_0 may be greater than L/2
                if (tmp_dif > Lx[pos] / 2.0) {
                    A_e[e](ii, pos) = tmp_dif - Lx[pos];
                } else if (tmp_dif < -Lx[pos] / 2.0) {
                    A_e[e](ii, pos) = tmp_dif + Lx[pos];
                } else {
                    A_e[e](ii, pos) = x[pos] - x_0[pos];
                }
            }
        }

        //psuedo inverse
        A_inv[e] = Eigen::MatrixXd(GRID_DIM,A_e[e].rows());
        Eigen::MatrixXd AtA_inv = Eigen::MatrixXd(GRID_DIM,GRID_DIM);
        AtA_inv.setZero();
        Eigen::MatrixXd AtA = A_e[e].transpose()*A_e[e];
        Eigen::VectorXd a;
        Eigen::VectorXd e_vec = Eigen::VectorXd(AtA.rows());
        for (int ii=0; ii<AtA.rows(); ii++){
            e_vec.setZero();
            e_vec(ii) = 1;
            //a = AtA.householderQr().solve(e_vec);
            a = AtA.fullPivHouseholderQr().solve(e_vec);
            for (int jj=0; jj<GRID_DIM; jj++){
                AtA_inv(jj,ii) = a(jj);
            }
        }
        A_inv[e] = AtA_inv*A_e[e].transpose();
    }

    //initialize grid mappings
    generateMappings(job, driver);

    //print grid properties
    std::cout << "FiniteVolumeGrid properties:" << std::endl;
    std::cout << "    Lx =";
    for (int pos=0; pos<Lx.size(); pos++){
        std::cout << " " << Lx[pos];
    }
    std::cout << std::endl << "    Nx =";
    for (int pos=0; pos<Lx.size(); pos++){
        std::cout << " " << Nx[pos];
    }
    std::cout << std::endl;

    std::cout << "    Elements: " << element_count << std::endl;
    std::cout << "    Faces:    " << face_count << std::endl;
    std::cout << "    Quad Points: " << int_quad_count + ext_quad_count << std::endl;

    std::cout << "FiniteVolumeGrid initialized." << std::endl;
}

/*----------------------------------------------------------------------------*/

std::vector<int> FVMLinear1DNonUniform::e_to_ijk(int e){
    std::vector<int> ijk = std::vector<int>(GRID_DIM);
    int tmp = e;
    for (int i=0;i<ijk.size();i++){
        ijk[i] = tmp % Nx[i];
        tmp = tmp/Nx[i];
    }
    return ijk;
}

int FVMLinear1DNonUniform::ijk_to_e(std::vector<int> ijk){
    int tmp = 0;
    for (int i=ijk.size(); i>0; i--){
        if (ijk[i-1] < 0 || ijk[i-1] >= Nx[i-1]){
            return -1;
        }
        tmp = tmp*Nx[i-1] + ijk[i-1];
    }
    return tmp;
}

std::vector<int> FVMLinear1DNonUniform::n_to_ijk(int n){
    std::vector<int> ijk = std::vector<int>(GRID_DIM);
    int tmp = n;
    for (int i=0;i<ijk.size();i++){
        ijk[i] = tmp % (Nx[i]+1);
        tmp = tmp/(Nx[i]+1);
    }
    return ijk;
}

int FVMLinear1DNonUniform::ijk_to_n(std::vector<int> ijk){
    int tmp = 0;
    for (int i=ijk.size(); i>0; i--){
        if (ijk[i-1] < 0 || ijk[i-1] >= (Nx[i-1]+1)){
            return -1;
        }
        tmp = tmp*(Nx[i-1]+1) + ijk[i-1];
    }
    return tmp;
}

/*---------------------------------------------------------------------------*/

void FVMLinear1DNonUniform::writeHeader(std::ofstream& file, int SPEC){
    if (SPEC != Serializer::VTK){
        std::cerr << "ERROR: Unknown file SPEC in writeHeader: " << SPEC  << "! Exiting." << std::endl;
        exit(0);
    }

    int nlen = node_count;

    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    file << "POINTS " << nlen << " double\n";

    std::vector<int> ijk, tmp;
    for (int i=0;i<nlen;i++){
        //vtk files require x,y,z
        ijk = n_to_ijk(i);
        for (int pos = 0; pos < 3; pos++){
            if (pos < GRID_DIM){
                file << x_f(i,pos) << " ";
            } else {
                file << "0 ";
            }
        }
        file << "\n";
    }

    if (GRID_DIM == 1) {
        //use lines
        file << "CELLS " << element_count << " " << 3 * element_count << "\n";
        for (int e = 0; e < element_count; e++) {
            ijk = e_to_ijk(e);
            tmp = ijk; tmp[0] += 1;
            file << "2 " << ijk_to_n(ijk) << " " << ijk_to_n(tmp) << "\n";
        }

        file << "CELL_TYPES " << element_count << "\n";
        for (int e = 0; e < element_count; e++) {
            file << "3\n";
        }
    } else if (GRID_DIM == 2){
        file << "CELLS " << element_count << " " << 5 * element_count << "\n";
        for (int e = 0; e < element_count; e++){
            ijk = e_to_ijk(e);
            file << "4 " << ijk_to_n(ijk);
            ijk[0] += 1;
            file << " " << ijk_to_n(ijk);
            ijk[0] -= 1; ijk[1] += 1;
            file << " " << ijk_to_n(ijk);
            ijk[0] += 1;
            file << " " << ijk_to_n(ijk) << "\n";
        }

        file << "CELL_TYPES " << element_count << "\n";
        for (int e=0; e<element_count; e++){
            file << "8\n";
        }
    } else if (GRID_DIM == 3){
        file << "CELLS " << element_count << " " << 9 * element_count << "\n";
        for (int e = 0; e < element_count; e++){
            file << "8 ";
            ijk = e_to_ijk(e);
            file << ijk_to_n(ijk) << " ";
            ijk[0] += 1;
            file << ijk_to_n(ijk) << " ";
            ijk[0] -= 1; ijk[1] += 1;
            file << ijk_to_n(ijk) << " ";
            ijk[0] += 1;
            file << ijk_to_n(ijk) << " ";
            ijk[0] -= 1; ijk[1] -= 1; ijk[2] += 1;
            file << ijk_to_n(ijk) << " ";
            ijk[0] += 1;
            file << ijk_to_n(ijk) << " ";
            ijk[0] -= 1; ijk[1] += 1;
            file << ijk_to_n(ijk) << " ";
            ijk[0] += 1;
            file << ijk_to_n(ijk) << "\n";
        }

        file << "CELL_TYPES " << element_count << "\n";
        for (int e=0; e<element_count; e++){
            file << "11\n";
        }
    }

    file << "CELL_DATA " << element_count << "\n";
    return;
}

/*----------------------------------------------------------------------------*/
void FVMLinear1DNonUniform::constructMomentumField(Job* job, FiniteVolumeDriver* driver){
    if (driver->ORDER == 1){
        driver->fluid_body->p_x.setZero(); //let momentum be constant within an element
    } else if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMLinear1DNonUniform does not currently implement higher ORDER momentum reconstruction."
                      << std::endl;
        }

        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        KinematicVector p_0, p, p_max, p_min, min_dif;
        min_dif = KinematicVector(job->JOB_TYPE);
        double tmp_dif, rho_0;
        for (int e = 0; e<element_count; e++){
            //least squares fit of u_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            p_0 = driver->fluid_body->p[e] / driver->fluid_body->n_e(e);
            rho_0 = driver->fluid_body->rho(e) / driver->fluid_body->n_e(e);
            p_max = p_0; p_min = p_0;
            min_dif.setZero();
            //loop over components of momentum
            for (int mom_index = 0; mom_index<GRID_DIM; mom_index++){
                //create system of equations
                for (int ii=0; ii<element_neighbors[e].size(); ii++){
                    //just fill in b vector
                    p = driver->fluid_body->p[element_neighbors[e][ii]] / driver->fluid_body->n_e(element_neighbors[e][ii]);
                    b_e[e](ii) = p[mom_index] - p_0[mom_index];

                    //update maximum and minimum velocities
                    if (p[mom_index] > p_max[mom_index]){
                        p_max[mom_index] = p[mom_index];
                    } else if (p[mom_index] < p_min[mom_index]){
                        p_min[mom_index] = p[mom_index];
                    }
                }

                //solve for mom_pos component of gradient
                sol = A_inv[e]*b_e[e];
                for (int pos = 0; pos<GRID_DIM; pos++){
                    driver->fluid_body->p_x(e, mom_index, pos) = sol(pos);
                }

                //calculate minimum magnitude difference in velocity components
                min_dif[mom_index] = std::abs(p_0[mom_index] - p_min[mom_index]);
                if (std::abs(p_0[mom_index] - p_max[mom_index]) < min_dif[mom_index]) {
                    min_dif[mom_index] = std::abs(p_0[mom_index] - p_max[mom_index]);
                }
            }

            //limit gradient to ensure monotonicity
            for (int mom_index = 0; mom_index<GRID_DIM; mom_index++){
                //calculate maximum momentum change in cell
                tmp_dif = 0;
                for (int pos = 0; pos < GRID_DIM; pos++){
                    tmp_dif += v_e[e]*std::abs(driver->fluid_body->p_x(e, mom_index, pos)/2.0);
                }

                if (tmp_dif > min_dif[mom_index]){
                    //if maximum velocity change in cell is larger than maximum difference b/w neighbors
                    //need to scale gradient
                    for (int pos=0; pos<GRID_DIM; pos++){
                        driver->fluid_body->p_x(e, mom_index, pos) *= min_dif[mom_index]/tmp_dif;
                    }
                }
            }

            //transform gradient from true momentum space to effective momentum space
            driver->fluid_body->p_x[e] *= driver->fluid_body->n_e(e);
        }
    }
    return;
}

void FVMLinear1DNonUniform::constructDensityField(Job* job, FiniteVolumeDriver* driver){
    if (driver->ORDER == 1){
        driver->fluid_body->rho_x.setZero(); //let density be constant within an element
    } if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMLinear1DNonUniform does not currently implement higher ORDER momentum reconstruction."
                      << std::endl;
        }

        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        double rho_0, rho, rho_max, rho_min, min_dif;
        double tmp_dif;
        for (int e = 0; e<element_count; e++){
            //least squares fit of u_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            rho_0 = driver->fluid_body->rho(e) / driver->fluid_body->n_e(e);
            rho_max = rho_0;
            rho_min = rho_0;

            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                //fill in b vector
                rho = driver->fluid_body->rho(element_neighbors[e][ii]) / driver->fluid_body->n_e(element_neighbors[e][ii]);
                b_e[e](ii) = rho - rho_0;

                //update maximum and minimum velocities
                if (rho > rho_max){
                    rho_max = rho;
                } else if (rho < rho_min){
                    rho_min = rho;
                }
            }

            //solve for components of gradient
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                driver->fluid_body->rho_x(e, pos) = sol(pos);
            }

            //calculate minimum magnitude difference in velocity components
            min_dif = std::abs(rho_0 - rho_min);
            if (std::abs(rho_0 - rho_max) < min_dif) {
                min_dif = std::abs(rho_0 - rho_max);
            }

            //limit gradient to ensure monotonicity
            //calculate maximum velocity change in cell
            tmp_dif = 0;
            for (int pos = 0; pos < GRID_DIM; pos++){
                tmp_dif += v_e[e]*std::abs(driver->fluid_body->rho_x(e, pos)/2.0);
            }

            if (tmp_dif > min_dif){
                //if maximum velocity change in cell is larger than maximum difference b/w neighbors
                //need to scale gradient
                for (int pos=0; pos<GRID_DIM; pos++){
                    driver->fluid_body->rho_x(e, pos) *= min_dif/tmp_dif;
                }
            }

            //transform gradient from true density to effective density
            driver->fluid_body->rho_x[e] *=  driver->fluid_body->n_e(e);
        }
    }
    return;
}

void FVMLinear1DNonUniform::constructEnergyField(Job* job, FiniteVolumeDriver* driver){
    if (driver->ORDER == 1){
        driver->fluid_body->rhoE_x.setZero(); //let energy be constant within an element
    } if (driver->ORDER >= 2){
        if (driver->ORDER > 2) {
            std::cout << "ERROR: FVMLinear1DNonUniform does not currently implement higher ORDER momentum reconstruction."
                      << std::endl;
        }

        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        double rhoE_0, rhoE, rhoE_max, rhoE_min, min_dif;
        double tmp_dif;
        for (int e = 0; e<element_count; e++){
            //least squares fit of u_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            rhoE_0 = driver->fluid_body->rhoE(e) / driver->fluid_body->n_e(e);
            rhoE_max = rhoE_0;
            rhoE_min = rhoE_0;

            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                //fill in b vector
                rhoE = driver->fluid_body->rhoE(element_neighbors[e][ii]) / driver->fluid_body->n_e(element_neighbors[e][ii]);
                b_e[e](ii) = rhoE - rhoE_0;

                //update maximum and minimum velocities
                if (rhoE > rhoE_max){
                    rhoE_max = rhoE;
                } else if (rhoE < rhoE_min){
                    rhoE_min = rhoE;
                }
            }

            //solve for components of gradient
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                driver->fluid_body->rhoE_x(e, pos) = sol(pos);
            }

            //calculate minimum magnitude difference in velocity components
            min_dif = std::abs(rhoE_0 - rhoE_min);
            if (std::abs(rhoE_0 - rhoE_max) < min_dif) {
                min_dif = std::abs(rhoE_0 - rhoE_max);
            }

            //limit gradient to ensure monotonicity
            //calculate maximum velocity change in cell
            tmp_dif = 0;
            for (int pos = 0; pos < GRID_DIM; pos++){
                tmp_dif += v_e[e]*std::abs(driver->fluid_body->rhoE_x(e, pos)/2.0);
            }

            if (tmp_dif > min_dif){
                //if maximum velocity change in cell is larger than maximum difference b/w neighbors
                //need to scale gradient
                for (int pos=0; pos<GRID_DIM; pos++){
                    driver->fluid_body->rhoE_x(e, pos) *= min_dif/tmp_dif;
                }
            }

            //transform true energy gradient to effective energy gradient
            driver->fluid_body->rhoE_x[e] *= driver->fluid_body->n_e(e);
        }
    }
    return;
}


void FVMLinear1DNonUniform::constructPorosityField(Job* job, FiniteVolumeDriver* driver){
    //use porosity field gradients to adjust field gradients
    if (USE_LOCAL_POROSITY_CORRECTION) {
        Eigen::VectorXd gradn_star = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x, u, p, tmp_gradn_star;
        KinematicVector rho_x, rhoE_x;
        KinematicTensor p_x;
        double rho, rhoE, M, c;
        double n_0, n, n_max, n_min, min_dif;
        double tmp_dif;
        double dn_ds, P;
        KinematicVector u_s = KinematicVector(job->JOB_TYPE);

        //loop over elements
        for (int e = 0; e < element_count; e++) {
            //estimate gradn_star (reconstructed porosity gradient)
            if (driver->ORDER == 1) {
                //gradients are initially zero
                gradn_star.setZero();
            } else if (driver->ORDER >= 2) {
                if (driver->ORDER > 2) {
                    std::cout << "ERROR: FVMLinear1DNonUniform does not currently implement higher ORDER porosity reconstruction." << std::endl;
                }

                //least squares fit of n to neighbors of element e
                x_0 = getElementCentroid(job, e);
                n_0 = driver->fluid_body->n_e(e);
                n_max = n_0;
                n_min = n_0;

                //create system of equations
                for (int ii = 0; ii < element_neighbors[e].size(); ii++) {
                    //fill in b vector
                    n = driver->fluid_body->n_e(element_neighbors[e][ii]);
                    b_e[e](ii) = n - n_0;

                    //update maximum and minimum velocities
                    if (n > n_max) {
                        n_max = n;
                    } else if (n < n_min) {
                        n_min = n;
                    }
                }

                //solve for components of gradient
                gradn_star = A_inv[e] * b_e[e];

                //calculate minimum magnitude difference in velocity components
                min_dif = std::abs(n_0 - n_min);
                if (std::abs(n_0 - n_max) < min_dif) {
                    min_dif = std::abs(n_0 - n_max);
                }

                //limit gradient to ensure monotonicity
                //calculate maximum velocity change in cell
                tmp_dif = 0;
                for (int pos = 0; pos < GRID_DIM; pos++) {
                    tmp_dif += v_e[e] * std::abs(gradn_star(pos) / 2.0);
                }

                if (tmp_dif > min_dif) {
                    //need to scale gradient
                    for (int pos = 0; pos < GRID_DIM; pos++) {
                        gradn_star(pos) *= min_dif / tmp_dif;
                    }
                }
            }

            //update gradient estimates based on difference between n_e_x and gradn_star
            //pg 37-38 nb #7
            //get cell velocity
            p = driver->fluid_body->p(e);
            rho = driver->fluid_body->rho(e);
            rhoE = driver->fluid_body->rhoE(e);
            u = p/rho;
            n = driver->fluid_body->n_e(e);

            //get average solid velocity in cell
            u_s.setZero();
            for (int q=0; q<qpe; q++){
                u_s += v_sq[e*qpe+q] * w_q(e*qpe+q) / getElementVolume(e);
            }

            //convert p, rho, rhoE to solid reference frame
            p -= rho*u_s;
            rhoE += -p.dot(u_s) + 0.5*rho*u_s.dot(u_s);
            u -= u_s;

            //estimate Mach number
            c = driver->fluid_material->getSpeedOfSound(job, driver, rho, p, rhoE, n);
            M = u.norm()/c;

            if (u.norm() > 1e-10 && n > 1e-10 && M < (1.0 - delta)) {
                tmp_gradn_star = KinematicVector(job->JOB_TYPE);
                for (int pos = 0; pos < GRID_DIM; pos++){
                    tmp_gradn_star[pos] = gradn_star(pos);
                }

                //estimate dn/ds correction (amount not 'seen' by standard reconstruction operation)
                dn_ds = (driver->fluid_body->n_e_x[e] - tmp_gradn_star).dot(u) / u.norm();

                //calculate gradients of properties in solid frame
                rho_x = rho/n * dn_ds * u/u.norm() * (M*M)/(1.0 - M*M);
                p_x = -p.tensor(u)/(u.norm()*n) * dn_ds;

                if (driver->TYPE == driver->THERMAL) {
                    //get pressure
                    P = driver->fluid_material->getPressure(job, driver, rho, p, rhoE, n);
                    rhoE_x = -(rho*c*c - rhoE - P) / n * dn_ds * u / u.norm() * (M * M) / (1.0 - M * M);
                }

                //convert gradients back to spatial frame and update gradients
                driver->fluid_body->rho_x[e] += rho_x;
                driver->fluid_body->p_x[e] += p_x + u_s.tensor(rho_x);

                if (driver->TYPE == driver->THERMAL) {
                    driver->fluid_body->rhoE_x[e] += rhoE_x + p_x.transpose()*u_s + 0.5*u_s.dot(u_s)*rho_x;
                }
            } else {
                //flow too slow, too low porosity, or trans/supersonic
                //do nothing
            }
        }
    } else {
        //do nothing
    }

    return;
}
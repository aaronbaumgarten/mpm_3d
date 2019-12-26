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

void FVMCartesian::init(Job* job, FiniteVolumeDriver* driver){
    //assign grid dimension from job type
    if (job->JOB_TYPE == job->JOB_1D){
        GRID_DIM = 1;
    } else if (job->JOB_TYPE == job->JOB_2D){
        GRID_DIM = 2;
    } else if (job->JOB_TYPE == job->JOB_3D){
        GRID_DIM = 3;
    } else if (job->JOB_TYPE == job->JOB_2D_OOP){
        GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
    } else if (job->JOB_TYPE == job->JOB_AXISYM){
        GRID_DIM = 2; //this is important, job->DIM =/= job->grid->GRID_DIM
    } else {
        std::cerr << "Job doesn't have defined type for input " << job->JOB_TYPE << "." << std::endl;
    }

    //check size of properties passed to driver object
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

    if (job->JOB_TYPE > 3){
        std::cerr << "WARNING: Cannot use FVMCartesian with JOB_TYPE of " << job->JOB_TYPE << "!" << std::endl;
    }

    //initialize tags and values vectors;
    bc_tags = std::vector<int>(2*GRID_DIM);
    bc_values = std::vector<KinematicVector>(2*GRID_DIM);
    for (int i=0; i<2*GRID_DIM; i++){
        //set kinematic vector type from MPM
        bc_values[i] = KinematicVector(job->JOB_TYPE);
    }

    if (int_props.size() == GRID_DIM){
        //no bc tags given
        //set to zero dirichlet
        for (int i=0; i<bc_tags.size(); i++){
            bc_tags[i] = FiniteVolumeGrid::DIRICHLET;
            bc_values[i].setZero();
        }
    } else if (int_props.size() >= 3*GRID_DIM && fp64_props.size() >= GRID_DIM + 2*GRID_DIM*GRID_DIM){
        //bc_tags given
        for (int i=0; i<bc_tags.size(); i++){
            bc_tags[i] = int_props[i+GRID_DIM];
        }

        //check bc_tags for consistency
        for (int i=0; i<GRID_DIM; i++){
            if ((bc_tags[2*i] == FiniteVolumeGrid::PERIODIC
                 || bc_tags[2*i+1] == FiniteVolumeGrid::PERIODIC)
                   && bc_tags[2*i] != bc_tags[2*i+1]){
                std::cerr << "ERROR: Boundary conditions defined on FVMCartesian grid do not match! Exiting." << std::endl;
                exit(0);
            }
        }

        //fill in bc_values
        for (int i=0; i<2*GRID_DIM; i++){
            for (int pos=0; pos<GRID_DIM; pos++){
                bc_values[i][pos] = fp64_props[GRID_DIM + i*GRID_DIM + pos];
            }
        }
    }

    //initialize grid information
    if (GRID_DIM == 1){
        //number of grid elements
        element_count = Nx[0];
        node_count = Nx[0] + 1;
        face_count = Nx[0] + 1;

        //volumes and areas
        element_volumes = hx[0];
        face_areas = std::vector<double>(GRID_DIM);
        face_areas[0] = 1.0;

        //face normals
        face_normals = std::vector<int>(face_count);
        for (int f=0; f<face_count; f++){
            face_normals[f] = 0; //+x
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

        //specify boundary behavior on periodic ends
        if (bc_tags[0] == FiniteVolumeGrid::PERIODIC){
            //remove right face and replace left face
            face_elements[0][0] = element_count-1;
            face_elements[face_count-1][0] = -1;
        }

        //define which boundary is associated with each face
        face_bcs = std::vector<int>(face_count);
        face_bcs[0] = 0;                                    //-x
        face_bcs[face_count-1] = 1;                         //+x
        for (int i=1; i<face_count-1; i++){
            face_bcs[i] = -1;                               //not bc
        }

    } else if (GRID_DIM == 2){
        //number of grid elements
        element_count = Nx[0]*Nx[1];
        node_count = (Nx[0] + 1)*(Nx[1] + 1);
        face_count = 2*Nx[1]*Nx[0] + Nx[0] + Nx[1];

        //volumes and areas
        element_volumes = hx[0]*hx[1];
        face_areas = std::vector<double>(GRID_DIM);
        face_areas[0] = hx[1];
        face_areas[1] = hx[0];

        //face normals
        face_normals = std::vector<int>(face_count);
        for (int i=0; i<Nx[1]*Nx[0]; i++){
            //-x, -y faces of ith element
            face_normals[2*i] = 0; //+x
            face_normals[2*i+1] = 1; //+y
        }
        for (int i=0; i<Nx[1]; i++){
            //+x faces of domain
            face_normals[2*Nx[1]*Nx[0] + i] = 0; //+x
        }
        for (int i=0; i<Nx[0]; i++){
            //+y faces of domain
            face_normals[2*Nx[1]*Nx[0] + Nx[1] + i] = 1; //+y
        }

        //define default face to element definitions.
        face_elements = std::vector<std::array<int,2>>(face_count);
        face_bcs = std::vector<int>(face_count);
        std::vector<int> ijk = std::vector<int>(2);
        std::vector<int> tmp = ijk;
        int e_iminus, e_jminus;
        for (int e=0; e<element_count; e++){
            ijk = e_to_ijk(e);
            //assign element
            face_elements[2*e][1] = e;
            face_elements[2*e+1][1] = e;

            //get element in the -x direction
            tmp = ijk;
            tmp[0] -= 1;
            e_iminus = ijk_to_e(tmp);
            face_elements[2*e][0] = e_iminus;

            //get element in the -y direction
            tmp[0] = ijk[0]; tmp[1] -= 1;
            e_jminus = ijk_to_e(tmp);
            face_elements[2*e+1][0] = e_jminus;

            //check bc_tags
            tmp = ijk;
            if (ijk[0] == 0){
                face_bcs[2*e] = 0;                                  //-x
                if (bc_tags[0] == FiniteVolumeGrid::PERIODIC) {
                    tmp[0] = Nx[0] - 1;
                    face_elements[2 * e][0] = ijk_to_e(tmp);
                }
            } else {
                face_bcs[2*e] = -1;
            }

            tmp = ijk;
            if (ijk[1] == 0){
                face_bcs[2*e+1] = 2;                                //-y
                if (bc_tags[2] == FiniteVolumeGrid::PERIODIC){
                    tmp[1] = Nx[1] - 1;
                    face_elements[2 * e + 1][0] = ijk_to_e(tmp);
                }
            } else {
                face_bcs[2*e+1] = -1;
            }
        }
        for (int i=0; i<Nx[1]; i++){
            //+x faces of domain
            face_bcs[2*Nx[1]*Nx[0] + i] = 1;                    //+x
            if (bc_tags[1] == FiniteVolumeGrid::PERIODIC){
                face_elements[2*Nx[1]*Nx[0] + i][0] = -1;
                face_elements[2*Nx[1]*Nx[0] + i][1] = -1;
            } else {
                ijk[0] = Nx[0] - 1;
                ijk[1] = i;
                face_elements[2 * Nx[1] * Nx[0] + i][0] = ijk_to_e(ijk); //+x
                face_elements[2 * Nx[1] * Nx[0] + i][1] = -1;
            }
        }
        for (int i=0; i<Nx[0]; i++){
            //+y faces of domain
            face_bcs[2*Nx[1]*Nx[0] + Nx[1] + i] = 3;                    //+y
            if (bc_tags[3] == FiniteVolumeGrid::PERIODIC){
                face_elements[2*Nx[1]*Nx[0] + Nx[1] + i][0] = -1;
                face_elements[2*Nx[1]*Nx[0] + Nx[1] + i][1] = -1;
            } else {
                ijk[0] = i;
                ijk[1] = Nx[1] - 1;
                face_elements[2 * Nx[1] * Nx[0] + Nx[1] + i][0] = ijk_to_e(ijk); //+y
                face_elements[2 * Nx[1] * Nx[0] + Nx[1] + i][1] = -1; //+y
            }
        }
    } else if (GRID_DIM == 3){
        //number of grid elements
        element_count = Nx[0]*Nx[1]*Nx[2];
        node_count = (Nx[0] + 1)*(Nx[1] + 1)*(Nx[2] + 1);
        face_count = 3*Nx[0]*Nx[1]*Nx[2] + Nx[1]*Nx[2] + Nx[0]*Nx[2] + Nx[0]*Nx[1];

        //volumes and areas
        element_volumes = hx[0]*hx[1]*hx[2];
        face_areas = std::vector<double>(GRID_DIM);
        face_areas[0] = hx[1]*hx[2];
        face_areas[1] = hx[0]*hx[2];
        face_areas[2] = hx[0]*hx[1];

        //face normals
        face_normals = std::vector<int>(face_count);
        for (int i=0; i<Nx[2]*Nx[1]*Nx[0]; i++){
            //-x, -y faces of ith element
            face_normals[3*i] = 0; //+x
            face_normals[3*i+1] = 1; //+y
            face_normals[3*i+2] = 2; //+z
        }
        for (int i=0; i<Nx[1]*Nx[2]; i++){
            //+x faces of domain
            face_normals[3*Nx[2]*Nx[1]*Nx[0] + i] = 0; //+x
        }
        for (int i=0; i<Nx[0]*Nx[2]; i++){
            //+y faces of domain
            face_normals[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + i] = 1; //+y
        }
        for (int i=0; i<Nx[0]*Nx[1]; i++){
            //+z faces of domain
            face_normals[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + Nx[0]*Nx[2] + i] = 2; //+z
        }

        //define default face to element definitions.
        face_elements = std::vector<std::array<int,2>>(face_count);
        face_bcs = std::vector<int>(face_count);
        std::vector<int> ijk = std::vector<int>(2);
        std::vector<int> tmp = ijk;
        int e_iminus, e_jminus, e_kminus;
        for (int e=0; e<element_count; e++){
            ijk = e_to_ijk(e);
            face_elements[3*e][1] = e;
            face_elements[3*e+1][1] = e;
            face_elements[3*e+2][1] = e;

            //get element in the -x direction
            tmp = ijk;
            tmp[0] -= 1;
            e_iminus = ijk_to_e(tmp);
            face_elements[3*e][0] = e_iminus;

            //get element in the -y direction
            tmp[0] = ijk[0]; tmp[1] -= 1;
            e_jminus = ijk_to_e(tmp);
            face_elements[3*e+1][0] = e_jminus;

            //get element in the -z direction
            tmp[1] = ijk[1]; tmp[2] -= 1;
            e_kminus = ijk_to_e(tmp);
            face_elements[3*e+2][0] = e_kminus;

            //check bc_tags
            tmp = ijk;
            if (ijk[0] == 0){
                face_bcs[3*e] = 0;                                  //-x
                if (bc_tags[0] == FiniteVolumeGrid::PERIODIC) {
                    tmp[0] = Nx[0] - 1;
                    face_elements[3 * e][0] = ijk_to_e(tmp);
                }
            } else {
                face_bcs[3*e] = -1;
            }

            tmp = ijk;
            if (ijk[1] == 0){
                face_bcs[3*e+1] = 2;                                //-y
                if (bc_tags[2] == FiniteVolumeGrid::PERIODIC){
                    tmp[1] = Nx[1] - 1;
                    face_elements[3 * e + 1][0] = ijk_to_e(tmp);
                }
            } else {
                face_bcs[3*e+1] = -1;
            }

            tmp = ijk;
            if (ijk[1] == 0){
                face_bcs[3*e+2] = 4;                                //-y
                if (bc_tags[2] == FiniteVolumeGrid::PERIODIC){
                    tmp[2] = Nx[2] - 1;
                    face_elements[3 * e + 2][0] = ijk_to_e(tmp);
                }
            } else {
                face_bcs[3*e+2] = -1;
            }
        }
        for (int i=0; i<Nx[1]*Nx[2]; i++){
            //+x faces of domain
            face_bcs[3*Nx[2]*Nx[1]*Nx[0] + i] = 1;                    //+x
            if (bc_tags[1] == FiniteVolumeGrid::PERIODIC){
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + i][0] = -1;
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + i][1] = -1;
            } else {
                tmp[0] = Nx[0] - 1;
                tmp[1] = i%Nx[1];
                tmp[2] = ((i - i%Nx[1])/Nx[1]);
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + i][0] = ijk_to_e(tmp); //+x
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + i][1] = -1;
            }
        }
        for (int i=0; i<Nx[0]*Nx[2]; i++){
            //+y faces of domain
            face_bcs[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + i] = 3;                    //+y
            if (bc_tags[3] == FiniteVolumeGrid::PERIODIC){
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + i][0] = -1;
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + i][1] = -1;
            } else {
                tmp[0] = i%Nx[0];
                tmp[1] = Nx[1] - 1;
                tmp[2] = (i - i%Nx[0])/Nx[0];
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + i][0] = ijk_to_e(tmp); //+y
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + i][1] = -1;
            }
        }
        for (int i=0; i<Nx[0]*Nx[1]; i++){
            //+z faces of domain
            face_bcs[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + Nx[0]*Nx[2] + i] = 5;                    //+z
            if (bc_tags[5] == FiniteVolumeGrid::PERIODIC){
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + Nx[0]*Nx[2] + i][0] = -1;
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + Nx[0]*Nx[2] + i][1] = -1;
            } else {
                tmp[0] = i%Nx[0];
                tmp[1] = (i - i%Nx[0])/Nx[0];
                tmp[2] = Nx[2] - 1;
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + Nx[0]*Nx[2] + i][0] = ijk_to_e(tmp); //+z
                face_elements[3*Nx[2]*Nx[1]*Nx[0] + Nx[1]*Nx[2] + Nx[0]*Nx[2] + i][1] = -1;
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
    std::vector<int> tmp = ijk;
    int tmp_e;
    for (int e=0; e<element_count; e++){
        ijk = e_to_ijk(e);
        //loop over i
        for (int i=-1; i<2; i++){
            tmp[0] = ijk[0] + i;
            //check for periodic BC
            for (int f=0; f<element_faces[e].size(); f++){
                //loop -x
                if (face_bcs[f] == 0 && bc_tags[face_bcs[f]] == PERIODIC){
                    if (ijk[0] == 0 && i==-1) {
                        tmp[0] = Nx[0] - 1;
                    } else if (ijk[0] == Nx[0]-1 && i==1){
                        tmp[0] = 0;
                    }
                }
            }
            if (GRID_DIM > 1){
                //loop over j
                for (int j=-1; j<2; j++){
                    tmp[1] = ijk[1] + j;
                    //check for periodic BC
                    for (int f=0; f<element_faces[e].size(); f++){
                        //loop -y
                        if (face_bcs[f] == 2 && bc_tags[face_bcs[f]] == PERIODIC){
                            if (ijk[1] == 0 && j == -1) {
                                tmp[1] = Nx[1] - 1;
                            } else if (ijk[1] == Nx[1] -1 && j == 1) {
                                tmp[1] = 0;
                            }
                        }
                    }
                    if (GRID_DIM > 2){
                        for (int k=-1; k<2; k++){
                            //loop over k
                            tmp[2] = ijk[2] + k;
                            //check for periodic BC
                            for (int f=0; f<element_faces[e].size(); f++){
                                //loop -z
                                if (face_bcs[f] == 4 && bc_tags[face_bcs[f]] == PERIODIC){
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
                    } else {
                        tmp_e = ijk_to_e(tmp);
                        if (tmp_e >= 0 && tmp_e != e && tmp_e < element_count){
                            element_neighbors[e].push_back(tmp_e);
                        }
                    }
                }
            } else {
                tmp_e = ijk_to_e(tmp);
                if (tmp_e >= 0 && tmp_e != e && tmp_e < element_count){
                    element_neighbors[e].push_back(tmp_e);
                }
            }
        }
    }


    //consistency check!!!
    num_neighbors = 3;
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
        for (int j = 0; j < element_faces[e].size(); j++) {
            //only add BCs for dirichlet conditions
            int f = element_faces[e][j];
            if (face_bcs[f] >= 0 && bc_tags[face_bcs[f]] == FiniteVolumeGrid::DIRICHLET) {
                length_of_A++; //add dirichlet conditions to A
            }
        }

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

        //add dirichlet boundary conditions
        int i = element_neighbors[e].size();
        for (int j = 0; j < element_faces[e].size(); j++) {
            //only add BCs for dirichlet conditions
            int f = element_faces[e][j];
            if (face_bcs[f] >= 0 && bc_tags[face_bcs[f]] == FiniteVolumeGrid::DIRICHLET) {
                x = getFaceCentroid(job, f);
                for (int pos = 0; pos < GRID_DIM; pos++) {
                    tmp_dif = x[pos] - x_0[pos];
                    //if periodic, then x - x_0 may be greater than L/2
                    if (tmp_dif > Lx[pos] / 2.0) {
                        A_e[e](i, pos) = tmp_dif - Lx[pos];
                    } else if (tmp_dif < -Lx[pos] / 2.0) {
                        A_e[e](i, pos) = tmp_dif + Lx[pos];
                    } else {
                        A_e[e](i, pos) = x[pos] - x_0[pos];
                    }
                }
                //increment counter
                i++;
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
            a = AtA.householderQr().solve(e_vec);
            for (int jj=0; jj<GRID_DIM; jj++){
                AtA_inv(jj,ii) = a(jj);
            }
        }
        A_inv[e] = AtA_inv*A_e[e].transpose();
    }

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
    std::cout << "    -x: " << bc_tags[0] << " : ";
    for (int pos=0; pos<GRID_DIM; pos++){
        std::cout << bc_values[0][pos] << " ";
    }
    std::cout << " +x: " << bc_tags[1] << " :";
    for (int pos=0; pos<GRID_DIM; pos++){
        std::cout << bc_values[1][pos] << " ";
    }
    std::cout << std::endl;
    if (GRID_DIM > 1){
        std::cout << "    -y: " << bc_tags[2] << " :";
        for (int pos=0; pos<GRID_DIM; pos++){
            std::cout << bc_values[2][pos] << " ";
        }
        std::cout << " +y: " << bc_tags[3] << " :";
        for (int pos=0; pos<GRID_DIM; pos++){
            std::cout << bc_values[3][pos] << " ";
        }
        std::cout << std::endl;
    }
    if (GRID_DIM > 2){
        std::cout << "-z: " << bc_tags[4] << " :";
        for (int pos=0; pos<GRID_DIM; pos++){
            std::cout << bc_values[4][pos];
        }
        std::cout << " +z: " << bc_tags[5] << " :";
        for (int pos=0; pos<GRID_DIM; pos++){
            std::cout << bc_values[5][pos];
        }
        std::cout << std::endl;
    }

    std::cout << "FiniteVolumeGrid initialized." << std::endl;

    /*
    //error checking
    for (int f=0; f<face_elements.size(); f++){
        std::cout << f << " : " << face_elements[f][0] << " | " << face_elements[f][1] << std::endl;
    }
    std::cout << std::endl;
    for (int e=0; e<element_count; e++){
        std::cout << e << " : ";
        for (int i=0; i<element_neighbors[e].size(); i++){
            std::cout << element_neighbors[e][i] << ", ";
        }
        std::cout << std::endl;
    }

    //check periodic faces for double counting
    for (int f=0; f<face_count; f++){
        if (face_bcs[f] > -1 && bc_tags[face_bcs[f]] == PERIODIC){
            std::cout << f << " : " << face_elements[f][0] << " | " << face_elements[f][1] << std::endl;
            std::cout << "  " << face_elements[f][0] << " : ";
            for (int i=0; i<element_neighbors[face_elements[f][0]].size(); i++){
                std::cout << element_neighbors[face_elements[f][0]][i] << " ";
            }
            std::cout << std::endl << "  " << face_elements[f][1] << " : ";
            for (int i=0; i<element_neighbors[face_elements[f][1]].size(); i++){
                std::cout << element_neighbors[face_elements[f][1]][i] << " ";
            }
            std::cout << std::endl;
        }
    }
     */
}

/*----------------------------------------------------------------------------*/

std::vector<int> FVMCartesian::e_to_ijk(int e){
    std::vector<int> ijk = std::vector<int>(GRID_DIM);
    int tmp = e;
    for (int i=0;i<ijk.size();i++){
        ijk[i] = tmp % Nx[i];
        tmp = tmp/Nx[i];
    }
    return ijk;
}

int FVMCartesian::ijk_to_e(std::vector<int> ijk){
    int tmp = 0;
    for (int i=ijk.size(); i>0; i--){
        if (ijk[i-1] < 0 || ijk[i-1] >= Nx[i-1]){
            return -1;
        }
        tmp = tmp*Nx[i-1] + ijk[i-1];
    }
    return tmp;
}

std::vector<int> FVMCartesian::n_to_ijk(int n){
    std::vector<int> ijk = std::vector<int>(GRID_DIM);
    int tmp = n;
    for (int i=0;i<ijk.size();i++){
        ijk[i] = tmp % (Nx[i]+1);
        tmp = tmp/(Nx[i]+1);
    }
    return ijk;
}

int FVMCartesian::ijk_to_n(std::vector<int> ijk){
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

void FVMCartesian::writeHeader(std::ofstream& file, int SPEC){
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
                file << ijk[pos]*hx[pos] << " ";
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

double FVMCartesian::getElementVolume(int e){
    return element_volumes;
}

int FVMCartesian::getElementTag(int e){
    //get ijk position of element id
    //std::vector<int> ijk = e_to_ijk(e);

    //do nothing.
    return 0;
}

double FVMCartesian::getFaceArea(int f){
    //get normal definition
    int i = face_normals[f];
    return face_areas[i];
}

int FVMCartesian::getFaceTag(int f){
    //return boundary tag of faces
    return face_bcs[f];
}

KinematicVector FVMCartesian::getFaceNormal(Job* job, int f){
    //return oriented normal associated with face
    KinematicVector tmp = KinematicVector(job->JOB_TYPE);
    int i = face_normals[f];
    tmp[i] = 1;
    return tmp;
}

std::vector<int> FVMCartesian::getElementFaces(int e){
    //get faces associated with each element
    return element_faces[e];
}

std::array<int,2> FVMCartesian::getOrientedElementsByFace(int f){
    //get oriented elements associated with given face
    // A -|-> B
    return face_elements[f];
}

std::vector<int> FVMCartesian::getElementNeighbors(int e){
    //get list of elements that share a node with e
    return element_neighbors[e];
}


KinematicVector FVMCartesian::getElementCentroid(Job* job, int e){
    std::vector<int> tmp = e_to_ijk(e);
    KinematicVector x = KinematicVector(job->JOB_TYPE);
    for (int i=0; i<x.size(); i++){
        x[i] = hx[i]*tmp[i] + hx[i]/2.0;
    }
    return x;
}

KinematicVector FVMCartesian::getFaceCentroid(Job* job, int f){
    KinematicVector x = KinematicVector(job->JOB_TYPE);
    if (f < GRID_DIM*element_count){
        //face is in first set, therefore element B defines location
        int e = face_elements[f][1];
        std::vector<int> tmp = e_to_ijk(e);
        for (int i=0; i<x.size(); i++){
            x[i] = hx[i]*tmp[i] + hx[i]/2.0; //centroid of element
        }
        x[face_normals[f]] -= hx[face_normals[f]]/2.0; //move to face
    } else {
        //face is in second set, therefore element A defines location
        int e = face_elements[f][0];
        if (e >= 0) {
            //element exists
            std::vector<int> tmp = e_to_ijk(e);
            for (int i = 0; i < x.size(); i++) {
                x[i] = hx[i] * tmp[i] + hx[i] / 2.0; //centroid of element
            }
            x[face_normals[f]] += hx[face_normals[f]] / 2.0; //move to face
        } else {
            //face is not used on current grid (must be periodic)
            x.setZero();
        }
    }
    return x;
}

/*----------------------------------------------------------------------------*/
//functions for solving system of equations
void FVMCartesian::generateMappings(Job* job, FiniteVolumeDriver* driver){
    //for now do nothing
    std::cout << "WARNING: generateMappings function not implemented in FVMCartesian." << std::endl;
    return;
}

void FVMCartesian::constructMomentumField(Job* job, FiniteVolumeDriver* driver){
    if (driver->order == 1){
        driver->fluid_body->p_x.setZero(); //let momentum be constant within an element
    } else if (driver->order >= 2){
        if (driver->order > 2) {
            std::cout << "ERROR: FVMCartesian does not currently implement higher order momentum reconstruction."
                      << std::endl;
        }
        //initialize matrix A and vector B
        //Eigen::MatrixXd A = Eigen::MatrixXd(num_neighbors, GRID_DIM);
        //Eigen::VectorXd b = Eigen::VectorXd(num_neighbors);
        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        KinematicVector p_0, p, p_max, p_min, min_dif;
        min_dif = KinematicVector(job->JOB_TYPE);
        double tmp_dif, rho_0;
        for (int e = 0; e<element_count; e++){
            //least squares fit of u_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            p_0 = driver->fluid_body->p[e];
            rho_0 = driver->fluid_body->rho(e);
            p_max = p_0; p_min = p_0;
            min_dif.setZero();
            for (int mom_index = 0; mom_index<GRID_DIM; mom_index++){
                //create system of equations
                for (int ii=0; ii<element_neighbors[e].size(); ii++){
                    //x = getElementCentroid(job, element_neighbors[e][ii]);
                    /*
                    for (int pos = 0; pos<GRID_DIM; pos++) {
                        tmp_dif = x[pos] - x_0[pos];
                        //if periodic, then x - x_0 may be greater than L/2
                        if (tmp_dif > Lx[pos]/2.0){
                            A(ii, pos) = tmp_dif - Lx[pos];
                        } else if (tmp_dif < -Lx[pos]/2.0){
                            A(ii, pos) = tmp_dif + Lx[pos];
                        } else {
                            A(ii, pos) = x[pos] - x_0[pos];
                        }
                    }
                     */

                    //just fill in b vector
                    p = driver->fluid_body->p[element_neighbors[e][ii]];
                    b_e[e](ii) = p[mom_index] - p_0[mom_index];

                    //update maximum and minimum velocities
                    if (p[mom_index] > p_max[mom_index]){
                        p_max[mom_index] = p[mom_index];
                    } else if (p[mom_index] < p_min[mom_index]){
                        p_min[mom_index] = p[mom_index];
                    }
                }

                //add boundary conditions where applicable (total number of eq'ns still less than rows of A)
                int i=element_neighbors[e].size();
                for (int j = 0; j<element_faces[e].size(); j++){
                    //only add BCs for dirichlet conditions
                    int f = element_faces[e][j];
                    if (face_bcs[f] >= 0 && bc_tags[face_bcs[f]] == FiniteVolumeGrid::DIRICHLET) {
                        /*
                        x = getFaceCentroid(job, f);
                        p = rho_0 * bc_values[face_bcs[f]];
                        for (int pos = 0; pos < GRID_DIM; pos++) {
                            tmp_dif = x[pos] - x_0[pos];
                            //if periodic, then x - x_0 may be greater than L/2
                            if (tmp_dif > Lx[pos] / 2.0) {
                                A(i, pos) = tmp_dif - Lx[pos];
                            } else if (tmp_dif < -Lx[pos] / 2.0) {
                                A(i, pos) = tmp_dif + Lx[pos];
                            } else {
                                A(i, pos) = x[pos] - x_0[pos];
                            }
                        }
                         */
                        //fill in b vector
                        p = rho_0 * bc_values[face_bcs[f]];
                        b_e[e](i) = p[mom_index] - p_0[mom_index];

                        //update maximum and minimum velocities
                        if (p[mom_index] > p_max[mom_index]) {
                            p_max[mom_index] = p[mom_index];
                        } else if (p[mom_index] < p_min[mom_index]) {
                            p_min[mom_index] = p[mom_index];
                        }
                        //increment i
                        i++;
                    }
                }

                /*
                //zero remainder of system of eq'ns
                for (int ii=i; ii<num_neighbors; ii++){
                    for (int pos=0; pos<GRID_DIM; pos++){
                        A(ii,pos) = 0;
                        b(ii) = 0;
                    }
                }
                 */

                //solve for u_pos component of gradient
                //sol = A_e[e].householderQr().solve(b_e[e]);
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
                    tmp_dif += hx[pos]*std::abs(driver->fluid_body->p_x(e, mom_index, pos)/2.0);
                }

                if (tmp_dif > min_dif[mom_index]){
                    //if maximum velocity change in cell is larger than maximum difference b/w neighbors
                    //need to scale gradient
                    for (int pos=0; pos<GRID_DIM; pos++){
                        driver->fluid_body->p_x(e, mom_index, pos) *= min_dif[mom_index]/tmp_dif;
                    }
                }
            }
        }
    }
    return;
}

void FVMCartesian::constructDensityField(Job* job, FiniteVolumeDriver* driver){
    if (driver->order == 1){
        driver->fluid_body->rho_x.setZero(); //let density be constant within an element
    } if (driver->order >= 2){
        if (driver->order > 2) {
            std::cout << "ERROR: FVMCartesian does not currently implement higher order momentum reconstruction."
                      << std::endl;
        }
        //initialize matrix A and vector B
        //Eigen::MatrixXd A = Eigen::MatrixXd(num_neighbors, GRID_DIM);
        //Eigen::VectorXd b = Eigen::VectorXd(num_neighbors);
        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x;
        double rho_0, rho, rho_max, rho_min, min_dif;
        double tmp_dif;
        for (int e = 0; e<element_count; e++){
            //least squares fit of u_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            rho_0 = driver->fluid_body->rho(e);
            rho_max = rho_0;
            rho_min = rho_0;

            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                /*
                x = getElementCentroid(job, element_neighbors[e][ii]);
                for (int pos = 0; pos<GRID_DIM; pos++) {
                    tmp_dif = x[pos] - x_0[pos];
                    //if periodic, then x - x_0 may be greater than L/2
                    if (tmp_dif > Lx[pos]/2.0){
                        A(ii, pos) = tmp_dif - Lx[pos];
                    } else if (tmp_dif < -Lx[pos]/2.0){
                        A(ii, pos) = tmp_dif + Lx[pos];
                    } else {
                        A(ii, pos) = x[pos] - x_0[pos];
                    }
                }
                 */

                //fill in b vector
                rho = driver->fluid_body->rho(element_neighbors[e][ii]);
                b_e[e](ii) = rho - rho_0;

                //update maximum and minimum velocities
                if (rho > rho_max){
                    rho_max = rho;
                } else if (rho < rho_min){
                    rho_min = rho;
                }
            }

            //add boundary conditions where applicable (total number of eq'ns still less than rows of A)
            int i=element_neighbors[e].size();
            for (int j = 0; j<element_faces[e].size(); j++){
                //only add BCs for dirichlet conditions
                int f = element_faces[e][j];
                if (face_bcs[f] >= 0 && bc_tags[face_bcs[f]] == FiniteVolumeGrid::DIRICHLET) {
                    //x = getFaceCentroid(job, f);
                    //rho = rho_0;
                    /*
                    for (int pos = 0; pos < GRID_DIM; pos++) {
                        tmp_dif = x[pos] - x_0[pos];
                        //if periodic, then x - x_0 may be greater than L/2
                        if (tmp_dif > Lx[pos] / 2.0) {
                            A(i, pos) = tmp_dif - Lx[pos];
                        } else if (tmp_dif < -Lx[pos] / 2.0) {
                            A(i, pos) = tmp_dif + Lx[pos];
                        } else {
                            A(i, pos) = x[pos] - x_0[pos];
                        }
                    }
                     */
                    b_e[e](i) = 0;

                    //increment i
                    i++;
                }
            }

            /*
            //zero remainder of system of eq'ns
            for (int ii=i; ii<num_neighbors; ii++){
                for (int pos=0; pos<GRID_DIM; pos++){
                    A(ii,pos) = 0;
                    b(ii) = 0;
                }
            }
             */

            //solve for u_pos component of gradient
            //sol = A_e[e].householderQr().solve(b_e[e]);
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
                tmp_dif += hx[pos]*std::abs(driver->fluid_body->rho_x(e, pos)/2.0);
            }

            if (tmp_dif > min_dif){
                //if maximum velocity change in cell is larger than maximum difference b/w neighbors
                //need to scale gradient
                for (int pos=0; pos<GRID_DIM; pos++){
                    driver->fluid_body->rho_x(e, pos) *= min_dif/tmp_dif;
                }
            }
            /*
            std::cout << e << std::endl;
            std::cout << A << std::endl;
             */
        }
    }
    return;
}

KinematicTensorArray FVMCartesian::getVelocityGradients(Job* job, FiniteVolumeDriver* driver){
    //reconstruct velocity field
    KinematicTensorArray u_x = KinematicTensorArray(element_count, job->JOB_TYPE);
    KinematicVector u = KinematicVector(job->JOB_TYPE);
    KinematicVector p = KinematicVector(job->JOB_TYPE);
    double rho, rho_0, tmp_dif;

    //initialize matrix A and vector B
    //Eigen::MatrixXd A = Eigen::MatrixXd(num_neighbors, GRID_DIM);
    //Eigen::VectorXd b = Eigen::VectorXd(num_neighbors);
    Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
    KinematicVector x_0, x;
    KinematicVector u_0;

    for (int e=0; e<element_count; e++){
        p = driver->fluid_body->p[e];
        rho_0 = driver->fluid_body->rho(e);
        u_0 = p/rho_0;

        //least squares fit of u_x to neighbors of element e
        x_0 = getElementCentroid(job, e);
        for (int dir = 0; dir<GRID_DIM; dir++){
            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                /*
                x = getElementCentroid(job, element_neighbors[e][ii]);
                p = driver->fluid_body->p[element_neighbors[e][ii]];
                rho = driver->fluid_body->rho(element_neighbors[e][ii]);
                for (int pos = 0; pos<GRID_DIM; pos++) {
                    tmp_dif = x[pos] - x_0[pos];
                    //if periodic, then x - x_0 may be greater than L/2
                    if (tmp_dif > Lx[pos]/2.0){
                        A(ii, pos) = tmp_dif - Lx[pos];
                    } else if (tmp_dif < -Lx[pos]/2.0){
                        A(ii, pos) = tmp_dif + Lx[pos];
                    } else {
                        A(ii, pos) = x[pos] - x_0[pos];
                    }
                }
                 */
                rho = driver->fluid_body->rho(element_neighbors[e][ii]);
                p = driver->fluid_body->p[element_neighbors[e][ii]];
                b_e[e](ii) = p[dir]/rho - u_0[dir];
            }

            //add boundary conditions where applicable (total number of eq'ns still less than rows of A)
            int i=element_neighbors[e].size();
            for (int j = 0; j<element_faces[e].size(); j++){
                //only add BCs for dirichlet conditions
                int f = element_faces[e][j];
                if (face_bcs[f] >= 0 && bc_tags[face_bcs[f]] == FiniteVolumeGrid::DIRICHLET) {
                    /*
                    x = getFaceCentroid(job, f);
                    for (int pos = 0; pos < GRID_DIM; pos++) {
                        tmp_dif = x[pos] - x_0[pos];
                        //if periodic, then x - x_0 may be greater than L/2
                        if (tmp_dif > Lx[pos] / 2.0) {
                            A(i, pos) = tmp_dif - Lx[pos];
                        } else if (tmp_dif < -Lx[pos] / 2.0) {
                            A(i, pos) = tmp_dif + Lx[pos];
                        } else {
                            A(i, pos) = x[pos] - x_0[pos];
                        }
                    }
                     */
                    u = bc_values[face_bcs[f]];
                    b_e[e](i) = u[dir] - u_0[dir];

                    //increment i
                    i++;
                }
            }

            /*
            //zero remainder of system of eq'ns
            for (int ii=i; ii<num_neighbors; ii++){
                for (int pos=0; pos<GRID_DIM; pos++){
                    A(ii,pos) = 0;
                    b(ii) = 0;
                }
            }
             */

            //solve for u_pos component of gradient
            //sol = A_e[e].householderQr().solve(b_e[e]);
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                u_x(e, dir, pos) = sol(pos);
            }
        }

    }
    return u_x;
}

/*----------------------------------------------------------------------------*/
//functions to compute element-wise fluxes of field variables using reconstructed velocity field
Eigen::VectorXd FVMCartesian::calculateElementFluxIntegrals(Job* job, FiniteVolumeDriver* driver, Eigen::VectorXd& phi){
    if (phi.rows() != element_count){
        std::cerr << "ERROR! Length of vector passed to calculateElementFluxIntegrals is wrong size! ";
        std::cerr << phi.rows() << " != " << element_count << std::endl;
        exit(0);
    }

    //flux rates
    Eigen::VectorXd result = Eigen::VectorXd(element_count);
    result.setZero();


    double v_plus, v_minus; //value
    KinematicVector u_plus, u_minus, normal, u_bar; //velocity
    int e_plus, e_minus;    //elements
    double rho_plus, rho_minus;
    double flux, area;
    double a_1, a_2, a_6, c, rho_bar, phi_bar, theta_bar;
    double lambda_1, lambda_2;

    //number of quadrature points per face
    int num_quad = 1;
    if (GRID_DIM == 2){
        num_quad = 2;
    } else if (GRID_DIM == 3){
        num_quad = 4;
    }

    //vector of quad point positions relative to face center
    std::array<std::array<double,2>,4> quad_points;
    quad_points[0][0] = -1.0/std::sqrt(3);  quad_points[0][1] = -1.0/std::sqrt(3);
    quad_points[1][0] = 1.0/std::sqrt(3);   quad_points[1][1] = -1.0/std::sqrt(3);
    quad_points[2][0] = -1.0/std::sqrt(3);  quad_points[2][1] = 1.0/std::sqrt(3);
    quad_points[3][0] = 1.0/std::sqrt(3);   quad_points[3][1] = 1.0/std::sqrt(3);

    //different flux calculation for different order approximations
    if (driver->order == 1){
        //rho, u, value constant within element
        //loop over faces
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);
        KinematicVector x;

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (face_bcs[f] == -1 || (bc_tags[face_bcs[f]] == PERIODIC && face_elements[f][0] > -1)){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];
                v_minus = phi(e_minus);
                v_plus  = phi(e_plus);
                rho_minus = driver->fluid_body->rho(e_minus);
                rho_plus = driver->fluid_body->rho(e_plus);
                u_minus = driver->fluid_body->p[e_minus]/rho_minus; //u = p/rho
                u_plus  = driver->fluid_body->p[e_plus]/rho_plus;

                //approximate Roe advective rate
                u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                rho_bar = std::sqrt(rho_plus*rho_minus);
                phi_bar = (v_plus*std::sqrt(rho_minus) + v_minus*std::sqrt(rho_plus))/(std::sqrt(rho_minus) + std::sqrt(rho_plus));
                theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                c = driver->fluid_material->getSpeedOfSound(job, driver, x, rho_bar, theta_bar);

                //roe eigenvalues
                lambda_1 = std::abs(u_bar.dot(normal) - c);
                lambda_2 = std::abs(u_bar.dot(normal) + c);
                if (lambda_1 < delta*c){
                    lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                }
                if (lambda_2 < delta*c){
                    lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                }

                //calculate Roe eigenvector coefficients
                a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                a_2 = (rho_plus - rho_minus) - a_1;
                a_6 = v_plus - v_minus - (rho_plus - rho_minus)*phi_bar/rho_bar;

                //flux in n direction
                flux = area*0.5*(v_plus*u_plus.dot(normal) + v_minus*u_minus.dot(normal)
                                                - a_1*lambda_1*phi_bar/rho_bar
                                                - a_2*lambda_2*phi_bar/rho_bar
                                                - a_6*std::abs(u_bar.dot(normal)));

                //add flux to element integrals
                result(e_minus) -= flux;
                result(e_plus) += flux;
            } else if (bc_tags[face_bcs[f]] == DIRICHLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and assigned velocity
                    v_minus = phi(e_minus);
                    flux = v_minus*area*bc_values[face_bcs[f]].dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and assigned velocity
                    v_plus = phi(e_plus);
                    flux = v_plus*area*bc_values[face_bcs[f]].dot(normal);
                    result(e_plus) += flux;
                }
            } else if (bc_tags[face_bcs[f]] == NEUMANN){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and cell velocity
                    v_minus = phi(e_minus);
                    u_minus = driver->fluid_body->p[e_minus]/driver->fluid_body->rho(e_minus);
                    flux = v_minus*area*u_minus.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and cell velocity
                    v_plus = phi(e_plus);
                    u_plus = driver->fluid_body->p[e_plus]/driver->fluid_body->rho(e_plus);
                    flux = v_plus*area*u_plus.dot(normal);
                    result(e_plus) += flux;
                }
            }
        }
    } else if (driver->order >= 2){
        if (driver->order > 2) {
            std::cout << "ERROR: FVMCartesian does not currently implement higher order flux reconstruction."
                      << std::endl;
        }

        //first reconstruct input field
        KinematicVectorArray phi_x = KinematicVectorArray(element_count, job->JOB_TYPE);

        //initialize matrix A and vector B
        //Eigen::MatrixXd A = Eigen::MatrixXd(num_neighbors, GRID_DIM);
        //Eigen::VectorXd b = Eigen::VectorXd(num_neighbors);
        Eigen::VectorXd sol = Eigen::VectorXd(GRID_DIM);
        KinematicVector x_0, x, x_face, x_quad;
        x_quad = KinematicVector(job->JOB_TYPE);
        double value_0, value, value_max, value_min, min_dif;
        double tmp_dif;
        for (int e = 0; e<element_count; e++){
            //least squares fit of u_x to neighbors of element e
            x_0 = getElementCentroid(job, e);
            value_0 = phi(e);
            min_dif = 0;

            //create system of equations
            for (int ii=0; ii<element_neighbors[e].size(); ii++){
                /*
                x = getElementCentroid(job, element_neighbors[e][ii]);
                value = phi(element_neighbors[e][ii]);
                for (int pos = 0; pos<GRID_DIM; pos++) {
                    tmp_dif = x[pos] - x_0[pos];
                    //if periodic, then x - x_0 may be greater than L/2
                    if (tmp_dif > Lx[pos]/2.0){
                        A(ii, pos) = tmp_dif - Lx[pos];
                    } else if (tmp_dif < -Lx[pos]/2.0){
                        A(ii, pos) = tmp_dif + Lx[pos];
                    } else {
                        A(ii, pos) = x[pos] - x_0[pos];
                    }
                }
                 */

                value = phi(element_neighbors[e][ii]);
                b_e[e](ii) = value - value_0;

                //update maximum and minimum velocities
                if (value > value_max){
                    value_max = value;
                } else if (value < value_min){
                    value_min = value;
                }
            }

            //zero remainder of system of eq'ns
            for (int ii=element_neighbors[e].size(); ii<b_e[e].rows(); ii++){
                for (int pos=0; pos<GRID_DIM; pos++){
                    //A(ii,pos) = 0;
                    b_e[e](ii) = 0;
                }
            }

            //solve for u_pos component of gradient
            //sol = A_e[e].householderQr().solve(b_e[e]);
            sol = A_inv[e]*b_e[e];
            for (int pos = 0; pos<GRID_DIM; pos++){
                phi_x(e, pos) = sol(pos);
            }

            //calculate minimum magnitude difference in velocity components
            min_dif = std::abs(value_0 - value_min);
            if (std::abs(value_0 - value_max) < min_dif) {
                min_dif = std::abs(value_0 - value_max);
            }


            //limit gradient to ensure monotonicity
            //calculate maximum velocity change in cell
            tmp_dif = 0;
            for (int pos = 0; pos < GRID_DIM; pos++){
                tmp_dif += hx[pos]*std::abs(phi_x(e, pos)/2.0);
            }

            if (tmp_dif > min_dif){
                //if maximum velocity change in cell is larger than maximum difference b/w neighbors
                //need to scale gradient
                for (int pos=0; pos<GRID_DIM; pos++){
                    phi_x(e, pos) *= min_dif/tmp_dif;
                }
            }
        }

        //loop over faces and use quadrature to reconstruct flux integral
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x_face = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (face_bcs[f] == -1 || (bc_tags[face_bcs[f]] == PERIODIC && face_elements[f][0] > -1)){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    int ii=0;
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = normal[pos] * hx[pos]/2.0;
                            x_quad[pos] = x_face[pos];
                        } else {
                            //move along face
                            x[pos] = hx[pos]/2.0 * quad_points[q][ii];
                            x_quad[pos] = x_face[pos] + x[pos];
                            ii++;
                        }
                    }

                    //calculate A properties
                    v_minus = phi(e_minus) + phi_x(e_minus).dot(x);
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    u_minus = (driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x)/rho_minus;

                    //relative position to centroid of B
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = -normal[pos] * hx[pos]/2.0;
                        }
                    }

                    //calculate B properties
                    v_plus = phi(e_plus) + phi_x(e_plus).dot(x);
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x)/rho_plus;

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = std::sqrt(rho_plus*rho_minus);
                    phi_bar = (v_plus*std::sqrt(rho_minus) + v_minus*std::sqrt(rho_plus))/(std::sqrt(rho_minus) + std::sqrt(rho_plus));
                    theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                    c = driver->fluid_material->getSpeedOfSound(job, driver, x_quad, rho_bar, theta_bar);

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    if (lambda_1 < delta*c){
                        lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                    }
                    if (lambda_2 < delta*c){
                        lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                    }

                    //calculate Roe eigenvector coefficients
                    a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                    a_2 = (rho_plus - rho_minus) - a_1;
                    a_6 = v_plus - v_minus - (rho_plus - rho_minus)*phi_bar/rho_bar;

                    //flux in n direction
                    flux = area/num_quad * 0.5 * (v_plus*u_plus.dot(normal) + v_minus*u_minus.dot(normal)
                                                                            - a_1*lambda_1*phi_bar/rho_bar
                                                                            - a_2*lambda_2*phi_bar/rho_bar
                                                                            - a_6*std::abs(u_bar.dot(normal)));

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else if (bc_tags[face_bcs[f]] == DIRICHLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    int ii=0;
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = normal[pos] * hx[pos]/2.0;
                        } else {
                            //move along face
                            x[pos] = hx[pos]/2.0 * quad_points[q][ii];
                            ii++;
                        }
                    }

                    if (e_minus > -1) {
                        //calculate A properties
                        v_minus = phi(e_minus) + phi_x(e_minus).dot(x);
                        flux = v_minus*area/num_quad*bc_values[face_bcs[f]].dot(normal);
                        result(e_minus) -= flux;
                    }

                    //relative position to centroid of B
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = -normal[pos] * hx[pos]/2.0;
                        }
                    }

                    if (e_plus > -1) {
                        //calculate B properties
                        v_plus = phi(e_plus) + phi_x(e_plus).dot(x);
                        flux = v_plus*area/num_quad*bc_values[face_bcs[f]].dot(normal);
                        result(e_plus) += flux;
                    }
                }
            } else if (bc_tags[face_bcs[f]] == NEUMANN){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    int ii=0;
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = normal[pos] * hx[pos]/2.0;
                        } else {
                            //move along face
                            x[pos] = hx[pos]/2.0 * quad_points[q][ii];
                            ii++;
                        }
                    }

                    if (e_minus > -1) {
                        //calculate A properties
                        v_minus = phi(e_minus) + phi_x(e_minus).dot(x);
                        rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        u_minus = (driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x)/rho_minus;
                        flux = v_minus*area/num_quad*u_minus.dot(normal);
                        result(e_minus) -= flux;
                    }

                    //relative position to centroid of B
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = -normal[pos] * hx[pos]/2.0;
                        }
                    }

                    if (e_plus > -1) {
                        //calculate B properties
                        v_plus = phi(e_plus) + phi_x(e_plus).dot(x);
                        rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x)/rho_plus;
                        flux = v_plus*area/num_quad*u_plus.dot(normal);
                        result(e_plus) += flux;
                    }
                }
            }
        }
    }
    return result;
}

//functions to compute element mass flux
Eigen::VectorXd FVMCartesian::calculateElementMassFluxes(Job* job, FiniteVolumeDriver* driver){
    //flux rates
    Eigen::VectorXd result = Eigen::VectorXd(element_count);
    result.setZero();

    KinematicVector u_plus, u_minus, normal, u_bar; //velocity
    int e_plus, e_minus;    //elements
    double rho_plus, rho_minus;
    double flux, area;
    double lambda_1, lambda_2, a_1, a_2, c, rho_bar, theta_bar;

    //number of quadrature points per face
    int num_quad = 1;
    if (GRID_DIM == 2){
        num_quad = 2;
    } else if (GRID_DIM == 3){
        num_quad = 4;
    }

    //vector of quad point positions relative to face center
    std::array<std::array<double,2>,4> quad_points;
    quad_points[0][0] = -1.0/std::sqrt(3);  quad_points[0][1] = -1.0/std::sqrt(3);
    quad_points[1][0] = 1.0/std::sqrt(3);   quad_points[1][1] = -1.0/std::sqrt(3);
    quad_points[2][0] = -1.0/std::sqrt(3);  quad_points[2][1] = 1.0/std::sqrt(3);
    quad_points[3][0] = 1.0/std::sqrt(3);   quad_points[3][1] = 1.0/std::sqrt(3);

    //different flux calculation for different order approximations
    if (driver->order == 1){
        //rho, u, value constant within element
        //loop over faces
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);
        KinematicVector x;

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (face_bcs[f] == -1 || (bc_tags[face_bcs[f]] == PERIODIC && face_elements[f][0] > -1)){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];
                rho_minus = driver->fluid_body->rho(e_minus);
                rho_plus = driver->fluid_body->rho(e_plus);
                u_minus = driver->fluid_body->p[e_minus]/rho_minus; //u = p/rho
                u_plus  = driver->fluid_body->p[e_plus]/rho_plus;

                //approximate Roe advective rate
                u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                rho_bar = std::sqrt(rho_plus*rho_minus);
                theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                c = driver->fluid_material->getSpeedOfSound(job, driver, x, rho_bar, theta_bar);

                //roe eigenvalues
                lambda_1 = std::abs(u_bar.dot(normal) - c);
                lambda_2 = std::abs(u_bar.dot(normal) + c);
                if (lambda_1 < delta*c){
                    lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                }
                if (lambda_2 < delta*c){
                    lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                }

                //calculate Roe eigenvector coefficients
                a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                a_2 = (rho_plus - rho_minus) - a_1;

                //flux in n direction
                flux = area * 0.5 * (rho_plus*u_plus.dot(normal) + rho_minus*u_minus.dot(normal)
                                                                 - a_1*lambda_1 - a_2*lambda_2);

                //add flux to element integrals
                result(e_minus) -= flux;
                result(e_plus) += flux;
            } else if (bc_tags[face_bcs[f]] == DIRICHLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and assigned velocity
                    rho_minus = driver->fluid_body->rho(e_minus);
                    flux = rho_minus*area*bc_values[face_bcs[f]].dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and assigned velocity
                    rho_plus = driver->fluid_body->rho(e_plus);
                    flux = rho_plus*area*bc_values[face_bcs[f]].dot(normal);
                    result(e_plus) += flux;
                }
            } else if (bc_tags[face_bcs[f]] == NEUMANN){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and cell velocity
                    rho_minus = driver->fluid_body->rho(e_minus);
                    u_minus = driver->fluid_body->p[e_minus]/driver->fluid_body->rho(e_minus);
                    flux = rho_minus*area*u_minus.dot(normal);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and cell velocity
                    rho_plus = driver->fluid_body->rho(e_plus);
                    u_plus = driver->fluid_body->p[e_plus]/driver->fluid_body->rho(e_plus);
                    flux = rho_plus*area*u_plus.dot(normal);
                    result(e_plus) += flux;
                }
            }
        }
    } else if (driver->order >= 2){
        if (driver->order > 2) {
            std::cout << "ERROR: FVMCartesian does not currently implement higher order flux reconstruction."
                      << std::endl;
        }

        KinematicVector x = KinematicVector(job->JOB_TYPE);
        KinematicVector x_face = KinematicVector(job->JOB_TYPE);
        KinematicVector x_quad = KinematicVector(job->JOB_TYPE);

        //loop over faces and use quadrature to reconstruct flux integral
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x_face = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (face_bcs[f] == -1 || (bc_tags[face_bcs[f]] == PERIODIC && face_elements[f][0] > -1)){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    int ii=0;
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = normal[pos] * hx[pos]/2.0;
                            x_quad[pos] = x_face[pos];
                        } else {
                            //move along face
                            x[pos] = hx[pos]/2.0 * quad_points[q][ii];
                            x_quad[pos] = x_face[pos] + x[pos];
                            ii++;
                        }
                    }

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    u_minus = (driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x)/rho_minus;

                    //relative position to centroid of B
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = -normal[pos] * hx[pos]/2.0;
                        }
                    }

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x)/rho_plus;

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = std::sqrt(rho_plus*rho_minus);
                    theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                    c = driver->fluid_material->getSpeedOfSound(job, driver, x_quad, rho_bar, theta_bar);

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    if (lambda_1 < delta*c){
                        lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                    }
                    if (lambda_2 < delta*c){
                        lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                    }

                    //calculate Roe eigenvector coefficients
                    a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                    a_2 = (rho_plus - rho_minus) - a_1;

                    //flux in n direction
                    flux = area/num_quad * 0.5 * (rho_plus*u_plus.dot(normal) + rho_minus*u_minus.dot(normal)
                                         - a_1*lambda_1 - a_2*lambda_2);

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else if (bc_tags[face_bcs[f]] == DIRICHLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    int ii=0;
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = normal[pos] * hx[pos]/2.0;
                        } else {
                            //move along face
                            x[pos] = hx[pos]/2.0 * quad_points[q][ii];
                            ii++;
                        }
                    }

                    if (e_minus > -1) {
                        //calculate A properties
                        rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        flux = rho_minus*area/num_quad*bc_values[face_bcs[f]].dot(normal);
                        result(e_minus) -= flux;
                    }

                    //relative position to centroid of B
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = -normal[pos] * hx[pos]/2.0;
                        }
                    }

                    if (e_plus > -1) {
                        //calculate B properties
                        rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        flux = rho_plus*area/num_quad*bc_values[face_bcs[f]].dot(normal);
                        result(e_plus) += flux;
                    }
                }
            } else if (bc_tags[face_bcs[f]] == NEUMANN){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    int ii=0;
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = normal[pos] * hx[pos]/2.0;
                        } else {
                            //move along face
                            x[pos] = hx[pos]/2.0 * quad_points[q][ii];
                            ii++;
                        }
                    }

                    if (e_minus > -1) {
                        //calculate A properties
                        rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        u_minus = (driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x)/rho_minus;
                        flux = rho_minus*area/num_quad*u_minus.dot(normal);
                        result(e_minus) -= flux;
                    }

                    //relative position to centroid of B
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = -normal[pos] * hx[pos]/2.0;
                        }
                    }

                    if (e_plus > -1) {
                        //calculate B properties
                        rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        u_plus = (driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x)/rho_plus;
                        flux = rho_plus*area/num_quad*u_plus.dot(normal);
                        result(e_plus) += flux;
                    }
                }
            }
        }
    }
    return result;
}

//functions to compute element momentum fluxes (including tractions)
KinematicVectorArray FVMCartesian::calculateElementMomentumFluxes(Job* job, FiniteVolumeDriver* driver){
    //flux rates
    KinematicVectorArray result = KinematicVectorArray(element_count, job->JOB_TYPE);
    result.setZero();

    //intermediate variables
    KinematicVector u_plus, u_minus, normal, u_bar, flux; //velocity
    KinematicVector p_plus, p_minus;
    int e_plus, e_minus;    //elements
    double rho_plus, rho_minus;
    double area;
    double rho_bar, c, a_1, a_2, lambda_1, lambda_2;
    KinematicVector x_0 = KinematicVector(job->JOB_TYPE);
    KinematicVector x_face = KinematicVector(job->JOB_TYPE);
    KinematicVector x_quad = KinematicVector(job->JOB_TYPE);
    KinematicVector a_3 = KinematicVector(job->JOB_TYPE);

    //traction calculatons
    KinematicTensorArray L = getVelocityGradients(job, driver);
    KinematicTensor L_tmp = KinematicTensor(job->JOB_TYPE);
    MaterialTensor tau_plus, tau_minus;
    double P_plus, P_minus, theta_bar;

    //number of quadrature points per face
    int num_quad = 1;
    if (GRID_DIM == 2){
        num_quad = 2;
    } else if (GRID_DIM == 3){
        num_quad = 4;
    }

    //vector of quad point positions relative to face center
    std::array<std::array<double,2>,4> quad_points;
    quad_points[0][0] = -1.0/std::sqrt(3);  quad_points[0][1] = -1.0/std::sqrt(3);
    quad_points[1][0] = 1.0/std::sqrt(3);   quad_points[1][1] = -1.0/std::sqrt(3);
    quad_points[2][0] = -1.0/std::sqrt(3);  quad_points[2][1] = 1.0/std::sqrt(3);
    quad_points[3][0] = 1.0/std::sqrt(3);   quad_points[3][1] = 1.0/std::sqrt(3);

    //different flux calculation for different order approximations
    if (driver->order == 1){
        //rho, u, value constant within element
        //loop over faces
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x_face = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (face_bcs[f] == -1 || (bc_tags[face_bcs[f]] == PERIODIC && face_elements[f][0] > -1)){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];
                rho_minus = driver->fluid_body->rho(e_minus);
                rho_plus = driver->fluid_body->rho(e_plus);
                p_minus = driver->fluid_body->p[e_minus];
                p_plus = driver->fluid_body->p[e_plus];
                u_minus = p_minus/rho_minus; //u = p/rho
                u_plus  = p_plus/rho_plus;

                //approximate Roe advective rate
                u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                rho_bar = std::sqrt(rho_plus*rho_minus);
                theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                c = driver->fluid_material->getSpeedOfSound(job, driver, x_face, rho_bar, theta_bar);

                //roe eigenvalues
                lambda_1 = std::abs(u_bar.dot(normal) - c);
                lambda_2 = std::abs(u_bar.dot(normal) + c);
                if (lambda_1 < delta*c){
                    lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                }
                if (lambda_2 < delta*c){
                    lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                }

                //calculate Roe eigenvector coefficients
                a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                a_2 = (rho_plus - rho_minus) - a_1;
                a_3 = p_plus - p_minus - (rho_plus-rho_minus)*u_bar;
                a_3 = a_3 - a_3.dot(normal)*normal; //remove normal component of a_3 vector

                //flux in n direction
                flux = area*0.5*(p_plus*u_plus.dot(normal) + p_minus*u_minus.dot(normal)
                                                           - a_1*lambda_1*(u_bar - c*normal)
                                                             - a_2*lambda_2*(u_bar + c*normal)
                                                               - a_3*std::abs(u_bar.dot(normal)));

                //add tractions to momentum flux
                //for now use simple reconstruction of theta
                tau_minus = driver->fluid_material->getShearStress(job, driver, x_face, L[e_minus], rho_bar, theta_bar);
                tau_plus  = driver->fluid_material->getShearStress(job, driver, x_face, L[e_plus], rho_bar, theta_bar);
                P_plus = c*c*(rho_plus - rho_bar) + driver->fluid_material->getPressure(job, driver, x_face, rho_bar, theta_bar);
                P_minus = c*c*(rho_minus - rho_plus) + P_plus;

                flux += area * 0.5 * ((P_plus + P_minus)*normal - KinematicVector((tau_plus + tau_minus)*normal, job->JOB_TYPE));

                //add flux to element integrals
                result(e_minus) -= flux;
                result(e_plus) += flux;
            } else if (bc_tags[face_bcs[f]] == DIRICHLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and assigned velocity
                    p_minus = driver->fluid_body->p[e_minus];

                    //tractions
                    rho_bar = driver->fluid_body->rho(e_minus);
                    theta_bar = driver->fluid_body->theta(e_minus);

                    //estimate L
                    x_0 = getElementCentroid(job, e_minus);
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_values[face_bcs[f]][ii] - p_minus[ii]/rho_bar)/(x_face - x_0).dot(normal)*normal[jj];
                        }
                    }

                    tau_minus = driver->fluid_material->getShearStress(job, driver, x_face, L_tmp, rho_bar, theta_bar);
                    P_minus = driver->fluid_material->getPressure(job, driver, x_face, rho_bar, theta_bar);

                    flux = p_minus*area*bc_values[face_bcs[f]].dot(normal)
                           + area*P_minus*normal
                           - area*KinematicVector(tau_minus*normal, job->JOB_TYPE);
                    result(e_minus) -= flux;
                }

                if (e_plus > -1){
                    //flux defined by cell value and assigned velocity
                    p_plus = driver->fluid_body->p[e_plus];

                    //tractions
                    rho_bar = driver->fluid_body->rho(e_plus);
                    theta_bar = driver->fluid_body->theta(e_plus);

                    //estimate L
                    x_0 = getElementCentroid(job, e_plus);
                    for (int ii=0; ii<GRID_DIM; ii++){
                        for (int jj=0; jj<GRID_DIM; jj++){
                            L_tmp(ii,jj) = (bc_values[face_bcs[f]][ii] - p_plus[ii]/rho_bar)/(x_face - x_0).dot(normal)*normal[jj];
                        }
                    }

                    tau_plus = driver->fluid_material->getShearStress(job, driver, x_face, L_tmp, rho_bar, theta_bar);
                    P_plus = driver->fluid_material->getPressure(job, driver, x_face, rho_bar, theta_bar);

                    flux = p_plus*area*bc_values[face_bcs[f]].dot(normal)
                           + area*P_plus*normal
                           - area*KinematicVector(tau_plus*normal, job->JOB_TYPE);

                    result(e_plus) += flux;
                }
            } else if (bc_tags[face_bcs[f]] == NEUMANN){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                if (e_minus > -1){
                    //flux defined by cell value and cell velocity
                    rho_minus = driver->fluid_body->rho(e_minus);
                    p_minus = driver->fluid_body->p[e_minus];
                    u_minus = p_minus/rho_minus;
                    flux = p_minus*area*u_minus.dot(normal);

                    //add traction directly to flux integral
                    result(e_minus) -= flux;
                    result(e_minus) += area*bc_values[face_bcs[f]];
                }

                if (e_plus > -1){
                    //flux defined by cell value and cell velocity
                    rho_plus = driver->fluid_body->rho(e_plus);
                    p_plus = driver->fluid_body->p[e_plus];
                    u_plus = p_plus/rho_plus;
                    flux = p_plus*area*u_plus.dot(normal);

                    //add traction directly to flux integral
                    result(e_plus) += flux;
                    result(e_minus) += area*bc_values[face_bcs[f]];
                }
            }
        }
    } else if (driver->order >= 2){
        if (driver->order > 2) {
            std::cout << "ERROR: FVMCartesian does not currently implement higher order flux reconstruction."
                      << std::endl;
        }

        KinematicVector x = KinematicVector(job->JOB_TYPE);

        //loop over faces and use quadrature to reconstruct flux integral
        u_plus = KinematicVector(job->JOB_TYPE);
        u_minus = KinematicVector(job->JOB_TYPE);

        for (int f=0; f<face_count; f++){
            //face dimensions
            area = getFaceArea(f);
            normal = getFaceNormal(job, f);
            x_face = getFaceCentroid(job, f);

            //flux calculation depends on whether face is on boundary
            if (face_bcs[f] == -1 || (bc_tags[face_bcs[f]] == PERIODIC && face_elements[f][0] > -1)){
                //face is interior to domain; no BCs (or periodic)
                e_minus = face_elements[f][0];
                e_plus  = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    int ii=0;
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = normal[pos] * hx[pos]/2.0;
                            x_quad[pos] = x_face[pos];
                        } else {
                            //move along face
                            x[pos] = hx[pos]/2.0 * quad_points[q][ii];
                            x_quad[pos] = x_face[pos] + x[pos];
                            ii++;
                        }
                    }

                    //calculate A properties
                    rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                    p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;
                    u_minus = p_minus/rho_minus;

                    //relative position to centroid of B
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = -normal[pos] * hx[pos]/2.0;
                        }
                    }

                    //calculate B properties
                    rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                    p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;
                    u_plus = p_plus/rho_plus;

                    //approximate Roe advective rate
                    u_bar = (std::sqrt(rho_plus)*u_plus + std::sqrt(rho_minus)*u_minus)/(std::sqrt(rho_plus) + std::sqrt(rho_minus));
                    rho_bar = std::sqrt(rho_plus*rho_minus);
                    theta_bar = std::sqrt(driver->fluid_body->theta(e_plus) * driver->fluid_body->theta(e_minus));
                    c = driver->fluid_material->getSpeedOfSound(job, driver, x_quad, rho_bar, theta_bar);

                    //roe eigenvalues
                    lambda_1 = std::abs(u_bar.dot(normal) - c);
                    lambda_2 = std::abs(u_bar.dot(normal) + c);
                    if (lambda_1 < delta*c){
                        lambda_1 = 0.5*((lambda_1*lambda_1)/(delta*c) + (delta*c));
                    }
                    if (lambda_2 < delta*c){
                        lambda_2 = 0.5*((lambda_2*lambda_2)/(delta*c) + (delta*c));
                    }

                    //calculate Roe eigenvector coefficients
                    a_1 = 0.5*(rho_plus - rho_minus) - 1.0/(2.0*c)*(rho_plus*u_plus - rho_minus*u_minus - (rho_plus-rho_minus)*u_bar).dot(normal);
                    a_2 = (rho_plus - rho_minus) - a_1;
                    a_3 = p_plus - p_minus - (rho_plus-rho_minus)*u_bar;
                    a_3 = a_3 - a_3.dot(normal)*normal; //remove normal component of a_3 vector

                    /*
                    std::cout << "alpha_int: " << a_1 << " " << a_2 << std::endl;
                    std::cout << "quadrature_location: " << x[0] << ", " << x[1] << std::endl;
                    std::cout << "rho: " << rho_plus << ", " << rho_minus << std::endl;
                    std::cout << ", u: " << u_plus[0] << ", " << u_minus[0] << std::endl;
                    std::cout << "u_bar: " << u_bar[0] << ", c_bar: " << c << std::endl;
                     */

                    //flux in n direction
                    flux = area/num_quad * 0.5 * (p_plus*u_plus.dot(normal) + p_minus*u_minus.dot(normal)
                                                    - a_1*lambda_1*(u_bar - c*normal)
                                                    - a_2*lambda_2*(u_bar + c*normal)
                                                    - a_3*std::abs(u_bar.dot(normal)));

                    //add tractions to momentum flux
                    //for now use simple reconstruction of theta
                    tau_minus = driver->fluid_material->getShearStress(job, driver, x_quad, L[e_minus], rho_bar, theta_bar);
                    tau_plus  = driver->fluid_material->getShearStress(job, driver, x_quad, L[e_plus], rho_bar, theta_bar);
                    P_plus = c*c*(rho_plus - rho_bar) + driver->fluid_material->getPressure(job, driver, x_quad, rho_bar, theta_bar);
                    P_minus = c*c*(rho_minus - rho_plus) + P_plus;

                    flux += area/num_quad * 0.5 * ((P_plus + P_minus)*normal - KinematicVector((tau_plus + tau_minus)*normal, job->JOB_TYPE));

                    //add flux to element integrals
                    result(e_minus) -= flux;
                    result(e_plus) += flux;
                }
            } else if (bc_tags[face_bcs[f]] == DIRICHLET){
                //face has prescribed velocity
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    int ii=0;
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = normal[pos] * hx[pos]/2.0;
                            x_quad[pos] = x_face[pos];
                        } else {
                            //move along face
                            x[pos] = hx[pos]/2.0 * quad_points[q][ii];
                            x_quad[pos] = x_face[pos] + x[pos];
                            ii++;
                        }
                    }

                    if (e_minus > -1) {
                        //calculate A properties
                        p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;

                        //tractions
                        rho_bar = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        theta_bar = driver->fluid_body->theta(e_minus);

                        //estimate L
                        x_0 = getElementCentroid(job, e_minus);
                        for (int ii=0; ii<GRID_DIM; ii++){
                            for (int jj=0; jj<GRID_DIM; jj++){
                                L_tmp(ii,jj) = (bc_values[face_bcs[f]][ii] - driver->fluid_body->p[e_minus][ii]/rho_bar)
                                               /(x_face - x_0).dot(normal)*normal[jj];
                            }
                        }

                        tau_minus = driver->fluid_material->getShearStress(job, driver, x_quad, L_tmp, rho_bar, theta_bar);
                        P_minus = driver->fluid_material->getPressure(job, driver, x_quad, rho_bar, theta_bar);

                        flux = p_minus*area/num_quad*bc_values[face_bcs[f]].dot(normal)
                               + area/num_quad*P_minus*normal
                               - area/num_quad*KinematicVector(tau_minus*normal, job->JOB_TYPE);

                        result(e_minus) -= flux;
                    }

                    //relative position to centroid of B
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = -normal[pos] * hx[pos]/2.0;
                        }
                    }

                    if (e_plus > -1) {
                        //calculate B properties
                        p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;

                        //tractions
                        rho_bar = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        theta_bar = driver->fluid_body->theta(e_plus);

                        //estimate L
                        x_0 = getElementCentroid(job, e_plus);
                        for (int ii=0; ii<GRID_DIM; ii++){
                            for (int jj=0; jj<GRID_DIM; jj++){
                                L_tmp(ii,jj) = (bc_values[face_bcs[f]][ii] - driver->fluid_body->p[e_plus][ii]/rho_bar)
                                               /(x_face - x_0).dot(normal)*normal[jj];
                            }
                        }

                        tau_plus = driver->fluid_material->getShearStress(job, driver, x_quad, L_tmp, rho_bar, theta_bar);
                        P_plus = driver->fluid_material->getPressure(job, driver, x_quad, rho_bar, theta_bar);

                        flux = p_plus*area/num_quad*bc_values[face_bcs[f]].dot(normal)
                               + area/num_quad*P_plus*normal
                               - area/num_quad*KinematicVector(tau_plus*normal, job->JOB_TYPE);

                        result(e_plus) += flux;
                    }
                }
            } else if (bc_tags[face_bcs[f]] == NEUMANN){
                //face has prescribed traction
                e_minus = face_elements[f][0];
                e_plus = face_elements[f][1];

                //loop over quadrature points
                for (int q=0; q<num_quad; q++) {
                    //relative position to centroid of A
                    int ii=0;
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = normal[pos] * hx[pos]/2.0;
                        } else {
                            //move along face
                            x[pos] = hx[pos]/2.0 * quad_points[q][ii];
                            ii++;
                        }
                    }

                    if (e_minus > -1) {
                        //calculate A properties
                        rho_minus = driver->fluid_body->rho(e_minus) + driver->fluid_body->rho_x[e_minus].dot(x);
                        p_minus = driver->fluid_body->p[e_minus] + driver->fluid_body->p_x[e_minus]*x;
                        u_minus = p_minus/rho_minus;
                        flux = p_minus*area/num_quad*u_minus.dot(normal);

                        //add traction directly
                        result(e_minus) -= flux;
                        result(e_minus) += area/num_quad*bc_values[face_bcs[f]];
                    }

                    //relative position to centroid of B
                    for (int pos=0; pos<GRID_DIM; pos++){
                        if (std::abs(normal[pos]) > 0.5) {
                            //move to face
                            x[pos] = -normal[pos] * hx[pos]/2.0;
                        }
                    }

                    if (e_plus > -1) {
                        //calculate B properties
                        rho_plus = driver->fluid_body->rho(e_plus) + driver->fluid_body->rho_x[e_plus].dot(x);
                        p_plus = driver->fluid_body->p[e_plus] + driver->fluid_body->p_x[e_plus]*x;
                        u_plus = p_plus/rho_plus;
                        flux = p_plus*area/num_quad*u_plus.dot(normal);

                        //add traction directly
                        result(e_plus) += flux;
                        result(e_minus) += area/num_quad*bc_values[face_bcs[f]];
                    }
                }
            }
        }
    }
    return result;
}
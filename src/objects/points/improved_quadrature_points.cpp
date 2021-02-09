//
// Created by aaron on 7/31/19.
// improved_quadrature_points.cpp
//


#include <iostream>
#include <iomanip>
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
#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "job.hpp"

#include "points.hpp"
#include "objects/bodies/bodies.hpp"
#include "objects/grids/grids.hpp"

static bool DEBUG_IMPROVED_QUAD = false;

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
/*
    * static const int STANDARD = 0;
    * static const int UGIMP = 1;
    * static const int CPGIMP = 2;
    * static const int CPDI = 3;
    * static const int CPDI2 = 4;
    *
    * static const int NO_CORRECTION = 0;
    * static const int AVAV = 1;
    * static const int SPH_LIKE = 2;
    * static const int DELTA_STRAIN = 3;
    * static const int DELTA_DISP = 4;
    * static const int DELTA_NEWTON = 5;
     */

/*--------*/
//search cell functions for avoid-a-void
int ImprovedQuadraturePoints::pos_to_cell(Job* job, KinematicVector x){
    KinematicVector tmpX = x - x_min;
    return CartesianLinear::cartesianWhichElement(job, tmpX, Lx, hx, Nx, job->grid->GRID_DIM);
}

std::vector<int> ImprovedQuadraturePoints::cell_to_ijk(Job* job, int i){
    std::vector<int> ijk = std::vector<int>(job->grid->GRID_DIM);
    int tmp = i;
    //find i,j,k count for element position
    for (int i=0;i<ijk.size();i++){
        ijk[i] = tmp % Nx(i);
        tmp = tmp / Nx(i);
    }

    return ijk;
}

int ImprovedQuadraturePoints::ijk_to_cell(Job* job, std::vector<int> ijk){
    int cell = 0;
    for (int i=0;i<job->grid->GRID_DIM;i++) {
        //n = i + j*imax + k*imax*jmax
        //hardcode
        if (i==0) {
            cell += ijk[i];
        } else if (i==1) {
            cell += ijk[i] * Nx(0);
        } else if (i==2) {
            cell += ijk[i] * Nx(0) * Nx(1);
        }
    }
    return cell;
}

/*----------------------------------------------------------------------------*/
//initialize point state (assumes that readFromFile has been called)
//no safety check on this, so be careful please
void ImprovedQuadraturePoints::init(Job* job, Body* body){
    int int_prop_counter = 2;
    int fp64_prop_counter = 0;
    if (int_props.size() < 2 || str_props.size() < 1){
        std::cout << int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 3 properties defined (QUADRULE, outputFolder).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //assign quadrature rule
        QUADRULE = int_props[0];
        //assing position rule
        POSITIONRULE = int_props[1];
        //assign output folder
        outputFolder = Parser::makeDirectory(str_props[0]);
        //counter for mixing properties from different flags
        int_prop_counter = 2;
        fp64_prop_counter = 0;
    }

    //all points require v0 initialization
    v0 = v;

    //initialize error measures
    V_i.resize(job->grid->node_count,1);
    v_i.resize(job->grid->node_count,1);
    e.resize(job->grid->node_count,1);
    grad_e = KinematicVectorArray(body->points->x.size(),job->JOB_TYPE);
    eonV.resize(job->grid->node_count,1);
    grad_eonV = KinematicVectorArray(body->points->x.size(),job->JOB_TYPE);
    H.resize(job->grid->node_count,1);
    grad_H = KinematicVectorArray(body->points->x.size(),job->JOB_TYPE);

    for (int i=0; i<job->grid->node_count;i++){
        V_i(i) = job->grid->nodeVolume(job,i);
    }

    //all points will need to know the potential number of corners
    int cp = 1;
    for (int i = 0; i < job->grid->GRID_DIM; i++) {
        cp *= 2; //square or cube
    }
    cpmp = cp;

    //initialize 'extent' of particle boxes
    if (QUADRULE == UGIMP || QUADRULE == CPGIMP || QUADRULE == CPDI || QUADRULE == CPDI2) {

        //adjust width of integration region (only for uGIMP)
        double scale_factor = 1.0;
        if (fp64_props.size() > (1+fp64_prop_counter) && QUADRULE == UGIMP) {
            scale_factor = fp64_props[0];
            fp64_prop_counter += 1;
        }

        //initialize extent of GIMP box
        if (job->grid->GRID_DIM == 1) {
            for (int i = 0; i < v.rows(); i++) {
                extent[i] = scale_factor * v[i];
            }
        } else if (job->JOB_TYPE == job->JOB_AXISYM) {
            for (int i = 0; i < v.rows(); i++) {
                extent[i] = scale_factor * std::sqrt(v[i] / x(i, 0));
            }
        } else if (job->grid->GRID_DIM == 2) {
            for (int i = 0; i < v.rows(); i++) {
                extent[i] = scale_factor * std::sqrt(v[i]);
            }
        } else if (job->grid->GRID_DIM == 3) {
            for (int i = 0; i < v.rows(); i++) {
                extent[i] = scale_factor * std::cbrt(v[i]);
            }
        }

        //initialize A matrix for corner mapping
        //for mapping position in cube to id
        //0 -> -1,-1,-1
        //1 -> +1,-1,-1
        //...
        //8 -> +1,+1,+1
        Eigen::VectorXi onoff = -1 * Eigen::VectorXi::Ones(job->grid->GRID_DIM);
        A = Eigen::MatrixXi(cp, job->grid->GRID_DIM);
        for (int c = 0; c < cp; c++) {
            for (int i = 0; i < onoff.rows(); i++) {
                A(c, i) = onoff(i);
            }
            for (int i = 0; i < onoff.rows(); i++) {
                if (onoff(i) == -1) {
                    onoff(i) = 1;
                    break;
                } else {
                    onoff(i) = -1;
                }
            }
        }
    }

    //CPGIMP and CPDI require F matrix for determining corner positions
    if (QUADRULE == CPGIMP || QUADRULE == CPDI){
        F = KinematicTensorArray(x.size(), x.VECTOR_TYPE);
        for (int i=0; i<F.size(); i++){
            //initialize F to identity
            F[i] = KinematicTensor::Identity(x.VECTOR_TYPE);
        }
    }

    //CPDI2 requires full list of corner positions and associated mapping matrix
    if (QUADRULE == CPDI2){
        //CPDI2 only defined on 2 dimensional grids
        if (job->grid->GRID_DIM != 2){
            std::cerr << "CPDI2 intergation scheme only implemented in 2D (GRID_DIM = " << job->grid->GRID_DIM << "). Exiting." << std::endl;
            exit(0);
        }

        //initialize list of corner position vectors
        corner_positions = KinematicVectorArray(x.size()*cpmp,x.VECTOR_TYPE);
        for (int i = 0; i < x.size(); i++) {
            for (int c = 0; c < A.rows(); c++) {
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++) {
                    corner_positions(i*cpmp + c,pos) = x(i,pos) + 0.5 * extent(i) * A(c, pos);
                }
                for (int pos = job->grid->GRID_DIM; pos < x.DIM; pos++) {
                    corner_positions(i*cpmp + c,pos) = 0;
                }
            }
        }
    }

    if (POSITIONRULE == DELTA_STRAIN
        || POSITIONRULE == DELTA_DISP
        || POSITIONRULE == SPH_LIKE
        || POSITIONRULE == DELTA_NEWTON) {
        //initialize position correction measure
        del_pos = KinematicVectorArray(body->points->x.size(), job->JOB_TYPE);
        del_pos.setZero();
    }

    //avoid-a-void algorithm requires radius, merge_dist
    if (POSITIONRULE == AVAV){
        if (fp64_props.size() < (fp64_prop_counter+3) || int_props.size() < (1+int_prop_counter)){
            std::cout << fp64_props.size() << ", " << int_props.size() << "\n";
            fprintf(stderr,
                    "%s:%s: Need at least 4 properties defined (r, merge_factor, buffer_scale) and skip value.\n",
                    __FILE__, __func__);
            exit(0);
        } else {
            //assign values
            r = fp64_props[fp64_prop_counter];
            merge_dist = r*fp64_props[fp64_prop_counter + 1];
            buffer_scale = fp64_props[fp64_prop_counter + 2];
            fp64_prop_counter += 3;

            skip_value = int_props[int_prop_counter];
            skip_counter = 0;
            int_prop_counter += 1;
        }

        //initialize grid variables
        d = Eigen::VectorXd(job->grid->node_count);
        node_to_cell_map = std::vector<int>(job->grid->node_count);

        //edge_list
        edge_list = job->grid->getEdgeList(job);

        //increase point size by buffer
        int old_len = x.size();
        int len = (int)(x.size()*(1.0 + buffer_scale));

        //allocate point to cell map
        point_to_cell_map = std::vector<int>(len);

        //allocate space for extra points
        x.resize(len);
        u.resize(len);
        x_t.resize(len);
        mx_t.resize(len);
        b.resize(len);

        //size scalar vectors (need to be careful with these)
        m.conservativeResize(len);
        v.conservativeResize(len);
        v0.conservativeResize(len);
        active.conservativeResize(len);
        extent.conservativeResize(len);

        //size tensor arrays
        T.resize(len);
        L.resize(len);

        //initialize point variables
        point_dist = Eigen::VectorXd(len);

        //if using CPDI2, do this:
        if (QUADRULE == CPDI2) {
            corner_positions.resize(x.size() * cpmp);
        }

        //fill active with zeros and update buffer indices
        buffer_list.resize(len - old_len);
        for (int i=old_len; i<len; i++){
            buffer_list[i-old_len] = i;
            active(i) = 0;
        }

        //generate search cell definitions
        x_min = job->grid->nodeIDToPosition(job,0);
        x_max = job->grid->nodeIDToPosition(job,0);
        for (int i=1; i<job->grid->node_count; i++){
            for (int pos=0; pos<job->grid->GRID_DIM; pos++){
                if (job->grid->nodeIDToPosition(job,i)[pos] < x_min[pos]){
                    x_min[pos] = job->grid->nodeIDToPosition(job,i)[pos];
                } else if (job->grid->nodeIDToPosition(job,i)[pos] > x_max[pos]){
                    x_max[pos] = job->grid->nodeIDToPosition(job,i)[pos];
                }
            }
        }
        Lx = KinematicVector(job->JOB_TYPE);
        hx = KinematicVector(job->JOB_TYPE);
        Nx = Eigen::VectorXi(job->grid->GRID_DIM);
        double cell_width = 2.0*r;
        for (int pos=0; pos<job->grid->GRID_DIM; pos++){
            Lx(pos) = x_max[pos] - x_min[pos];
            Nx(pos) = (int)(Lx(pos) / cell_width);
            hx(pos) = Lx(pos)/Nx(pos);
        }

        //allocate space for cell to point map
        int cell_count = 1;
        for (int pos=0; pos<job->grid->GRID_DIM; pos++){
            cell_count *= Nx(pos);
        }
        cell_to_point_map = std::vector<std::vector<int>>(cell_count);
        cell_to_node_map = std::vector<std::vector<int>>(cell_count);
        for (int cell=0; cell<cell_count; cell++){
            cell_to_point_map[cell] = std::vector<int>(0);
            cell_to_node_map[cell] = std::vector<int>(0);
        }

        //get list of search cell ids
        //for (int i=0; i<job->grid->node_count; i++){
        //    node_to_cell_map[i] = pos_to_cell(job, body->nodes->x[i]);
        //}
        int tmp_id = -1;
        KinematicVector tmpPos;
        for (int i=0; i<job->grid->node_count; i++){
            tmp_id = pos_to_cell(job, job->grid->nodeIDToPosition(job,i));
            if (tmp_id < 0){
                //check if point is merely too close to boundary
                tmpPos = job->grid->nodeIDToPosition(job,i);
                for (int pos = 0; pos<job->grid->GRID_DIM; pos++){
                    tmpPos[pos] *= (1.0-1e-10);
                }
                tmp_id = pos_to_cell(job, tmpPos);
            }
            if (tmp_id < 0){
                node_to_cell_map[i] = tmp_id;
            } else {
                node_to_cell_map[i] = tmp_id;
                cell_to_node_map[tmp_id].push_back(i);
            }
        }

        //set number of neighbors for search
        number_of_neighbors = 1;
        number_of_second_neighbors = 1;
        for (int dir=0; dir<job->grid->GRID_DIM; dir++){
            number_of_neighbors *= 3;
            number_of_second_neighbors *= 5;
        }

    }


    //sph-like algorithm requires radius, constant, skip counter
    if (POSITIONRULE == SPH_LIKE){
        if (fp64_props.size() < (fp64_prop_counter+2) || int_props.size() < (1+int_prop_counter)){
            std::cout << fp64_props.size() << ", " << int_props.size() << "\n";
            fprintf(stderr,
                    "%s:%s: Need at least 3 properties defined (r, C) and skip value.\n",
                    __FILE__, __func__);
            exit(0);
        } else {
            //assign values
            r = fp64_props[fp64_prop_counter];
            sph_const = r*fp64_props[fp64_prop_counter + 1];
            fp64_prop_counter += 2;

            skip_value = int_props[int_prop_counter];
            skip_counter = 0;
            int_prop_counter += 1;
        }

        //allocate point to cell map
        point_to_cell_map = std::vector<int>(x.size());

        //generate search cell definitions
        x_min = job->grid->nodeIDToPosition(job,0);
        x_max = job->grid->nodeIDToPosition(job,0);
        for (int i=1; i<job->grid->node_count; i++){
            for (int pos=0; pos<job->grid->GRID_DIM; pos++){
                if (job->grid->nodeIDToPosition(job,i)[pos] < x_min[pos]){
                    x_min[pos] = job->grid->nodeIDToPosition(job,i)[pos];
                } else if (job->grid->nodeIDToPosition(job,i)[pos] > x_max[pos]){
                    x_max[pos] = job->grid->nodeIDToPosition(job,i)[pos];
                }
            }
        }
        Lx = KinematicVector(job->JOB_TYPE);
        hx = KinematicVector(job->JOB_TYPE);
        Nx = Eigen::VectorXi(job->grid->GRID_DIM);
        double cell_width = 2.0*r;
        for (int pos=0; pos<job->grid->GRID_DIM; pos++){
            Lx(pos) = x_max[pos] - x_min[pos];
            Nx(pos) = (int)(Lx(pos) / cell_width);
            hx(pos) = Lx(pos)/Nx(pos);
        }

        //allocate space for cell to point map
        int cell_count = 1;
        for (int pos=0; pos<job->grid->GRID_DIM; pos++){
            cell_count *= Nx(pos);
        }
        cell_to_point_map = std::vector<std::vector<int>>(cell_count);
        for (int cell=0; cell<cell_count; cell++){
            cell_to_point_map[cell] = std::vector<int>(0);
        }

        int tmp_id = -1;
        KinematicVector tmpPos;

        //set number of neighbors for search
        number_of_neighbors = 1;
        number_of_second_neighbors = 1;
        for (int dir=0; dir<job->grid->GRID_DIM; dir++){
            number_of_neighbors *= 3;
            number_of_second_neighbors *= 5;
        }

    }

    //delta correction for constant strain and constant displacement requires h and alpha
    if (POSITIONRULE == DELTA_STRAIN || POSITIONRULE == DELTA_DISP){
        if (fp64_props.size() < 2+fp64_prop_counter){
            std::cout << fp64_props.size() << "\n";
            fprintf(stderr,
                    "%s:%s: Need at least 2 properties defined (h, alpha).\n",
                    __FILE__, __func__);
            exit(0);
        } else {
            //assign grid scale and rate scale
            h = fp64_props[fp64_prop_counter];
            alpha = fp64_props[fp64_prop_counter+1];
            fp64_prop_counter += 2;
        }
    }

    if (POSITIONRULE == DELTA_NEWTON){
        if (fp64_props.size() < 3+fp64_prop_counter || int_props.size() < 1+int_prop_counter){
            std::cout << fp64_props.size() << ", " << int_props.size() << "\n";
            fprintf(stderr,
                    "%s:%s: Need at least 3 properties defined (h, alpha, max_iter).\n",
                    __FILE__, __func__);
            exit(0);
        } else {
            //assign grid scale and rate scale
            h = fp64_props[fp64_prop_counter];
            alpha = fp64_props[fp64_prop_counter+1];
            REL_TOL = fp64_props[fp64_prop_counter+2];
            fp64_prop_counter += 3;

            max_iter = int_props[int_prop_counter];
            int_prop_counter += 1;
        }
    }

    //Avoid-a-void not implemented yet
    if ((QUADRULE != STANDARD) &&
            (QUADRULE != UGIMP) &&
            (QUADRULE != CPGIMP) &&
            (QUADRULE != CPDI) &&
            (QUADRULE != CPDI2)){
        std::cerr << "ERROR: Quadrature Rule " << QUADRULE << " not implemented in ImprovedQuadraturePoints yet. Exiting." << std::endl;
        exit(0);
    }

    if ((POSITIONRULE != NO_CORRECTION)
        && (POSITIONRULE != AVAV)
        && (POSITIONRULE != SPH_LIKE)
        && (POSITIONRULE != DELTA_STRAIN)
        && (POSITIONRULE != DELTA_DISP)
        && (POSITIONRULE != DELTA_NEWTON)
        ){
        std::cerr << "ERROR: Position Rule " << POSITIONRULE << " not implemented in ImprovedQuadraturePoints yet. Exiting." << std::endl;
        exit(0);
    }

    if (int_props.size() < int_prop_counter+1){
        //standard behavior
        use_elem = false;
    } else if (int_props[int_prop_counter] == 1){
        //create point-wise element list
        use_elem = true;
        elem = Eigen::MatrixXi(x.size(),cp);
        elem.setConstant(-1);
        std::cout << object_name << " using point-element history." << std::endl;
    }

    std::cout << "Points Initialized: [" << file << "]. Using Quadrature Rule " << QUADRULE << ".";
    std::cout << " Using Position Rule " << POSITIONRULE << "." << std::endl;
    std::cout << " * 0 -- standard MPM\n";
    std::cout << " * 1 -- uGIMP\n";
    std::cout << " * 2 -- cpGIMP\n";
    std::cout << " * 3 -- CPDI\n";
    std::cout << " * 4 -- CPDI2\n";
    std::cout << " * \n";
    std::cout << " * 0 -- standard MPM\n";
    std::cout << " * 1 -- Avoid-a-void\n";
    std::cout << " * 2 -- SPH-like shifting\n";
    std::cout << " * 3 -- \\delta-correction, strain" << std::endl;
    std::cout << " * 4 -- \\delta-correction, displacement" << std::endl;

    return;
}

void ImprovedQuadraturePoints::generateMap(Job* job, Body* body, int SPEC) {
    //the default body will defer to points to generate the mapping

    body->S.clear();
    body->gradS.clear();

    //calculate phi and grad phi
    std::vector<int> nvec(0);
    std::vector<double> valvec(0);
    KinematicVectorArray gradvec(0,job->JOB_TYPE);
    KinematicVector tmpGrad(job->JOB_TYPE);
    KinematicVector tmpGrad2(job->JOB_TYPE);
    KinematicVector tmpVec(job->JOB_TYPE);

    //cpdi2 variables
    double x1, x2, x3, x4, y1, y2, y3, y4, vp;
    double a1, a2, a3, b1, b2, b3;

    int ith_cpdi;

    if (QUADRULE == STANDARD) {
        /*---------------------------------------------------------------------------------------*/
        //standard mpm
        for (int i = 0; i < x.size(); i++) {
            //store point position
            KinematicVector::Map x_i = x[i];

            //check whether point is active:
            if (active(i) == 0) {
                continue;
            } else if (!job->grid->inDomain(job, x_i)) {
                active(i) = 0;
                continue;
            } else {
                //initialize std::vectors
                nvec.resize(0);
                valvec.resize(0);

                //if using element histroy,
                if (use_elem) {
                    job->grid->evaluateBasisFnValue(job, x_i, nvec, valvec, elem(i, 0));
                } else {
                    job->grid->evaluateBasisFnValue(job, x_i, nvec, valvec);
                }
                for (int j = 0; j < nvec.size(); j++) {
                    body->S.push_back(nvec[j], i, valvec[j]); //node, point, value
                }

                nvec.resize(0);
                gradvec.resize(0);
                if (use_elem) {
                    job->grid->evaluateBasisFnGradient(job, x_i, nvec, gradvec, elem(i, 0));
                } else {
                    job->grid->evaluateBasisFnGradient(job, x_i, nvec, gradvec);
                }
                for (int j = 0; j < nvec.size(); j++) {
                    body->gradS.push_back(nvec[j], i, gradvec[j]); //node, point, value
                }
            }
        }
    } else if (QUADRULE == UGIMP || QUADRULE == CPGIMP || QUADRULE == CPDI){
        /*---------------------------------------------------------------------------------------*/
        //ugimp, cpgimp, cpdi
        for (int i = 0; i < x.size(); i++) {
            //store point position
            KinematicVector::Map x_i = x[i];

            //initialize cpdi flag
            ith_cpdi = DefaultBody::CPDI_ON;

            //check whether point is active:
            if (active(i) == 0) {
                continue;
            } else if (!job->grid->inDomain(job, x_i)) {
                active(i) = 0;
                continue;
            } else {
                //check that corners are in domain
                for (int c = 0; c < A.rows(); c++) {
                    //determine corner positions
                    if (QUADRULE == UGIMP){
                        //corners are fixed distance from point center
                        for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                            tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos);
                        }
                    } else if (QUADRULE == CPGIMP){
                        //corners are stretched along cartesian axes
                        for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                            tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos)*F(i,pos,pos);
                        }
                    } else if (QUADRULE == CPDI){
                        //corners are defined by x_c = F \xi_c
                        for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                            tmpVec[pos] = x_i[pos];
                            for (int ii=0; ii<job->grid->GRID_DIM; ii++){
                                tmpVec[pos] += 0.5*extent(i)*A(c,ii)*F(i,pos,ii);
                            }
                        }
                    }
                    for (int pos = job->grid->GRID_DIM; pos < tmpVec.DIM; pos++){
                        tmpVec(pos) = 0;
                    }
                    if ((use_elem && !job->grid->inDomain(job, tmpVec, elem(i,c))) || !job->grid->inDomain(job, tmpVec)) {
                        //corner out of domain
                        ith_cpdi = DefaultBody::CPDI_OFF;
                        break;
                    }
                }

                //initialize std::vectors
                nvec.resize(0);
                valvec.resize(0);

                //calculate map for ith point
                if (ith_cpdi == DefaultBody::CPDI_OFF) {
                    //use standard mpm scheme for this point
                    nvec.resize(0);
                    valvec.resize(0);
                    if (use_elem) {
                        job->grid->evaluateBasisFnValue(job, x_i, nvec, valvec, elem(i,0));
                        for (int c=1; c<A.rows(); c++){
                            //assign all corners to centroid value
                            elem(i,c) = elem(i,0);
                        }
                    } else {
                        job->grid->evaluateBasisFnValue(job, x_i, nvec, valvec);
                    }
                    for (int j = 0; j < nvec.size(); j++) {
                        body->S.push_back(nvec[j], i, valvec[j]); //node, point, value
                    }

                    nvec.resize(0);
                    gradvec.resize(0);
                    if (use_elem) {
                        job->grid->evaluateBasisFnGradient(job, x_i, nvec, gradvec, elem(i,0));
                    } else {
                        job->grid->evaluateBasisFnGradient(job, x_i, nvec, gradvec);
                    }
                    for (int j = 0; j < nvec.size(); j++) {
                        body->gradS.push_back(nvec[j], i, gradvec[j]); //node, point, value
                    }
                } else if (ith_cpdi == DefaultBody::CPDI_ON) {
                    //initialize storage vectors
                    nvec.resize(0);
                    valvec.resize(0);

                    //add contribution of each corner to S
                    for (int c = 0; c < A.rows(); c++) {
                        //determine corner positions
                        if (QUADRULE == UGIMP){
                            //corners are fixed distance from point center
                            for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                                tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos);
                            }
                        } else if (QUADRULE == CPGIMP){
                            //corners are stretched along cartesian axes
                            for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                                tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos)*F(i,pos,pos);
                            }
                        } else if (QUADRULE == CPDI){
                            //corners are defined by x_c = F \xi_c
                            for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                                tmpVec[pos] = x_i[pos];
                                for (int ii=0; ii<job->grid->GRID_DIM; ii++){
                                    tmpVec[pos] += 0.5*extent(i)*A(c,ii)*F(i,pos,ii);
                                }
                            }
                        }

                        //zero out of plane dimension
                        for (int pos = job->grid->GRID_DIM; pos < tmpVec.DIM; pos++){
                            tmpVec(pos) = 0;
                        }

                        if (use_elem) {
                            job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec, elem(i,c));
                        } else {
                            job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
                        }
                    }

                    //insert mappings into mapping matrix with appropriate weights
                    /*
                    double detF = 1;
                    if (QUADRULE == UGIMP) {
                        detF = 1;
                    } else if (QUADRULE == CPGIMP){
                        for (int ii=0; ii<job->grid->GRID_DIM; ii++){
                            detF *= F(i,ii,ii);
                        }
                    } else if (QUADRULE == CPDI){
                        detF = F(i).det();
                    }*/

                    for (int j = 0; j < nvec.size(); j++) {
                        //body->S.push_back(nvec[j], i, detF * valvec[j] / A.rows()); //node, point, value
                        body->S.push_back(nvec[j], i, valvec[j] / A.rows()); //node, point, value
                    }

                    //determine gradient using quadrature rules
                    nvec.resize(0);
                    valvec.resize(0);
                    gradvec.resize(0);
                    for (int c = 0; c < A.rows(); c++) {
                        //find shape function value at offset point location
                        //add node ids to nodevec
                        //add values to valvec
                        valvec.resize(0);

                        //determine corner positions
                        if (QUADRULE == UGIMP){
                            //corners are fixed distance from point center
                            for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                                tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos);
                            }
                        } else if (QUADRULE == CPGIMP){
                            //corners are stretched along cartesian axes
                            for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                                tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos)*F(i,pos,pos);
                            }
                        } else if (QUADRULE == CPDI){
                            //corners are defined by x_c = F \xi_c
                            for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                                tmpVec[pos] = x_i[pos];
                                for (int ii=0; ii<job->grid->GRID_DIM; ii++){
                                    tmpVec[pos] += 0.5*extent(i)*A(c,ii)*F(i,pos,ii);
                                }
                            }
                        }

                        if (use_elem) {
                            job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec, elem(i,c));
                        } else {
                            job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
                        }

                        for (int v = 0; v < valvec.size(); v++) {
                            //gradient contribution from corner
                            //G(x) = (S(x+a) - S(x-a))/(2a)
                            for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                                tmpGrad[pos] = valvec[v] / (A.rows() * 0.5 * extent(i)) * A(c,pos);
                            }
                            for (int pos = job->grid->GRID_DIM; pos < tmpGrad.DIM; pos++){
                                tmpGrad[pos] = 0;
                            }

                            //adjust gradient
                            if (QUADRULE == UGIMP){
                                //do nothing
                            } else if (QUADRULE == CPGIMP){
                                for (int pos=0; pos<job->grid->GRID_DIM; pos++){
                                    tmpGrad[pos] /= F(i,pos,pos);
                                }
                            } else if (QUADRULE == CPDI){
                                tmpGrad2 = tmpGrad;
                                tmpGrad = (F(i).inverse())*tmpGrad2;
                            }

                            //add to point gradient
                            gradvec.push_back(tmpGrad);
                        }
                    }

                    for (int j = 0; j < nvec.size(); j++) {
                        body->gradS.push_back(nvec[j], i, gradvec[j]); //node, point, value
                    }
                }
            }
        }

    } else if (QUADRULE == CPDI2) {
        /*---------------------------------------------------------------------------------------*/
        //cpdi2
        for (int i=0; i<x.size(); i++){
            //store point centroid position
            KinematicVector::Map x_i = x[i];

            //initialize cpdi flag
            ith_cpdi = DefaultBody::CPDI_ON;

            //check whether point is active:
            if (active(i) == 0) {
                continue;
            } else if (!job->grid->inDomain(job, x_i)) {
                active(i) = 0;
                continue;
            } else {
                //check that corners are in domain
                for (int c = 0; c < cpmp; c++) {
                    //corner positions are stored over time
                    tmpVec = corner_positions(i * cpmp + c);

                    if ((use_elem && !job->grid->inDomain(job, tmpVec, elem(i, c))) ||
                        !job->grid->inDomain(job, tmpVec)) {
                        //corner out of domain
                        ith_cpdi = DefaultBody::CPDI_OFF;
                        break;
                    }
                }

                //initialize std::vectors
                nvec.resize(0);
                valvec.resize(0);

                //calculate map for ith point
                if (ith_cpdi == DefaultBody::CPDI_OFF) {
                    //use standard mpm scheme for this point
                    nvec.resize(0);
                    valvec.resize(0);
                    if (use_elem) {
                        job->grid->evaluateBasisFnValue(job, x_i, nvec, valvec, elem(i,0));
                        for (int c=1; c<A.rows(); c++){
                            //assign all corners to centroid value
                            elem(i,c) = elem(i,0);
                        }
                    } else {
                        job->grid->evaluateBasisFnValue(job, x_i, nvec, valvec);
                    }
                    for (int j = 0; j < nvec.size(); j++) {
                        body->S.push_back(nvec[j], i, valvec[j]); //node, point, value
                    }

                    nvec.resize(0);
                    gradvec.resize(0);
                    if (use_elem) {
                        job->grid->evaluateBasisFnGradient(job, x_i, nvec, gradvec, elem(i,0));
                    } else {
                        job->grid->evaluateBasisFnGradient(job, x_i, nvec, gradvec);
                    }
                    for (int j = 0; j < nvec.size(); j++) {
                        body->gradS.push_back(nvec[j], i, gradvec[j]); //node, point, value
                    }
                } else if (ith_cpdi == DefaultBody::CPDI_ON) {

                    //store corner positions
                    x1 = corner_positions(i*cpmp + 0,0);
                    x2 = corner_positions(i*cpmp + 1,0);
                    x3 = corner_positions(i*cpmp + 3,0);
                    x4 = corner_positions(i*cpmp + 2,0);

                    y1 = corner_positions(i*cpmp + 0,1);
                    y2 = corner_positions(i*cpmp + 1,1);
                    y3 = corner_positions(i*cpmp + 3,1);
                    y4 = corner_positions(i*cpmp + 2,1);

                    //calculate a and b constants
                    a1 = 0.25*(-x1 + x2 + x3 - x4);
                    a2 = 0.25*(-x1 - x2 + x3 + x4);
                    a3 = 0.25*( x1 - x2 + x3 - x4);

                    b1 = 0.25*(-y1 + y2 + y3 - y4);
                    b2 = 0.25*(-y1 - y2 + y3 + y4);
                    b3 = 0.25*( y1 - y2 + y3 - y4);

                    //calculate volume
                    vp = 4.0*(a1*b2 - a2*b1);

                    if (vp <= 0){
                        std::cout << "WARNING: Point has zero volume!." << std::endl;
                    }

                    //determine corner weights
                    for (int c = 0; c < A.rows(); c++) {
                        /* CALCULATE S_ip */
                        //initialize storage vectors
                        nvec.resize(0);
                        valvec.resize(0);
                        gradvec.resize(0);

                        //determine corner positions
                        tmpVec = corner_positions(i * cpmp + c);

                        //zero out of plane dimension
                        for (int pos = job->grid->GRID_DIM; pos < tmpVec.DIM; pos++) {
                            tmpVec(pos) = 0;
                        }

                        if (use_elem) {
                            job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec, elem(i, c));
                        } else {
                            job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
                        }

                        //corner weighting function
                        double w = a1*b2 - a2*b1 + ((a1*b3 - a3*b1)*A(c,0) + (a3*b2 - a2*b3)*A(c,1))/3.0;
                        //adjust weight by volume of point
                        w /= vp;

                        //add corner contribution to S and gradS
                        for (int j = 0; j < nvec.size(); j++) {

                            //mapping contribution from corner
                            body->S.push_back(nvec[j], i, w * valvec[j]);// / A.rows()); //node, point, value

                            //gradient contribution from corner
                            tmpGrad[0] = valvec[j]*(A(c,0)*b2 - A(c,1)*(a2 + a3*A(c,0)/3.0))/vp;
                            tmpGrad[1] = valvec[j]*(-A(c,0)*(b1 + b3*A(c,1)/3.0) + A(c,1)*a1)/vp;
                            for (int pos = job->grid->GRID_DIM; pos < tmpGrad.DIM; pos++) {
                                tmpGrad[pos] = 0;
                            }

                            //add to point gradient
                            //gradvec.push_back(tmpGrad);
                            body->gradS.push_back(nvec[j], i, tmpGrad);
                        }
                    }
                }
            }
        }
    }

    if (DEBUG_IMPROVED_QUAD) {
        //check for partition of unity in mapping
        Eigen::VectorXd unity_nodes = Eigen::VectorXd::Ones(body->nodes->x.size());
        Eigen::VectorXd unity_points = body->S.operate(unity_nodes, MPMSparseMatrixBase::TRANSPOSED);
        std::cout << "min: " << unity_points.minCoeff() << ", max: " << unity_points.maxCoeff() << std::endl;
    }

    //for all methods, calculate error measures
    //calculate nodal volume integral
    if (job->JOB_TYPE == job->JOB_AXISYM) {
        //need to adjust point integration area
        Eigen::VectorXd A_tmp = body->points->v;
        for (int i = 0; i < A_tmp.rows(); i++) {
            A_tmp(i) /= body->points->x(i, 0);
        }
        v_i = body->S * A_tmp;
    } else {
        //otherwise integration area and volume are the same
        v_i = body->S * body->points->v;
    }

    double tmpNum;
    for (int i = 0; i < V_i.rows(); i++) {
        tmpNum = (v_i(i) - V_i(i)) / (V_i(i));
        e(i) = std::max(0.0, tmpNum);
        eonV(i) = e(i)*std::sqrt(V_i(i));
        H(i) = e(i)*V_i(i);
    }

    //calculate gradient of nodal overshoot
    grad_e = body->gradS.operate(e, MPMSparseMatrixBase::TRANSPOSED);
    grad_eonV = body->gradS.operate(eonV, MPMSparseMatrixBase::TRANSPOSED);
    grad_H = body->gradS.operate(H, MPMSparseMatrixBase::TRANSPOSED);

    return;
}


/*----------------------------------------------------------------------------*/
//write relavent point data to file
void ImprovedQuadraturePoints::writeFrame(Job* job, Body* body, Serializer* serializer){
    //call default function
    DefaultPoints::writeFrame(job, body, serializer);

    //custom file output
    serializer->writeVectorArray(grad_e, "grad_err");
    serializer->writeScalarArray(V_i, "grid_volume");
    serializer->writeScalarArray(e, "err");
    serializer->writeVectorArray(grad_H, "grad_H");
    serializer->writeScalarArray(H, "H");

    if (POSITIONRULE == DELTA_STRAIN
        || POSITIONRULE == DELTA_DISP
        || POSITIONRULE == SPH_LIKE
        || POSITIONRULE == DELTA_NEWTON) {
        serializer->writeVectorArray(del_pos, "del_pos");
    } else if (POSITIONRULE == AVAV){
        point_dist = body->S.operate(d,MPMSparseMatrixBase::TRANSPOSED); //map distance field to points
        serializer->writeScalarArray(d, "grid_dist");
        serializer->writeScalarArray(point_dist, "point_dist");
    }

    //create temporary kinematicvector
    KinematicVector tmpVec(x.VECTOR_TYPE);
    //for uGIMP, cpGIMP, CPDI, and CPDI2, write output file with boxes for point domains
    if (QUADRULE == UGIMP || QUADRULE == CPGIMP || QUADRULE == CPDI || QUADRULE == CPDI2) {
        //open point file
        std::stringstream ss;
        ss << "fpd." << body->id << "." << body->name << ".corners." << std::setw(10) << std::setfill('0') << (sampledFrames)
           << ".vtk";
        //increment sampled frames
        sampledFrames++;

        std::string pfilename = ss.str();
        //pfile = std::ofstream(frameDirectory+pfilename,std::ios::trunc);
        std::ofstream pfile = std::ofstream(outputFolder + pfilename, std::ios::trunc);

        int plen = body->points->x.size();

        pfile << "# vtk DataFile Version 3.0\n" << "Frame: 0, Time: 0\n";
        pfile << "ASCII\n";
        pfile << "DATASET UNSTRUCTURED_GRID\n";

        pfile << "POINTS " << (plen * cpmp) << " double\n";
        for (int i = 0; i < plen; i++) {
            //store point position
            KinematicVector::Map x_i = x[i];

            for (int c = 0; c < cpmp; c++) {
                //determine corner positions
                if (QUADRULE == UGIMP){
                    //corners are fixed distance from point center
                    for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                        tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos);
                    }
                } else if (QUADRULE == CPGIMP){
                    //corners are stretched along cartesian axes
                    for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                        tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos)*F(i,pos,pos);
                    }
                } else if (QUADRULE == CPDI){
                    //corners are defined by x_c = F \xi_c
                    for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                        tmpVec[pos] = x_i[pos];
                        for (int ii=0; ii<job->grid->GRID_DIM; ii++){
                            tmpVec[pos] += 0.5*extent(i)*A(c,ii)*F(i,pos,ii);
                        }
                    }
                } else if (QUADRULE == CPDI2){
                    //corner positions defined by stored value
                    tmpVec = corner_positions[i*cpmp + c];
                }

                //vtk files require x,y,z
                for (int pos = 0; pos < 3; pos++) {
                    if (pos < body->points->x.DIM && (body->points->active(i) != 0) &&
                        std::isfinite(tmpVec[pos])) {
                        pfile << tmpVec[pos] << " ";
                    } else {
                        pfile << "0 ";
                    }
                }
                pfile << "\n";
            }
        }

        pfile << "CELLS " << plen << " " << (1 + cpmp) * plen << "\n";
        if (job->grid->GRID_DIM == 1){
            for (int i=0; i<plen; i++) {
                pfile << "2 " << i*cpmp << " "
                              << i*cpmp+1 << "\n";
            }
        } else if(job->grid->GRID_DIM == 2){
            for (int i=0; i<plen; i++) {
                pfile << "4 " << i*cpmp << " "
                              << i*cpmp+1 << " "
                              << i*cpmp+3 << " "
                              << i*cpmp+2 << "\n";
            }
        } else {
            for (int i=0; i<plen; i++) {
                pfile << "8 " << i*cpmp << " "
                               << i*cpmp+1 << " "
                               << i*cpmp+3 << " "
                               << i*cpmp+2 << " "
                               << i*cpmp+4 << " "
                               << i*cpmp+5 << " "
                               << i*cpmp+7 << " "
                               << i*cpmp+6 << "\n";
            }
        }

        pfile << "CELL_TYPES " << plen << "\n";
        if (job->grid->GRID_DIM == 1) {
            for (int i = 0; i < plen; i++) {
                pfile << "3\n";
            }
        } else if (job->grid->GRID_DIM == 2) {
            for (int i = 0; i < plen; i++) {
                pfile << "9\n";
            }
        } else {
            for (int i = 0; i < plen; i++) {
                pfile << "12\n";
            }
        }

        pfile << "POINT_DATA " << plen*cpmp << "\n";

        //write to point file
        pfile << "SCALARS volume double 1\n";
        pfile << "LOOKUP_TABLE default\n";
        for (int i = 0; i < plen; i++){
            for (int c=0; c<cpmp; c++){
                if (active(i) == 1 && std::isfinite(v(i))) {
                    pfile << v(i) << "\n";
                } else {
                    pfile << "0" << "\n";
                }
            }
        }

        pfile.close();
    }

    return;
}

void ImprovedQuadraturePoints::updateIntegrators(Job* job, Body* body){

    if (QUADRULE == STANDARD || QUADRULE == UGIMP){
        //do nothing
    } else if (QUADRULE == CPGIMP) {
        //\dot{F_ii} = L_ii * F_ii
        for (int i=0; i<F.size(); i++){
            for (int pos=0; pos<job->grid->GRID_DIM; pos++){
                //F^{n+1} = F^{n}*exp(L*dt)
                F(i,pos,pos) *= exp(job->dt*L(i,pos,pos));
            }

            //update volume
            v[i] = F[i].det() * v0[i];
        }
    } else if (QUADRULE == CPDI){
        KinematicTensor tmpL = KinematicTensor(F.TENSOR_TYPE);
        //\dot{F} = L*F
        for (int i=0; i<F.size(); i++){
            //F^{n+1} = exp(L*dt)*F^{n}
            tmpL = job->dt*L[i];
            F[i] = tmpL.exp(5)*F[i];

            //update volume
            v[i] = F[i].det() * v0[i];
        }
    } else if (QUADRULE == CPDI2){

        //mapping containers
        int ith_cpdi;
        KinematicVector tmpVec(x.VECTOR_TYPE);
        std::vector<int> nvec(0);
        std::vector<double> valvec(0);

        //cpdi2 variables
        double x1, x2, x3, x4, y1, y2, y3, y4, vp;
        double a1, a2, a3, b1, b2, b3;

        //update corner position for each point
        for (int i=0; i<x.size(); i++) {
            //store point centroid position
            KinematicVector::Map x_i = x[i];

            //initialize cpdi flag
            ith_cpdi = DefaultBody::CPDI_ON;

            //check whether point is active:
            if (active(i) == 0) {
                continue;
            } else if (!job->grid->inDomain(job, x_i)) {
                active(i) = 0;
                continue;
            } else {
                //check that corners are in domain
                for (int c = 0; c < cpmp; c++) {
                    //corner positions are stored over time
                    tmpVec = corner_positions(i * cpmp + c);

                    if ((use_elem && !job->grid->inDomain(job, tmpVec, elem(i, c))) ||
                        !job->grid->inDomain(job, tmpVec)) {
                        //corner out of domain
                        ith_cpdi = DefaultBody::CPDI_OFF;
                        break;
                    }
                }

                //initialize std::vectors
                nvec.resize(0);
                valvec.resize(0);

                //ony update corners if in domain
                if (ith_cpdi == DefaultBody::CPDI_ON) {
                    //store corner positions
                    x1 = corner_positions(i * cpmp + 0, 0);
                    x2 = corner_positions(i * cpmp + 1, 0);
                    x3 = corner_positions(i * cpmp + 3, 0);
                    x4 = corner_positions(i * cpmp + 2, 0);

                    y1 = corner_positions(i * cpmp + 0, 1);
                    y2 = corner_positions(i * cpmp + 1, 1);
                    y3 = corner_positions(i * cpmp + 3, 1);
                    y4 = corner_positions(i * cpmp + 2, 1);

                    //calculate a and b constants
                    a1 = 0.25 * (-x1 + x2 + x3 - x4);
                    a2 = 0.25 * (-x1 - x2 + x3 + x4);
                    a3 = 0.25 * (x1 - x2 + x3 - x4);

                    b1 = 0.25 * (-y1 + y2 + y3 - y4);
                    b2 = 0.25 * (-y1 - y2 + y3 + y4);
                    b3 = 0.25 * (y1 - y2 + y3 - y4);

                    //calculate volume
                    vp = 4.0 * (a1 * b2 - a2 * b1);

                    if (vp <= 0) {
                        std::cout << "WARNING: Point has zero volume!." << std::endl;
                    }

                    //determine corner weights
                    for (int c = 0; c < A.rows(); c++) {
                        //initialize storage vectors
                        nvec.resize(0);
                        valvec.resize(0);

                        //determine corner positions
                        tmpVec = corner_positions(i * cpmp + c);

                        if (use_elem) {
                            job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec, elem(i, c));
                        } else {
                            job->grid->evaluateBasisFnValue(job, tmpVec, nvec, valvec);
                        }

                        //corner weighting function
                        double w = a1 * b2 - a2 * b1 +
                                   ((a1 * b3 - a3 * b1) * A(c, 0) + (a3 * b2 - a2 * b3) * A(c, 1)) / 3.0;
                        //adjust weight by volume of point
                        w /= vp;

                        //use nodal velocities to update corner position
                        for (int j = 0; j < nvec.size(); j++) {
                            corner_positions[i*cpmp + c] += job->dt*body->nodes->x_t[nvec[j]]*valvec[j];
                        }
                    }

                    //update corner positions
                    x1 = corner_positions(i * cpmp + 0, 0);
                    x2 = corner_positions(i * cpmp + 1, 0);
                    x3 = corner_positions(i * cpmp + 3, 0);
                    x4 = corner_positions(i * cpmp + 2, 0);

                    y1 = corner_positions(i * cpmp + 0, 1);
                    y2 = corner_positions(i * cpmp + 1, 1);
                    y3 = corner_positions(i * cpmp + 3, 1);
                    y4 = corner_positions(i * cpmp + 2, 1);

                    //update a and b constants
                    a1 = 0.25 * (-x1 + x2 + x3 - x4);
                    a2 = 0.25 * (-x1 - x2 + x3 + x4);
                    a3 = 0.25 * (x1 - x2 + x3 - x4);

                    b1 = 0.25 * (-y1 + y2 + y3 - y4);
                    b2 = 0.25 * (-y1 - y2 + y3 + y4);
                    b3 = 0.25 * (y1 - y2 + y3 - y4);

                    //update volume
                    vp = 4.0 * (a1 * b2 - a2 * b1);
                    v[i] = vp;

                    //zero centroid before update
                    x_i.setZero();

                    //update centroid position
                    for (int c=0; c<A.rows(); c++){
                        //corner weighting function
                        double w = a1 * b2 - a2 * b1 +
                                   ((a1 * b3 - a3 * b1) * A(c, 0) + (a3 * b2 - a2 * b3) * A(c, 1)) / 3.0;
                        //adjust weight by volume of point
                        w /= vp;

                        //add contribution to centroid
                        x_i += corner_positions[i*cpmp + c] * w;
                    }
                }
            }
        }
    }

    //apply position correction
    if (POSITIONRULE == NO_CORRECTION){
        //do nothing
    } else if (POSITIONRULE == AVAV) {
        //do something!

        //iterate counter
        skip_counter += 1;

        //check if counter % skip_value = 0
        if (skip_counter % skip_value == 0){
            //run avoid-a-void algorithms
            //step -1: create grid based b, u, x_t, T fields
            Eigen::VectorXd v_nodes = body->S*v;       //apparent volume at node positions
            KinematicVectorArray vb_points = b;
            KinematicVectorArray vu_points = u;
            MaterialTensorArray vT_points = T;
            for (int p=0; p<x.size(); p++){
                vb_points[p] *= v(p);
                vu_points[p] *= v(p);
                vT_points[p] *= v(p);
            }
            KinematicVectorArray b_nodes = body->S*vb_points; //integrated b field at node
            KinematicVectorArray u_nodes = body->S*vu_points; //integrated u field at node
            MaterialTensorArray T_nodes = body->S*vT_points;  //integrated T field at node
            for (int i=0; i<body->nodes->x.size(); i++){
                b_nodes[i] /= v_nodes(i);
                u_nodes[i] /= v_nodes(i);
                T_nodes[i] /= v_nodes(i);
            }

            //step 0: create cell/point map
            for (int i=0; i<cell_to_point_map.size(); i++){
                cell_to_point_map[i].clear(); //clear previous map
            }
            int tmp_id = -1;
            for (int p=0; p<x.size(); p++){
                if (active(p) == 0){
                    point_to_cell_map[p] = -1; //do not track this point
                } else {
                    tmp_id = pos_to_cell(job, x[p]);
                    point_to_cell_map[p] = tmp_id;          //add cell to point list
                    if (tmp_id >= 0 && tmp_id < cell_to_point_map.size()) {
                        cell_to_point_map[tmp_id].push_back(p); //add point to cell list
                    }
                }
            }

            //step 1: create signed dist. field
            KinematicVectorArray avavS = KinematicVectorArray(0, job->JOB_TYPE); //set of zero crossings
            std::vector<int> ijk, rst, iijjkk;
            double dist = 0;
            for (int i=0; i<job->grid->node_count; i++){
                d(i) = alg_inf; //initialize with p.d = inf

                ijk = cell_to_ijk(job, node_to_cell_map[i]); //get ijk position of node in search grid
                rst = ijk;
                iijjkk = std::vector<int>(ijk.size(),-1); //initialize offset counter

                for (int cell = 0; cell<number_of_neighbors; cell++){
                    //determine direction of next search cell
                    for (int dir = 0; dir<ijk.size(); dir++){
                        if (iijjkk[dir] < 1){
                            iijjkk[dir] += 1;
                            break;
                        } else {
                            iijjkk[dir] = -1;
                        }
                    }

                    //ijk position of search cell
                    for (int dir=0; dir<ijk.size(); dir++){
                        rst[dir] = ijk[dir]+iijjkk[dir];
                    }

                    //search cell id
                    tmp_id = ijk_to_cell(job, rst);

                    if (tmp_id >= 0 && tmp_id < cell_to_point_map.size()) {
                        //find min dist. of mpm points
                        for (int p = 0; p < cell_to_point_map[tmp_id].size(); p++) {
                            dist = (body->nodes->x[i] - body->points->x[cell_to_point_map[tmp_id][p]]).norm() - r;
                            if (dist < d(i)) {
                                d(i) = dist;
                            }
                        }
                    }
                }
            }

            int n0, n1;
            double a;
            KinematicVector zero_crossing_pos(job->JOB_TYPE);
            for (int e=0; e<(edge_list.size()/2); e++){
                n0 = edge_list[2*e];
                n1 = edge_list[2*e+1];
                if (d(n0)*d(n1) < 0){
                    //if dist fields have opposite signs, then find zero crossing
                    a = d(n1)/(d(n1) - d(n0));

                    //calculate position
                    zero_crossing_pos = a*body->nodes->x[n0] + (1.0-a)*body->nodes->x[n1];

                    //add to list
                    avavS.push_back(zero_crossing_pos);
                }
            }

            for (int i=0; i<job->grid->node_count; i++){
                d(i) = ((d(i) > 0) - (d(i) < 0))*alg_inf;
            }

            for (int c=0; c<avavS.size(); c++){
                tmp_id = pos_to_cell(job, avavS[c]); //find search cell of grid point

                ijk = cell_to_ijk(job, tmp_id); //get ijk position of zero crossing
                rst = ijk;
                iijjkk = std::vector<int>(ijk.size(),-2); //initialize offset counter

                for (int cell = 0; cell<number_of_second_neighbors; cell++){
                    //determine direction of next search cell
                    for (int dir = 0; dir<ijk.size(); dir++){
                        if (iijjkk[dir] < 2){
                            iijjkk[dir] += 1;
                            break;
                        } else {
                            iijjkk[dir] = -2;
                        }
                    }

                    //ijk position of search cell
                    for (int dir=0; dir<ijk.size(); dir++){
                        rst[dir] = ijk[dir]+iijjkk[dir];
                    }

                    //search cell id (overwrites tmp_id in previous scope, so be careful)
                    tmp_id = ijk_to_cell(job, rst);

                    if (tmp_id >= 0 && tmp_id < cell_to_node_map.size()){
                        //find min dist. of mpm points
                        for (int n=0; n<cell_to_node_map[tmp_id].size(); n++){
                            n0 = cell_to_node_map[tmp_id][n];
                            dist = (body->nodes->x[n0] - avavS[c]).norm();
                            if (dist < std::abs(d(n0))){
                                d(n0) = ((d(n0) > 0) - (d(n0) < 0))*dist;
                            }
                        }
                    }
                }
            }

            //step 2: add more points to interior (if they fit)
            //use search cell grid for insertion of points
            KinematicVector x_new = KinematicVector(job->JOB_TYPE);
            std::vector<int> n_list = std::vector<int>(0);
            std::vector<double> s_list = std::vector<double>(0);
            std::vector<int> neighbor_point_indices;
            int p0, p1;
            double theta, z, l, p_d, s_sum;
            KinematicVector v_sum;
            MaterialTensor t_sum;
            bool accept_point;
            for (int c=0; c<cell_to_point_map.size(); c++){
                if (cell_to_point_map[c].size() > 0){
                    //pick a random point from point to cell map
                    p1 = cell_to_point_map[c][0];                   //zero is 'random' enough to test the alg.
                    x_new = x(p1);
                    if (job->grid->GRID_DIM == 1){
                        std::cerr << "UNSURE HOW TO PROCEED WITH AVOID-A-VOID ALGORITHM! Exiting." << std::endl;
                    } else if (job->grid->GRID_DIM == 2){
                        theta = 2.0*M_PI*std::rand()/RAND_MAX; //random angle
                        x_new[0] += r*(std::sqrt(3.0)/2.0 + 0.01)*std::cos(theta);
                        x_new[1] += r*(std::sqrt(3.0)/2.0 + 0.01)*std::sin(theta);
                    } else {
                        z = 1.0 - 2.0*std::rand()/RAND_MAX;     //random height
                        theta = 2.0*M_PI*std::rand()/RAND_MAX;  //random angle
                        l = (1.0 - z*z);
                        x_new[0] += l*r*(std::sqrt(3.0)/2.0 + 0.01)*std::cos(theta);
                        x_new[1] += l*r*(std::sqrt(3.0)/2.0 + 0.01)*std::sin(theta);
                        x_new[2] += z*r*(std::sqrt(3.0)/2.0 + 0.01);
                    }
                } else {
                    //pick cell center
                    ijk = cell_to_ijk(job, c);
                    for (int dir=0; dir<ijk.size(); dir++){
                        x_new[dir] = (ijk[dir] + 0.5)*hx[dir];
                    }
                }

                //check that point is in domain
                accept_point = true;
                if (!job->grid->inDomain(job, x_new)){
                    accept_point = false;
                }

                //check that point is in domain
                if (accept_point) {
                    n_list.clear();
                    s_list.clear();
                    p_d = 0;
                    job->grid->evaluateBasisFnValue(job, x_new, n_list, s_list); //get shape function values
                    for (int ii = 0; ii < n_list.size(); ii++) {
                        if (n_list[ii] > -1) {
                            p_d += d(n_list[ii])*s_list[ii];
                        }
                    }

                    if (p_d >= -2.2*r){
                        //point isn't interior enough
                        accept_point = false;
                    }
                }

                if (accept_point) {
                    //check if point is within \alpha*r of any other points
                    ijk = cell_to_ijk(job, c); //get ijk position of search cell
                    rst = ijk;
                    iijjkk = std::vector<int>(ijk.size(), -1); //initialize offset counter

                    for (int cell = 0; cell < number_of_neighbors; cell++) {
                        //determine direction of next search cell
                        for (int dir = 0; dir < ijk.size(); dir++) {
                            if (iijjkk[dir] < 1) {
                                iijjkk[dir] += 1;
                                break;
                            } else {
                                iijjkk[dir] = -1;
                            }
                        }

                        //ijk position of search cell
                        for (int dir = 0; dir < ijk.size(); dir++) {
                            rst[dir] = ijk[dir] + iijjkk[dir];
                        }

                        //search cell id
                        tmp_id = ijk_to_cell(job, rst);

                        if (tmp_id >= 0 && tmp_id < cell_to_point_map.size()) {
                            //find min dist. of mpm points to new point
                            for (int p = 0; p < cell_to_point_map[tmp_id].size(); p++) {
                                dist = (x_new - body->points->x[cell_to_point_map[tmp_id][p]]).norm();
                                if (dist < (std::sqrt(3.0) / 2.0 + 0.01) * r) {
                                    accept_point = false;
                                    break;
                                }
                            }
                        }

                        //if accept_point false, not need to continue search
                        if (!accept_point) {
                            break;
                        }
                    }
                }

                //if new point is acceptable and there are points in the buffer, we need to create it
                if (accept_point && !buffer_list.empty()){
                    //find nearest neighbor points
                    neighbor_point_indices.clear();
                    for (int cell = 0; cell < number_of_neighbors; cell++) {
                        //determine direction of next search cell
                        for (int dir = 0; dir < ijk.size(); dir++) {
                            if (iijjkk[dir] < 1) {
                                iijjkk[dir] += 1;
                                break;
                            } else {
                                iijjkk[dir] = -1;
                            }
                        }

                        //ijk position of search cell
                        for (int dir = 0; dir < ijk.size(); dir++) {
                            rst[dir] = ijk[dir] + iijjkk[dir];
                        }

                        //search cell id
                        tmp_id = ijk_to_cell(job, rst);

                        //add points to neighbor list
                        if (tmp_id >= 0 && tmp_id < cell_to_point_map.size()) {
                            for (int p = 0; p < cell_to_point_map[tmp_id].size(); p++) {
                                neighbor_point_indices.push_back(cell_to_point_map[tmp_id][p]);
                            }
                        }
                    }

                    //get new point index from buffer
                    tmp_id = buffer_list.back();
                    //delete index from buffer
                    buffer_list.pop_back();

                    //initialize point
                    point_to_cell_map[tmp_id] = c;
                    cell_to_point_map[c].push_back(tmp_id);
                    x[tmp_id] = x_new;
                    active(tmp_id) = 1;

                    //assign values as weighted sum of neighbors
                    //extent is averaged
                    s_sum = 0;
                    for (int p=0; p<neighbor_point_indices.size(); p++){
                        s_sum += extent(neighbor_point_indices[p]);
                    }
                    extent(tmp_id) = s_sum/neighbor_point_indices.size();

                    //volume is conserved
                    s_sum = 0;
                    for (int p=0; p<neighbor_point_indices.size(); p++){
                        s_sum += v(neighbor_point_indices[p]);
                        v(neighbor_point_indices[p]) -= v(neighbor_point_indices[p])/(neighbor_point_indices.size()+1);
                    }
                    v(tmp_id) = s_sum/(neighbor_point_indices.size()+1);

                    //mass is conserved
                    s_sum = 0;
                    for (int p=0; p<neighbor_point_indices.size(); p++){
                        s_sum += m(neighbor_point_indices[p]);
                        m(neighbor_point_indices[p]) -= m(neighbor_point_indices[p])/(neighbor_point_indices.size()+1);
                    }
                    m(tmp_id) = s_sum/(neighbor_point_indices.size()+1);

                    //initial volume is conserved
                    s_sum = 0;
                    for (int p=0; p<neighbor_point_indices.size(); p++){
                        s_sum += v0(neighbor_point_indices[p]);
                        v0(neighbor_point_indices[p]) -= v0(neighbor_point_indices[p])/(neighbor_point_indices.size()+1);
                    }
                    v0(tmp_id) = s_sum/(neighbor_point_indices.size()+1);


                    //b, u, x_t, mx_t, T mapped from grid
                    n_list.clear();
                    s_list.clear();
                    b[tmp_id].setZero();
                    u[tmp_id].setZero();
                    x_t[tmp_id].setZero();
                    T[tmp_id].setZero();
                    job->grid->evaluateBasisFnValue(job, x_new, n_list, s_list); //get shape function values
                    for (int ii = 0; ii < n_list.size(); ii++) {
                        if (n_list[ii] > -1) {
                            b[tmp_id] += b_nodes[n_list[ii]]*s_list[ii];
                            u[tmp_id] += u_nodes[n_list[ii]]*s_list[ii];
                            x_t[tmp_id] += body->nodes->x_t[n_list[ii]]*s_list[ii];
                            T[tmp_id] += T_nodes[n_list[ii]]*s_list[ii];
                        }
                    }
                    mx_t[tmp_id] = m(tmp_id)*x_t[tmp_id];

                    if (p_d >= -2.2*r){
                        //point isn't interior enough
                        accept_point = false;
                    }

                } else if (buffer_list.empty()){
                    //no more buffer points
                    std::cout << "Buffer list empty. No more points pre-allocated. May need to implement dynamic allocation" << std::endl;
                    break; //exit creation loop
                }
            }

            //step 3: determine points which should be merged
            int p_other = -1;
            int c = -1;
            bool remove_point = false;
            for (int p=0; p<x.size(); p++){
                if (active(p) != 0) {
                    remove_point = false;
                    c = point_to_cell_map[p];
                    if (c < 0) {
                        //uh oh
                        std::cerr << "[" << p << "] is active, but hasn't identified a search cell. Exiting." << std::endl;
                        exit(0);
                    }

                    //check if point is within 0.03*r of any other points
                    ijk = cell_to_ijk(job, c); //get ijk position of search cell
                    rst = ijk;
                    iijjkk = std::vector<int>(ijk.size(), -1); //initialize offset counter

                    for (int cell = 0; cell < number_of_neighbors; cell++) {
                        //determine direction of next search cell
                        for (int dir = 0; dir < ijk.size(); dir++) {
                            if (iijjkk[dir] < 1) {
                                iijjkk[dir] += 1;
                                break;
                            } else {
                                iijjkk[dir] = -1;
                            }
                        }

                        //ijk position of search cell
                        for (int dir = 0; dir < ijk.size(); dir++) {
                            rst[dir] = ijk[dir] + iijjkk[dir];
                        }

                        //search cell id
                        tmp_id = ijk_to_cell(job, rst);

                        //check that tmp_id is valid
                        if (tmp_id < 0 || tmp_id >= cell_to_point_map.size()){
                            continue; //skip this cell
                        }

                        //find min dist. of mpm points to new point
                        for (int q = 0; q < cell_to_point_map[tmp_id].size(); q++) {
                            dist = (x[p] - x[cell_to_point_map[tmp_id][q]]).norm();
                            if (dist < merge_dist && cell_to_point_map[tmp_id][q] != p) {
                                p_other = cell_to_point_map[tmp_id][q];
                                remove_point = true;
                                break;
                            }
                        }

                        //if remove_point true, do not need to continue search
                        if (remove_point) {
                            break;
                        }
                    }

                    if (remove_point){
                        //need to merge point with neighbor
                        m(p_other) += m(p);
                        v(p_other) += v(p);
                        v0(p_other) += v0(p);
                        mx_t[p_other] += mx_t[p];
                        x_t[p_other] = mx_t[p_other]/m(p_other);
                        T[p_other] = ((v(p_other) - v(p))*T[p_other] + v(p)*T[p])/v(p_other);

                        //need to remove myself from simulation
                        active(p) = 0;
                        point_to_cell_map[p] = -1;

                        //need to remove myself from cell_to_point map
                        for (int pp=0; pp<cell_to_point_map[c].size(); pp++){
                            if (cell_to_point_map[c][pp] == p){
                                cell_to_point_map[c].erase(cell_to_point_map[c].begin()+pp);
                            }
                        }

                        //need to add myself to buffer
                        buffer_list.push_back(p);
                    }
                }
            }
        }
    } else if (POSITIONRULE == SPH_LIKE) {
        KinematicVectorArray delta = KinematicVectorArray(x.size(), x.VECTOR_TYPE);

        //iterate counter
        skip_counter += 1;

        //check if counter % skip_value = 0
        if (skip_counter % skip_value == 0) {
            //run sph-like algorithms
            //step -1: find max velocity
            u_max = 0;
            for (int p=0; p<x.size(); p++){
                if (x_t[p].norm() > u_max){
                    u_max = x_t[p].norm();
                }
            }

            //step 0: create cell/point map
            for (int i = 0; i < cell_to_point_map.size(); i++) {
                cell_to_point_map[i].clear(); //clear previous map
            }
            int tmp_id = -1;
            for (int p = 0; p < x.size(); p++) {
                if (active(p) == 0) {
                    point_to_cell_map[p] = -1; //do not track this point
                } else {
                    tmp_id = pos_to_cell(job, x[p]);
                    point_to_cell_map[p] = tmp_id;          //add cell to point list
                    if (tmp_id >= 0 && tmp_id < cell_to_point_map.size()) {
                        cell_to_point_map[tmp_id].push_back(p); //add point to cell list
                    }
                }
            }

            //step 1: calculate point-wise position correction
            int p_other = -1;
            int c = -1;
            std::vector<int> ijk, rst, iijjkk;
            double dist, r_bar_i;
            KinematicVector R_i = KinematicVector(job->JOB_TYPE);
            int M_i;
            std::vector<int> point_neighbor_list = std::vector<int>(0);

            for (int p=0; p<x.size(); p++){
                if (active(p) != 0) {
                    //zero out point neighbors list
                    point_neighbor_list.clear();

                    c = point_to_cell_map[p];
                    if (c < 0) {
                        //uh oh
                        std::cerr << "[" << p << "] is active, but hasn't identified a search cell. Exiting." << std::endl;
                        exit(0);
                    }

                    ijk = cell_to_ijk(job, c); //get ijk position of search cell
                    rst = ijk;
                    iijjkk = std::vector<int>(ijk.size(), -1); //initialize offset counter

                    for (int cell = 0; cell < number_of_neighbors; cell++) {
                        //determine direction of next search cell
                        for (int dir = 0; dir < ijk.size(); dir++) {
                            if (iijjkk[dir] < 1) {
                                iijjkk[dir] += 1;
                                break;
                            } else {
                                iijjkk[dir] = -1;
                            }
                        }

                        //ijk position of search cell
                        for (int dir = 0; dir < ijk.size(); dir++) {
                            rst[dir] = ijk[dir] + iijjkk[dir];
                        }

                        //search cell id
                        tmp_id = ijk_to_cell(job, rst);

                        //check that tmp_id is valid
                        if (tmp_id < 0 || tmp_id >= cell_to_point_map.size()){
                            continue; //skip this cell
                        }

                        //add neighbors to list
                        for (int q = 0; q < cell_to_point_map[tmp_id].size(); q++) {
                            dist = (x[p] - x[cell_to_point_map[tmp_id][q]]).norm();
                            if (cell_to_point_map[tmp_id][q] != p) {
                                if (dist < 2.0*r){
                                    //add point to neighbor list
                                    point_neighbor_list.push_back(cell_to_point_map[tmp_id][q]);
                                }
                            }
                        }
                    }

                    //calculate M_i, r_bar_i
                    M_i = point_neighbor_list.size();
                    if (M_i > 0){
                        //calculate delta
                        r_bar_i = 0;
                        for (int q=0; q<M_i; q++){
                            r_bar_i += (x[p] - x[point_neighbor_list[q]]).norm()/M_i;
                        }

                        R_i.setZero();
                        for (int q=0; q<M_i; q++){
                            dist = (x[p] - x[point_neighbor_list[q]]).norm();
                            R_i += r_bar_i*r_bar_i*(x[p] - x[point_neighbor_list[q]])/(dist*dist*dist);
                        }

                        delta[p] = sph_const*u_max*skip_value*job->dt*R_i;

                    } else {
                        //no position correction
                        delta[p].setZero();
                    }
                }
            }

            //step 2: adjust positions
            for (int p=0; p<x.size(); p++){
                body->points->x(p) += delta(p);
                body->points->u(p) += delta(p);
                del_pos(p) += delta(p);

                body->points->x_t(p) += body->points->L(p) * delta(p);
                body->points->mx_t(p) = body->points->m(p) * body->points->x_t(p);
            }
        }
    } else if (POSITIONRULE == DELTA_STRAIN){
        KinematicVector delta = KinematicVector(x.VECTOR_TYPE);
        for (int i=0; i<x.size(); i++) {
            double trD = L[i].trace();
            //delta = -alpha * job->dt * (L[i] - (trD / 3.0) * KinematicTensor::Identity(x.VECTOR_TYPE)).norm() * grad_e(i) * h * h;
            delta = -alpha * job->dt * v(i) * (L[i] - (trD / 3.0) * KinematicTensor::Identity(x.VECTOR_TYPE)).norm() * grad_eonV(i) * h * h;

            for (int pos = 0; pos < delta.rows(); pos++) {
                if (//((body->points->x(i, pos) + delta(pos)) >= Lx(pos)) ||
                    //((body->points->x(i, pos) + delta(pos)) <= 0) ||
                        !std::isfinite(delta(pos))) {
                    delta(pos) = 0;
                }
            }

            body->points->x(i) += delta;
            body->points->u(i) += delta;
            del_pos(i) += delta;

            body->points->x_t(i) += body->points->L(i) * delta;
            body->points->mx_t(i) = body->points->m(i) * body->points->x_t(i);
        }
    } else if (POSITIONRULE == DELTA_DISP){
        KinematicVector delta = KinematicVector(x.VECTOR_TYPE);
        for (int i=0; i<x.size(); i++) {
            delta = -alpha * job->dt * x_t[i].norm() * v(i) * grad_eonV(i) * h;

            for (int pos = 0; pos < delta.rows(); pos++) {
                if (//((body->points->x(i, pos) + delta(pos)) >= Lx(pos)) ||
                    //((body->points->x(i, pos) + delta(pos)) <= 0) ||
                        !std::isfinite(delta(pos))) {
                    delta(pos) = 0;
                }
            }

            body->points->x(i) += delta;
            body->points->u(i) += delta;
            del_pos(i) += delta;

            body->points->x_t(i) += body->points->L(i) * delta;
            body->points->mx_t(i) = body->points->m(i) * body->points->x_t(i);
        }
    } else if (POSITIONRULE == DELTA_NEWTON){
        //iteration counter for newton solve
        int iter_count = 0;

        //error measure for newton solve
        double relative_error = 0;

        //position update vector
        KinematicVector delta = KinematicVector(x.VECTOR_TYPE);

        //calculate relative error measure (||H/V_i||_\infty)
        for (int i=0; i<body->nodes->x.size(); i++){
            if (e(i) > relative_error){
                relative_error = e(i);
            }
        }

        //newton scheme
        double step_length;
        while (iter_count <= max_iter && relative_error > REL_TOL){
            //increment iter_count
            iter_count += 1;

            //calculate step length from pg 88 nb #8
            step_length = 0;
            for (int p=0; p<x.size(); p++){
                step_length += v(p) * v(p) * grad_H[p].dot(grad_H[p]);
            }
            step_length = H.squaredNorm()/(2.0*step_length);

            //make sure step_length is not nan
            if (!std::isfinite(step_length)){
                //do nothing
                break;
            }

            //apply limiter to make gradient descent
            step_length *= alpha;

            //step in that direction
            for (int i=0; i<x.size(); i++) {
                delta = -step_length * v(i) * grad_H(i);

                if (delta.norm() > h){
                    delta *= h/delta.norm();
                }

                //std::cout << delta.norm() << ", " << v(i) << ", " << grad_H(i).norm() << std::endl;

                for (int pos = 0; pos < delta.rows(); pos++) {
                    if (//((body->points->x(i, pos) + delta(pos)) >= Lx(pos)) ||
                        //((body->points->x(i, pos) + delta(pos)) <= 0) ||
                            !std::isfinite(delta(pos))) {
                        delta(pos) = 0;
                    }
                }

                body->points->x(i) += delta;
                body->points->u(i) += delta;
                del_pos(i) += delta;

                body->points->x_t(i) += body->points->L(i) * delta;
                body->points->mx_t(i) = body->points->m(i) * body->points->x_t(i);
            }

            //update mappings
            generateMap(job, body, 0);

            //calculate strain rate
            L = body->gradS.tensor_product_transpose(body->nodes->x_t, MPMSparseMatrixBase::TRANSPOSED);

            //calculate relative error measure (||H/V_i||_\infty)
            relative_error = 0;
            for (int i=0; i<body->nodes->x.size(); i++){
                if (e(i) > relative_error){
                    relative_error = e(i);
                }
            }

            /*
            std::cout << "t: " << job->t << ", i: " << iter_count;
            std::cout << ", step_length: " << step_length << ", err: " << relative_error;
            std::cout << ", H: " << H.norm() << std::endl;
             */
        }
    }

    return;
}

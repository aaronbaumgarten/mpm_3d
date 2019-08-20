//
// Created by aaron on 7/31/19.
// improved_quadrature_points.cpp
//

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
#include "parser.hpp"

#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

#include "job.hpp"

#include "points.hpp"
#include "objects/bodies/bodies.hpp"

/* input to this file denotes the type of improved quadrature to use
     * 0 -- No quadrature improvement
     * 1 -- uGIMP
     * 2 -- cpGIMP
     * 3 -- CPDI
     * 4 -- CPDI2
     * 5 -- Avoid-a-void
     * 6 -- \delta-correction, constant strain
     * 7 -- \delta-correction, constant time
     */

/*----------------------------------------------------------------------------*/
//initialize point state (assumes that readFromFile has been called)
//no safety check on this, so be careful please
void ImprovedQuadraturePoints::init(Job* job, Body* body){
    if (int_props.size() < 1){
        std::cout << int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 1 property defined (QUADRULE).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //assign quadrature rule
        QUADRULE = int_props[0];
    }

    //all points require v0 initialization
    v0 = v;

    //initialize error measures
    V_i.resize(job->grid->node_count,1);
    v_i.resize(job->grid->node_count,1);
    e.resize(job->grid->node_count,1);
    grad_e = KinematicVectorArray(body->points->x.size(),job->JOB_TYPE);

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
        if (fp64_props.size() > 1 && QUADRULE == UGIMP) {
            scale_factor = fp64_props[0];
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
        corner_positions = KinematicVectorArray(x.size()*cp,x.VECTOR_TYPE);
        for (int i = 0; i < x.size(); i++) {
            for (int c = 0; c < A.rows(); c++) {
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++) {
                    corner_positions[i,c*x.DIM + pos] = x[i,pos] + 0.5 * extent(i) * A(c, pos);
                }
                for (int pos = job->grid->GRID_DIM; pos < x.DIM; pos++) {
                    corner_positions[i,c*x.DIM + pos] = 0;
                }
            }
        }
    }

    //delta correction for constant strain requires h and alpha
    if (QUADRULE == DELTA_CS){
        if (fp64_props.size() < 2){
            std::cout << fp64_props.size() << "\n";
            fprintf(stderr,
                    "%s:%s: Need at least 2 properties defined (h, alpha).\n",
                    __FILE__, __func__);
            exit(0);
        } else {
            //assign grid scale and rate scale
            h = fp64_props[0];
            alpha = fp64_props[1];
        }

        //initialize position correction measure
        del_pos = KinematicVectorArray(body->points->x.size(),job->JOB_TYPE);
        del_pos.setZero();
    }

    //Avoid-a-void not implemented yet
    if (QUADRULE != STANDARD ||
        QUADRULE != UGIMP ||
        QUADRULE != CPGIMP ||
        QUADRULE != CPDI ||
        QUADRULE != CPDI2 ||
        QUADRULE != DELTA_CS
            ){
        std::cerr << "ERROR: Rule " << QUADRULE << " not implemented in ImprovedQuadraturePoints yet. Exiting." << std::endl;
        exit(0);
    }

    if (int_props.size() < 2){
        //standard behavior
        use_elem = false;
    } else if (int_props[1] == 1){
        //create point-wise element list
        use_elem = true;
        elem = Eigen::MatrixXi(x.size(),cp);
        elem.setConstant(-1);
        std::cout << object_name << " using point-element history." << std::endl;
    }

    std::cout << "Points Initialized: [" << file << "]. Using Quadrature Rule " << QUADRULE << "." << std::endl;
    std::cout << " * 0 -- standard MPM\n";
    std::cout << " * 1 -- uGIMP\n";
    std::cout << " * 2 -- cpGIMP\n";
    std::cout << " * 3 -- CPDI\n";
    std::cout << " * 4 -- CPDI2\n";
    std::cout << " * 5 -- Avoid-a-void\n";
    std::cout << " * 6 -- \\delta-correction, constant strain" << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//read point data from file
void DefaultPoints::readFromFile(Job *job, Body *body, std::string fileIN) {
    //first line lists number of points in file
    //subsequent lines contain m, v, x, x_t, active
    file = fileIN;

    std::string line;
    std::ifstream fin(file);
    std::stringstream ss;
    std::vector<std::string> s_vec;

    if (fin.is_open()) {
        std::getline(fin, line);
        int len = std::stoi(line);

        //size KinematicVectors
        x = KinematicVectorArray(len, job->JOB_TYPE);
        u = KinematicVectorArray(len, job->JOB_TYPE);
        x_t = KinematicVectorArray(len, job->JOB_TYPE);
        mx_t = KinematicVectorArray(len, job->JOB_TYPE);
        b = KinematicVectorArray(len, job->JOB_TYPE);

        //size scalar vectors
        m.resize(len);
        v.resize(len);
        v0.resize(len);
        active.resize(len);
        extent.resize(len);

        //size tensor arrays
        T = MaterialTensorArray(len);
        L = KinematicTensorArray(len, job->JOB_TYPE);

        //zero out all entries to start
        x.setZero();
        u.setZero();
        x_t.setZero();
        m.setZero();
        v.setZero();
        v0.setZero();
        mx_t.setZero();
        b.setZero();
        T.setZero();
        L.setZero();
        active.setZero();
        extent.setZero();

        for (int i = 0; i < len; i++) {
            std::getline(fin, line);
            s_vec = Parser::splitString(line, ' ');
            if (s_vec.size() == (1 + 1 + job->DIM + job->DIM + 1)) {
                m[i] = std::stod(s_vec[0]);                 //first column is mass
                v[i] = std::stod(s_vec[1]);                 //second column is volume
                for (int pos = 0; pos < job->DIM; pos++) {
                    x(i, pos) = std::stod(s_vec[2 + pos]);    //following cols are position
                }
                for (int pos = 0; pos < job->DIM; pos++) {
                    x_t(i, pos) = std::stod(s_vec[2 + job->DIM + pos]);     //following cols are velocity
                }
                active[i] = std::stod(s_vec[2 + 2 * job->DIM]);

            } else if (s_vec.size() == (1 + 1 + job->grid->GRID_DIM + job->grid->GRID_DIM + 1)) {
                m[i] = std::stod(s_vec[0]);                 //first column is mass
                v[i] = std::stod(s_vec[1]);                 //second column is volume
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++) {
                    x(i, pos) = std::stod(s_vec[2 + pos]);    //following cols are position
                }
                for (int pos = 0; pos < job->grid->GRID_DIM; pos++) {
                    x_t(i, pos) = std::stod(s_vec[2 + job->grid->GRID_DIM + pos]);     //following cols are velocity
                }
                active[i] = std::stod(s_vec[2 + 2 * job->grid->GRID_DIM]);

            } else {
                std::cerr << "ERROR: Unable to read line: " << file << std::endl;
                return;
            }

            //correct volume and mass for axisymmetric simulation
            if (job->JOB_TYPE == job->JOB_AXISYM){
                m[i] *= x(i,0);
                v[i] *= x(i,0);
            }
        }

    } else {
        std::cerr << "ERROR: Unable to open file: " << file << std::endl;
        return;
    }

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

    if (QUADRULE == STANDARD || QUADRULE == DELTA_CS) {
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
                    for (int pos = 0; pos < job->grid->GRID_DIM; pos++){
                        tmpVec[pos] = x_i[pos] + 0.5*extent(i)*A(c,pos);
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
                                for (int ii=0; ii<job->grid->GRID_DIM; ii){
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
                    double detF = 1;
                    if (QUADRULE == UGIMP) {
                        detF = 1;
                    } else if (QUADRULE == CPGIMP){
                        for (int ii=0; ii<job->grid->GRID_DIM; ii++){
                            detF *= F(i,ii,ii);
                        }
                    } else if (QUADRULE == CPDI){
                        detF = F(i).det();
                    }

                    for (int j = 0; j < nvec.size(); j++) {
                        body->S.push_back(nvec[j], i, detF * valvec[j] / A.rows()); //node, point, value
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
                                for (int ii=0; ii<job->grid->GRID_DIM; ii){
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

                    /* CALCULATE S_ip */
                    //initialize storage vectors
                    nvec.resize(0);
                    valvec.resize(0);
                    gradvec.resize(0);

                    //determine corner weights
                    for (int c = 0; c < A.rows(); c++) {
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
                            body->S.push_back(nvec[j], i, w * valvec[j] / A.rows()); //node, point, value

                            //gradient contribution from corner
                            tmpGrad[0] = (A(c,0)*b2 - A(c,1)*(a2 + a3*A(c,0)/3.0))/vp;
                            tmpGrad[1] = (-A(c,0)*(b1 + b3*A(c,1)/3.0) + A(c,1)*a1)/vp;
                            for (int pos = job->grid->GRID_DIM; pos < tmpGrad.DIM; pos++) {
                                tmpGrad[pos] = 0;
                            }

                            //add to point gradient
                            gradvec.push_back(tmpGrad);
                        }
                    }
                }
            }
        }
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

    double tmpNum
    for (int i = 0; i < V_i.rows(); i++) {
        tmpNum = (v_i(i) - V_i(i)) / (V_i(i));
        e(i) = std::max(0.0, tmpNum);
    }

    //calculate gradient of nodal overshoot
    grad_e = body->gradS.operate(e, MPMSparseMatrixBase::TRANSPOSED);

    return;
}


/*----------------------------------------------------------------------------*/
//write relavent point data to file
void ImprovedQuadraturePoints::writeFrame(Job* job, Body* body, Serializer* serializer){
    //call default function
    DefaultPoints::writeFrame(job, body, serializer);

    //custom file output

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
            F[i] = tmpL.exp(4)*F[i];

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

                    //initialize storage vectors
                    nvec.resize(0);
                    valvec.resize(0);

                    //determine corner weights
                    for (int c = 0; c < A.rows(); c++) {
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
                            corner_positions[i*cpmp + c] += w*job->dt*body->nodes->x_t[nvec[j]]*valvec[j];
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
    } else if (QUADRULE == DELTA_CS){
        KinematicVector delta = KinematicVector(x.VECTOR_TYPE);
        for (int i=0; i<x.size(); i++) {
            double trD = L[i].trace();
            delta = -alpha * job->dt * (L[i] - (trD / 3.0) * KinematicTensor::Identity(x.VECTOR_TYPE)).norm() * grad_e(i) * h * h;

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
    }

    return;
}


/*----------------------------------------------------------------------------*/
//not implemented yet
std::string DefaultPoints::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int DefaultPoints::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}

/*----------------------------------------------------------------------------*/
//
void DefaultPoints::generateLoads(Job* job, Body* body){
    //do nothing
    return;
}

void DefaultPoints::applyLoads(Job* job, Body* body){
    //do nothing
    return;
}
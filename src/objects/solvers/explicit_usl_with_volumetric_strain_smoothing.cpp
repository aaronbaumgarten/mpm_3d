//
// Created by aaron on 2/26/21.
// explicit_usl.cpp
//

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "solvers.hpp"
#include "objects/bodies/bodies.hpp"
#include "objects/grids/grids.hpp"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

void ExplicitUSLwithVolumetricStrainSmoothing::init(Job* job){
    //call parent
    ExplicitUSL::init(job);

    //loop over string vector for flags
    for (int s=0; s<str_props.size(); s++){
        if (str_props[s].compare("USE_MAST_CELLS")==0){
            //use Mast et al (2012) cell-wise antilocking
            use_mast_cells = true;
            use_zhang_cells = false;
            std::cout << "Using Mast et al. (2012) cell-wise anti-locking." << std::endl;
            continue;
        } else if (str_props[s].compare("USE_ZHANG_CELLS")==0){
            //use Mast et al (2012) cell-wise antilocking
            use_mast_cells = false;
            use_zhang_cells = true;
            std::cout << "Using Zhang et al. (2018) cell-wise anti-locking." << std::endl;
            continue;
        } else if (str_props[s].compare("ANTI_LOCKING_ONLY")==0){
            //do not smooth, only apply anti-locking
            anti_locking_only = true;
            std::cout << "Anti-locking only! No strain smoothing." << std::endl;
            continue;
        } else if (str_props[s].compare("USE_ZHANG_FLUX")==0){
            //use flux based smoothing from Zhang et al (2018)
            use_zhang_flux = true;
            use_zhang_cells = false;
            use_mast_cells = false;
            std::cout << "Using Zhang et al. (2018) flux-based anti-locking." << std::endl;
            if (job->JOB_TYPE != job->JOB_2D){
                std::cerr << "ERROR! Zhang et al. (2018) flux-based smoothing algorithm requires 2D problem! Exiting." << std::endl;
                exit(0);
            }
            continue;
        } else if (str_props[s].compare("WRITE_ERROR_METRICS") == 0){
            //write out error metrics
            write_error_metrics = true;
            //assume index 0 has output filename
            output_file = str_props[0];
            //assume end of int_props includes stride
            if (int_props.size() > 1) {
                stride = int_props[int_props.size() - 1];
            } else {
                stride = 1;
            }
            skip_counter = stride-1;
            //open file and write header
            std::ofstream file (output_file,std::ios::trunc);
            if (file.is_open()){
                //success!
                //write to file
                file << "Time, ||E||_2^2, max(E_i/V_i), KE, y_cm, ||v||_L2\n";

                file.close();
            } else {
                std::cerr << "ERROR! Cannot open " << output_file << "!" << std::endl;;
            }
        }
    }

    if (use_zhang_flux){
        //initialize grid properties
        Nx = Eigen::VectorXi(2);
        if (int_props.size() == 2){
            Nx(0) = int_props[0];
            Nx(1) = int_props[1];
        } else if (int_props.size() >= 4){
            Nx(0) = int_props[2];
            Nx(1) = int_props[3];
        } else {
            std::cerr << "ERROR! USE_ZHANG_FLUX given with insufficient int-props! Exiting: " << int_props.size() << std::endl;
            exit(0);
        }

        Lx = KinematicVector(job->JOB_TYPE);
        hx = Lx;
        if (fp64_props.size() == 2){
            Lx(0) = fp64_props[0];
            Lx(1) = fp64_props[1];
        } else {
            std::cerr << "ERROR! USE_ZHANG_FLUX given with insufficient properties! Exiting: " << fp64_props.size() << std::endl;
            exit(0);
        }
        hx(0) = Lx(0)/Nx(0);
        hx(1) = Lx(1)/hx(1);

        //initialize cell storage
        int cell_count = Nx(0)*Nx(1);
        vc = Eigen::VectorXd(cell_count);
        v0c = Eigen::VectorXd(cell_count);
        dc = Eigen::VectorXd(cell_count);
        ddc = Eigen::VectorXd(cell_count);
        Sc = Eigen::VectorXd(cell_count);
    }
}

void ExplicitUSLwithVolumetricStrainSmoothing::updateDensity(Job* job){
    //loop over bodies and apply strain smoothing
    Eigen::VectorXd vtmp;
    Eigen::VectorXd v_c = Eigen::VectorXd::Zero(job->grid->element_count);
    Eigen::VectorXd v0_c = Eigen::VectorXd::Zero(job->grid->element_count);
    Eigen::VectorXd d_c = Eigen::VectorXd(job->grid->element_count);
    std::vector<int> p_to_c;
    std::vector<std::vector<int>> c_to_p;
    MPMScalarSparseMatrix S;
    std::vector<int> nvec;
    std::vector<double> valvec;
    Eigen::VectorXd vn, dn, pvec, dp;
    int c = -1;
    KinematicVector tmpX;
    for (int b=0; b<job->bodies.size(); b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //check if anti-locking only
        if (anti_locking_only){
            //save material volumes
            vtmp = job->bodies[b]->points->v;
        }

        if (use_zhang_cells){
            //apply normal volume update
            for (int i=0;i<job->bodies[b]->points->v.rows();i++) {
                job->bodies[b]->points->v(i) *= std::exp(job->dt * job->bodies[b]->points->L(i).trace());
            }

            //intialize vector
            p_to_c = std::vector<int>(job->bodies[b]->points->v.rows());

            for (int p=0; p<job->bodies[b]->points->v.rows(); p++){
                //find the element id
                tmpX = job->bodies[b]->points->x[p];
                c = job->grid->whichElement(job,tmpX);

                //add id to map
                p_to_c[p] = c;

                if (c >= 0 && c < job->grid->element_count) {
                    //add point volumes to cell-wise sum
                    v_c(c) += job->bodies[b]->points->v(p);
                    v0_c(c) += job->bodies[b]->points->v0(p);
                }
            }

            //determine cell-wise strain
            for (int i=0; i<job->grid->element_count; i++){
                //may be nan, but only if no points map to this cell
                d_c(i) = v_c(i)/v0_c(i);
            }

            //apply smoothing
            for (int p=0; p<job->bodies[b]->points->v.rows(); p++){
                if (p_to_c[p] >= 0 && p_to_c[p] < job->grid->element_count) {
                    job->bodies[b]->points->v(p) = job->bodies[b]->points->v0(p) * d_c(p_to_c[p]);
                }
            }

        } else if (use_mast_cells){
            //apply update first
            for (int i=0;i<job->bodies[b]->points->v.rows();i++) {
                job->bodies[b]->points->v(i) *= std::exp(job->dt * job->bodies[b]->points->L(i).trace());
            }

            //initialize vector
            c_to_p = std::vector<std::vector<int>>(job->grid->element_count);

            //initialize mapping matrix
            S = MPMScalarSparseMatrix(job->bodies[b]->nodes->x.size(), job->bodies[b]->points->x.size());;

            //determine cell-wise mapping
            for (int p=0; p<job->bodies[b]->points->v.rows(); p++){
                //find the element id
                tmpX = job->bodies[b]->points->x[p];
                c = job->grid->whichElement(job,tmpX);

                //add id to map
                if (c >= 0 && c < job->grid->element_count) {
                    c_to_p[c].push_back(p);
                }
            }

            //weighted strain at points
            pvec = job->bodies[b]->points->v;
            for (int p=0; p<job->bodies[b]->points->v.rows(); p++){
                pvec(p) *= job->bodies[b]->points->v(p) / job->bodies[b]->points->v0(p);
            }

            //create local mapping (slooooooow, but maybe faster than using body mapping matrix?)
            for (int i=0; i<job->grid->element_count; i++){
                S.clear();
                for (int q = 0; q<c_to_p[i].size(); q++){
                    nvec.clear();
                    valvec.clear();
                    tmpX = job->bodies[b]->points->x[c_to_p[i][q]];
                    job->grid->evaluateBasisFnValue(job,tmpX,nvec,valvec);
                    for (int k=0; k<nvec.size(); k++){
                        S.push_back(nvec[k],c_to_p[i][q],valvec[k]);
                    }
                }

                //weighted node sum of strains
                vn = S*job->bodies[b]->points->v;
                dn = S*pvec;
                for (int j=0; j<job->bodies[b]->nodes->x.size(); j++){
                    dn(j) /= vn(j);
                }

                //map back to points
                dp = S.operate(dn,MPMSparseMatrixBase::TRANSPOSED);
                for (int q = 0; q<c_to_p[i].size(); q++){
                    job->bodies[b]->points->v(c_to_p[i][q]) = job->bodies[b]->points->v0(c_to_p[i][q]) * dp(c_to_p[i][q]);
                }
            }
        } else if (use_zhang_flux) {
            //apply normal volume update
            for (int i = 0; i < job->bodies[b]->points->v.rows(); i++) {
                job->bodies[b]->points->v(i) *= std::exp(job->dt * job->bodies[b]->points->L(i).trace());
            }

            //intialize vector
            p_to_c = std::vector<int>(job->bodies[b]->points->v.rows());

            for (int p = 0; p < job->bodies[b]->points->v.rows(); p++) {
                //find the element id
                tmpX = job->bodies[b]->points->x[p];
                c = CartesianLinear::cartesianWhichElement(job, tmpX, Lx, hx, Nx); //job->grid->whichElement(job, tmpX);

                //add id to map
                p_to_c[p] = c;

                if (c >= 0 && c < Nx(0)*Nx(1)) {
                    //add point volumes to cell-wise sum
                    vc(c) += job->bodies[b]->points->v(p);
                    v0c(c) += job->bodies[b]->points->v0(p);
                }
            }

            //determine cell-wise strain
            for (int i = 0; i <  Nx(0)*Nx(1); i++) {
                ddc(i) = 0;
                if (v0c(i) > 0) {
                    dc(i) = vc(i) / v0c(i);
                } else {
                    dc(i) = 0;
                }
            }

            //calculate fluxes in each direction, x-first, y-second
            int c_id = 0;
            bool plus_osc, minus_osc;
            int i_plus, i_minus, i_plus_plus;
            double g;
            //loop over each row of grid
            for (int j = 0; j < Nx(1); j++) {

                //determine Sc in each cell
                c_id = j*Nx(0);
                Sc(c_id) = (dc(c_id + 1) - dc(c_id))/hx(0);
                c_id = (j+1)*Nx(0) - 1;
                Sc(c_id) = (dc(c_id) - dc(c_id - 1))/hx(0);
                for (int i = 1; i < Nx(0)-1; i++) {
                    c_id = j*Nx(0) + i;
                    Sc(c_id) = (dc(c_id + 1) - dc(c_id - 1))/(2.0*hx(0));
                }

                //loop over edges
                plus_osc = true;
                for (int f=0; f<Nx(0)-1; f++){
                    minus_osc = plus_osc;
                    i_plus_plus = j*Nx(0) + f + 2;
                    i_plus = j*Nx(0) + f + 1;
                    i_minus = j*Nx(0) + f;

                    //check if plus side cell is oscillating
                    if (i_plus_plus == (j+1)*Nx(0) || (dc(i_plus) - dc(i_minus))*(dc(i_plus_plus) - dc(i_plus)) < 0){
                        plus_osc = true;
                    } else {
                        plus_osc = false;
                    }

                    //if oscillating, determine gamma
                    if (plus_osc && minus_osc){
                        g = (dc(i_minus) - dc(i_plus) + hx(0)*(Sc(i_minus) + Sc(i_plus))/2.0)/(dc(i_minus) - dc(i_plus));
                    } else {
                        g = 0;
                    }

                    //make sure g is valid
                    if (!std::isfinite(g)){
                        g = 0;
                    } else if (g < 0){
                        g = 0;
                    } else if (g > 1){
                        g = 1;
                    }

                    //adjust g by relevant scale:
                    g *= 0.25;

                    //calculate strain increment
                    ddc(i_minus) -= g * v0c(i_plus) * (dc(i_minus) - dc(i_plus)) / (v0c(i_minus) + v0c(i_plus));
                    ddc(i_plus) += g * v0c(i_minus) * (dc(i_minus) - dc(i_plus)) / (v0c(i_minus) + v0c(i_plus));
                }
            }

            //loop over each column of grid
            for (int i = 0; i < Nx(0); i++) {

                //determine Sc in each cell
                Sc(i) = (dc(Nx(0)+i) - dc(i))/hx(1);
                Sc((Nx(1)-1)*Nx(0) + i) = (dc((Nx(1)-1)*Nx(0) + i)) - dc((Nx(1)-2)*Nx(0) + i)/hx(1);
                for (int j = 1; j < Nx(1)-1; j++) {
                    c_id = j*Nx(0) + i;
                    Sc(c_id) = (dc(c_id + Nx(0)) - dc(c_id - Nx(0)))/(2.0*hx(1));
                }

                //loop over edges
                plus_osc = true;
                for (int f=0; f<Nx(1)-1; f++){
                    minus_osc = plus_osc;
                    i_plus_plus = (f+2)*Nx(0) + i;
                    i_plus = (f+1)*Nx(0) + i;
                    i_minus = f*Nx(0) + i;

                    //check if plus side cell is oscillating
                    if (i_plus_plus >= Nx(1)*Nx(0) || (dc(i_plus) - dc(i_minus))*(dc(i_plus_plus) - dc(i_plus)) < 0){
                        plus_osc = true;
                    } else {
                        plus_osc = false;
                    }

                    //if oscillating, determine gamma
                    if (plus_osc && minus_osc){
                        g = (dc(i_minus) - dc(i_plus) + hx(0)*(Sc(i_minus) + Sc(i_plus))/2.0)/(dc(i_minus) - dc(i_plus));
                    } else {
                        g = 0;
                    }

                    //make sure g is valid
                    if (!std::isfinite(g)){
                        g = 0;
                    } else if (g < 0){
                        g = 0;
                    } else if (g > 1){
                        g = 1;
                    }

                    //adjust g by relevant scale:
                    g *= 0.25;

                    //calculate strain increment
                    ddc(i_minus) -= g * v0c(i_plus) * (dc(i_minus) - dc(i_plus)) / (v0c(i_minus) + v0c(i_plus));
                    ddc(i_plus) += g * v0c(i_minus) * (dc(i_minus) - dc(i_plus)) / (v0c(i_minus) + v0c(i_plus));
                }
            }

            //update cell-wise averages
            for (int i=0; i<Nx(0)*Nx(1); i++){
                dc(i) += ddc(i);
            }

            //apply smoothing
            for (int p = 0; p < job->bodies[b]->points->v.rows(); p++) {
                if (p_to_c[p] >= 0 && p_to_c[p] < Nx(0)*Nx(1)) {
                    job->bodies[b]->points->v(p) = job->bodies[b]->points->v0(p) * dc(p_to_c[p]);
                }
            }
        }

        //update stress now
        job->bodies[b]->material->calculateStress(job, job->bodies[b].get(), Material::UPDATE);

        //if anti-locking only, reset volumes
        if (anti_locking_only){
            for (int p=0; p<job->bodies[b]->points->v.rows(); p++){
                job->bodies[b]->points->v(p) = vtmp(p) * std::exp(job->dt * job->bodies[b]->points->L(p).trace());
            }
        }
    }

    //update integrators as necessary
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        //this is new, but maybe useful
        job->bodies[b]->points->updateIntegrators(job,job->bodies[b].get());
    }
    return;
}

void ExplicitUSLwithVolumetricStrainSmoothing::updateStress(Job* job){
    //do not update stress here

    //but maybe write out error info if flag is set
    skip_counter += 1;
    if (write_error_metrics && skip_counter == stride){
        skip_counter = 0;
        //calculate and write error measures to output file for first body in simulation
        //BE CAREFUL!!!

        //initialize arrays
        Eigen::VectorXd V_i = Eigen::VectorXd(job->grid->node_count);
        Eigen::VectorXd v_i = Eigen::VectorXd(job->grid->node_count);
        Eigen::VectorXd e = Eigen::VectorXd(job->grid->node_count);
        Eigen::VectorXd H = Eigen::VectorXd(job->grid->node_count);

        //initialize figures of merit
        double H_norm = 0;
        double e_norm = 0;
        double v_L2 = 0;

        //get exact node volumes form grid
        for (int i=0; i<job->grid->node_count;i++){
            V_i(i) = job->grid->nodeVolume(job,i);
        }

        //get integrated node volume from points
        v_i = job->bodies[0]->S * job->bodies[0]->points->v;

        //calculate arrays
        double tmpNum;
        for (int i = 0; i < V_i.rows(); i++) {
            tmpNum = (v_i(i) - V_i(i));
            H(i) = std::max(0.0, tmpNum);
            e(i) = H(i) / V_i(i);
        }

        //calculate error measures
        KinematicVector tmpAcc = KinematicVector(job->JOB_TYPE);
        KinematicVector tmpVec = KinematicVector(job->JOB_TYPE);
        for (int i=0; i<V_i.rows(); i++){
            if (job->bodies[0]->nodes->m(i) > 0) {
                //||H||_2^2
                H_norm += H(i) * H(i);

                //||e||_\infty
                if (e(i) > e_norm) {
                    e_norm = e(i);
                }

                //||a^* - a^Q||_L2
                tmpVec = job->bodies[0]->nodes->x_t[i];
                v_L2 += tmpVec.dot(tmpVec) * V_i(i);
            }
        }
        //sqrt of ||v_err||_L2^2
        v_L2 = std::sqrt(v_L2);

        double ke = 0;
        double y_cm = 0;
        double m_cm = 0;
        for (int p=0; p<job->bodies[0]->points->x.size(); p++){
            ke += job->bodies[0]->points->m(p) * job->bodies[0]->points->x_t[p].dot(job->bodies[0]->points->x_t[p]);
            y_cm += job->bodies[0]->points->m(p) * job->bodies[0]->points->x(p,1);
            m_cm += job->bodies[0]->points->m(p);
        }
        ke *= 0.5;
        y_cm /= m_cm;

        //open and write to file
        std::ofstream file (output_file,std::ios::app);
        if (file.is_open()){
            //success!
            //write to file
            file << job->t << ", ";
            file << H_norm << ", ";
            file << e_norm << ", ";
            file << ke << ", ";
            file << y_cm << ", ";
            file << v_L2 << "\n";

            file.close();
        } else {
            std::cerr << "ERROR! Cannot open " << output_file << "!" << std::endl;;
        }
    }

    return;
}
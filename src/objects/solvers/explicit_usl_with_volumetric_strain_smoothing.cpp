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
        }
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

            //apply normal volume update
            for (int i=0;i<job->bodies[b]->points->v.rows();i++) {
                job->bodies[b]->points->v(i) *= std::exp(job->dt * job->bodies[b]->points->L(i).trace());
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
    return;
}
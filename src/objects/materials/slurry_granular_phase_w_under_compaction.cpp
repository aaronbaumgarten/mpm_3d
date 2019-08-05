//
// Created by aaron on 8/2/19.
// slurry_granular_phase_w_under_compaction.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <dlfcn.h>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "materials.hpp"

/*----------------------------------------------------------------------------*/
//
void SlurryGranularPhase_wUnderCompaction::init(Job* job, Body* body){
    //call main initialization
    SlurryGranularPhase::init(job,body);

    //add phi_c to definition
    if (fp64_props.size() < 14){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 1 more property defined ({phi_c}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        phi_c = fp64_props[13];

        printf("Additional material properties (phi_c = %g).\n",
               phi_c);
    }

    return;
}

/*----------------------------------------------------------------------------*/
//
void SlurryGranularPhase_wUnderCompaction::calculateStress(Job* job, Body* body, int SPEC){
    MaterialTensor T, T_tr, L, D, W;
    MaterialTensor tmpMat;

    double trD, tau_bar, tau_bar_tr, p, p_tr;
    double beta, mu, phi_eq, xi_dot_1, xi_dot_2;
    double tau_bar_k, p_k, tau_bar_kplus, p_kplus;
    double I_v_tr, I_tr, I_m_tr, gammap_dot_tr;
    double tmpVal;

    double dgamma_dtau, dI_dtau, dIv_dtau, dIm_dtau, dmu_dtau, dbeta_dtau, dxi_dtau;
    double dI_dp, dIv_dp, dIm_dp, dmu_dp, dbeta_dp, dxi_dp;
    double dr1_dtau, dr2_dtau, dr1_dp, dr2_dp;
    double dr_dp, r_p, r_p_tr, b_p;

    Eigen::Vector2d r, r_tr, b, upd;
    Eigen::Matrix2d dr;

    double p_tmp, tau_bar_tmp, lambda_tmp, gammap_dot_tmp;

    bool is_solved;
    int k;

    //calculate packing fraction
    for (int i=0;i<body->points->x.size();i++){
        phi(i) = 1/grains_rho * body->points->m(i) / body->points->v(i);
    }

    //interpolate eta from fluid points
    if (eta_0 > 0) {
        if (fluid_body_id >= 0) {
            Eigen::VectorXd pvec = job->bodies[fluid_body_id]->points->m * eta_0;
            Eigen::VectorXd nvec(job->bodies[fluid_body_id]->nodes->m.rows());
            nvec = job->bodies[fluid_body_id]->S * pvec;
            for (int i = 0; i < nvec.rows(); i++) {
                if (job->bodies[fluid_body_id]->nodes->m(i) > 0) {
                    nvec(i) = nvec(i) / job->bodies[fluid_body_id]->nodes->m(i);
                } else {
                    nvec(i) = 0;
                }
            }
            eta = body->S.operate(nvec, MPMSparseMatrixBase::TRANSPOSED);
        } else {
            for (int i=0;i<eta.rows();i++){
                eta(i) = eta_0;
            }
        }
    } else {
        eta.setZero();
    }

    //calculate stress state
    for (int i=0;i<body->points->x.size();i++) {
        if (body->points->active[i] == 0) {
            continue;
        }

        L = body->points->L(i);
        T = body->points->T(i);

        D = 0.5 * (L + L.transpose());
        W = 0.5 * (L - L.transpose());

        trD = D.trace();

        //trial stress
        T_tr = T + job->dt * ((2 * G * D) + (lambda * trD * MaterialTensor::Identity()) + (W * T) - (T * W));
        tau_bar_tr = (T_tr - (T_tr.trace() / 3.0) * MaterialTensor::Identity()).norm() / std::sqrt(2.0);
        p_tr = -T_tr.trace() / 3.0;

        //state determined?
        is_solved = false;

        //check admissibility ****************************************************************************************//
        if (!is_solved){
            beta = getBeta(phi(i), phi_m); //zero plastic flow limit

            if ((p_tr >= 0) && (tau_bar_tr <= ((mu_1 + beta)*p_tr))){
                gammap_dot_tr = 0;
                p = p_tr;
                tau_bar = tau_bar_tr;
                is_solved = true; //exit condition
            }
        }

        //check f2 yield surface *************************************************************************************//
        if (!is_solved){
            beta = getBeta(phi(i), 0.0);

            if ((p_tr + (K*tau_bar_tr/G)*beta) <= 0){
                gammap_dot_tr = tau_bar_tr/(G*job->dt);
                p = 0;
                tau_bar = 0;
                is_solved = true; //exit condition
            }
        }

        //check zero strength limit on f1 only ***********************************************************************//
        if (!is_solved){
            beta = getBeta(phi(i), phi_m); //high pressure limit

            //setup newton method to find pressure
            r_p = std::max(std::abs(p_tr),std::abs(p_tr + (K*tau_bar_tr/G)*beta));
            b_p = r_p;
            tau_bar_k = 0;
            p_k = 0;
            k = 0;

            //newton method
            while (std::abs(r_p) > ABS_TOL and std::abs(r_p)/std::abs(b_p) > REL_TOL){
                k += 1;

                //calculate gammap_dot
                gammap_dot_tr = tau_bar_tr / (G*job->dt);

                //calculate state for step
                calcState(gammap_dot_tr,p_k,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

                //calculate residual
                r_p = p_k - p_tr - K*job->dt*beta*gammap_dot_tr;

                //approx gradients
                p_tmp = getStep(p_k,NAN,0.0);
                calcState(gammap_dot_tr,p_tmp,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);
                dr_dp = (p_tmp - p_tr - K*job->dt*beta*gammap_dot_tr - r_p)/(p_tmp - p_k);

                //limit step length
                lambda_tmp = 1.0;
                r_p_tr = r_p;
                do{

                    //update guess
                    p_kplus = p_k - lambda_tmp * r_p / dr_dp;

                    if (p_kplus < 0){
                        p_kplus = 0;
                    }

                    //calculate new state
                    calcState(gammap_dot_tr,p_kplus,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

                    //test residual
                    r_p_tr = p_kplus - p_tr - K*job->dt*beta*gammap_dot_tr;

                    //set new lambda
                    lambda_tmp *= 0.5;

                } while(std::abs(r_p_tr) >= std::abs(r_p) && lambda_tmp > REL_TOL);

                if (k > 100 && std::abs(p_k - p_kplus) < ABS_TOL){
                    p_k = p_kplus;
                    r_p = r_p_tr;
                    break;
                }

                p_k = p_kplus;
                r_p = r_p_tr;
            }

            //check that zero strength assumption is true and admissible
            //new limit imposed on zero compressive strength material
            gammap_dot_tr = tau_bar_tr / (G*job->dt);
            if (((mu + beta) <= 0) and ((phi(i) >= phi_c) or (p_k <= (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*(gammap_dot_tr*gammap_dot_tr*grains_d*grains_d*grains_rho + 2.0*eta(i)*gammap_dot_tr) ))){
                p = p_k;
                tau_bar = tau_bar_k;
                is_solved = true;
            }
        }

        //check f1 yield surface *************************************************************************************//
        if (!is_solved){
            //setup initial conditions for newton solve
            beta = getBeta(phi(i),phi_m);
            tau_bar_k = 0;
            p_k = 0;
            //intial residual is arbitrary
            r(0) = tau_bar_tr;
            r(1) = -p_tr - (K*tau_bar_tr/G)*beta;
            b = r;

            k = 0;
            //newton method
            while (r.norm() > b.norm() * REL_TOL && r.norm() > ABS_TOL) {
                k+=1;

                //calculate equiv shear rate
                gammap_dot_tr = (tau_bar_tr - tau_bar_k)/(G*job->dt);

                //calculate state
                calcState(gammap_dot_tr, p_k, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                //calculate residual
                r(0) = tau_bar_k - (mu + beta)*p_k;
                r(1) = p_k - p_tr - K*job->dt*beta*gammap_dot_tr;

                //approx gradients
                tau_bar_tmp = getStep(tau_bar_k, tau_bar_tr, 0.0);
                gammap_dot_tmp = (tau_bar_tr - tau_bar_tmp)/(G*job->dt);
                calcState(gammap_dot_tmp,p_k,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);
                dr1_dtau = (tau_bar_tmp - (mu + beta)*p_k - r(0))/(tau_bar_tmp - tau_bar_k);
                dr2_dtau = (p_k - p_tr - K*job->dt*beta*gammap_dot_tmp - r(1))/(tau_bar_tmp - tau_bar_k);

                p_tmp = getStep(p_k, NAN, 0.0);
                calcState(gammap_dot_tr, p_tmp, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);
                dr1_dp = (tau_bar_k - (mu+beta)*p_tmp - r(0))/(p_tmp - p_k);
                dr2_dp = (p_tmp - p_tr - K*job->dt*beta*gammap_dot_tr - r(1))/(p_tmp - p_k);

                //form jacobian
                dr(0,0) = dr1_dtau;
                dr(0,1) = dr1_dp;
                dr(1,0) = dr2_dtau;
                dr(1,1) = dr2_dp;

                //calculate update and limit newton step
                upd = dr.inverse()*r;
                lambda_tmp = 1.0;
                do {
                    //update variables
                    tau_bar_kplus = tau_bar_k - lambda_tmp*upd(0);
                    p_kplus = p_k - lambda_tmp*upd(1);

                    if (p_kplus < 0){
                        p_kplus = 0;
                    }

                    if (tau_bar_kplus > tau_bar_tr){
                        tau_bar_kplus = tau_bar_tr;
                    }

                    //calculate strain rate
                    gammap_dot_tr = (tau_bar_tr - tau_bar_kplus)/(G*job->dt);

                    //calculate state
                    calcState(gammap_dot_tr,p_kplus,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

                    r_tr(0) = tau_bar_kplus - (mu+beta)*p_kplus;
                    r_tr(1) = p_kplus - p_tr - K*job->dt*beta*gammap_dot_tr;

                    lambda_tmp *= 0.5;
                } while (r_tr.norm() >= r.norm() && lambda_tmp > REL_TOL);

                //std::cout << p_k << ", " << tau_bar_k << ", " << r.norm() << ", " << r_tr.norm() << std::endl;

                if (k > 100 && std::abs(p_k - p_kplus) < ABS_TOL && std::abs(tau_bar_k - tau_bar_kplus) < ABS_TOL){
                    p_k = p_kplus;
                    tau_bar_k = tau_bar_kplus;
                    r = r_tr;
                    break;
                }

                p_k = p_kplus;
                tau_bar_k = tau_bar_kplus;
                r = r_tr;
            }
            //std::cout << gammap_dot_tr << ", " << p_kplus << ", " << eta(i) << ", " << phi(i) << std::endl;
            //std::cout << mu << ", " << beta << ", " << p_k << ", " << tau_bar_k << std::endl;

            //check that solution meets criteria for f1 yield ONLY
            gammap_dot_tr = (tau_bar_tr - tau_bar_k) / (G*job->dt);
            if ((phi(i) >= phi_c) or (p_k <= (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*(gammap_dot_tr*gammap_dot_tr*grains_d*grains_d*grains_rho + 2.0*eta(i)*gammap_dot_tr) )){
                p = p_k;
                tau_bar = tau_bar_k;
                is_solved = true;
                //std::cout << "success ^" << std::endl;
            }
        }

        //check zero strength f1,f3 yield condition ******************************************************************//
        if (!is_solved){
            //setup newton method to find pressure
            r_p = std::max(std::abs(p_tr),std::abs(p_tr + (K*tau_bar_tr/G)*beta));
            b_p = r_p;
            tau_bar_k = 0;
            p_k = 0;

            k = 0;
            //newton method
            while (std::abs(r_p) > ABS_TOL and std::abs(r_p)/std::abs(b_p) > REL_TOL){
                k += 1;

                //calculate gammap_dot
                gammap_dot_tr = tau_bar_tr / (G*job->dt);

                //calculate state for step
                calcState(gammap_dot_tr,p_k,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

                //calculate xi_dot_2
                xi_dot_2 = (p_k - p_tr)/(K*job->dt) - beta*gammap_dot_tr;

                //calculate residual
                r_p = p_k - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((gammap_dot_tr - K_4*xi_dot_2)*(gammap_dot_tr - K_4*xi_dot_2)*grains_d*grains_d*grains_rho + 2.0*eta(i)*(gammap_dot_tr - K_4*xi_dot_2));


                //approx gradients
                p_tmp = getStep(p_k, NAN, 0.0);
                calcState(gammap_dot_tr,p_tmp,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);
                xi_dot_2 = (p_tmp - p_tr)/(K*job->dt) - beta*gammap_dot_tr;
                dr_dp = (p_tmp - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((gammap_dot_tr - K_4*xi_dot_2)*(gammap_dot_tr - K_4*xi_dot_2)*grains_d*grains_d*grains_rho + 2.0*eta(i)*(gammap_dot_tr - K_4*xi_dot_2)) - r_p)/(p_tmp - p_k);

                //limit step length
                lambda_tmp = 1.0;
                r_p_tr = r_p;
                do{
                    //update guess
                    p_kplus = p_k - lambda_tmp * r_p / dr_dp;

                    if (p_kplus < 0){
                        p_kplus = 0;
                    }

                    //calculate new state
                    calcState(gammap_dot_tr,p_kplus,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

                    //test residual
                    r_p_tr = p_kplus - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((gammap_dot_tr - K_4*xi_dot_2)*(gammap_dot_tr - K_4*xi_dot_2)*grains_d*grains_d*grains_rho + 2.0*eta(i)*(gammap_dot_tr - K_4*xi_dot_2));

                    //set new lambda
                    lambda_tmp *= 0.5;

                } while(std::abs(r_p_tr) >= std::abs(r_p) && lambda_tmp > REL_TOL);

                if (k > 100 && std::abs(p_k - p_kplus) < ABS_TOL) {
                    p_k = p_kplus;
                    r_p = r_p_tr;
                    break;
                }

                p_k = p_kplus;
                r_p = r_p_tr;
            }

            //check that zero strength assumption is true
            gammap_dot_tr = tau_bar_tr / (G*job->dt);
            if ((mu + beta) <= 0){
                p = p_k;
                tau_bar = tau_bar_k;
                is_solved = true;
            }
        }

        //check f1,f3 yield surface **********************************************************************************//
        if (!is_solved){
            //setup initial conditions for newton solve
            beta = getBeta(phi(i),0.0);
            tau_bar_k = 0;
            p_k = 0;
            //intial residual is arbitrary
            r(0) = tau_bar_tr;
            r(1) = p_tr;
            b = r;

            k = 0;

            //newton method
            while (r.norm() > b.norm() * REL_TOL && r.norm() > ABS_TOL) {
                k+=1;

                //calculate equiv shear rate
                gammap_dot_tr = (tau_bar_tr - tau_bar_k)/(G*job->dt);

                //calculate state
                calcState(gammap_dot_tr, p_k, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                //calculate xi_dot_2
                xi_dot_2 = (p_k - p_tr)/(K*job->dt) - beta*gammap_dot_tr;

                //calculate residual
                r(0) = tau_bar_k - (mu + beta)*p_k;
                r(1) = p_k - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((gammap_dot_tr - K_4*xi_dot_2)*(gammap_dot_tr - K_4*xi_dot_2)*grains_d*grains_d*grains_rho + 2.0*eta(i)*(gammap_dot_tr - K_4*xi_dot_2));

                //approx gradients
                tau_bar_tmp = getStep(tau_bar_k, tau_bar_tr, 0.0);
                gammap_dot_tmp = (tau_bar_tr - tau_bar_tmp)/(G*job->dt);
                calcState(gammap_dot_tmp,p_k,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);
                xi_dot_2 = (p_k - p_tr)/(K*job->dt) - beta*gammap_dot_tmp;
                dr1_dtau = (tau_bar_tmp - (mu + beta)*p_k - r(0))/(tau_bar_tmp - tau_bar_k);
                dr2_dtau = (p_k - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((gammap_dot_tmp - K_4*xi_dot_2)*(gammap_dot_tmp - K_4*xi_dot_2)*grains_d*grains_d*grains_rho + 2.0*eta(i)*(gammap_dot_tmp - K_4*xi_dot_2)) - r(1))/(tau_bar_tmp - tau_bar_k);

                p_tmp = getStep(p_k, NAN, 0.0);
                calcState(gammap_dot_tr, p_tmp, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);
                xi_dot_2 = (p_tmp - p_tr)/(K*job->dt) - beta*gammap_dot_tr;
                dr1_dp = (tau_bar_k - (mu+beta)*p_tmp - r(0))/(p_tmp - p_k);
                dr2_dp = (p_tmp - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((gammap_dot_tr - K_4*xi_dot_2)*(gammap_dot_tr - K_4*xi_dot_2)*grains_d*grains_d*grains_rho + 2.0*eta(i)*(gammap_dot_tr - K_4*xi_dot_2)) - r(1))/(p_tmp - p_k);

                //form jacobian
                dr(0,0) = dr1_dtau;
                dr(0,1) = dr1_dp;
                dr(1,0) = dr2_dtau;
                dr(1,1) = dr2_dp;

                //calculate update and limit newton step
                upd = dr.inverse()*r;
                lambda_tmp = 1.0;
                do {
                    //update variables
                    tau_bar_kplus = tau_bar_k - lambda_tmp*upd(0);
                    p_kplus = p_k - lambda_tmp*upd(1);

                    if (p_kplus < 0){
                        p_kplus = 0;
                    }
                    if (tau_bar_kplus > tau_bar_tr){
                        tau_bar_kplus = tau_bar_tr;
                    }

                    //calculate strain rate
                    gammap_dot_tr = (tau_bar_tr - tau_bar_kplus)/(G*job->dt);

                    //calculate state
                    calcState(gammap_dot_tr,p_kplus,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);
                    xi_dot_2 = (p_kplus - p_tr)/(K*job->dt) - beta*gammap_dot_tr;

                    //calculate residual
                    r_tr(0) = tau_bar_kplus - (mu + beta)*p_kplus;
                    r_tr(1) = p_kplus - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((gammap_dot_tr - K_4*xi_dot_2)*(gammap_dot_tr - K_4*xi_dot_2)*grains_d*grains_d*grains_rho + 2.0*eta(i)*(gammap_dot_tr - K_4*xi_dot_2));

                    lambda_tmp *= 0.5;
                } while (r_tr.norm() >= r.norm() && lambda_tmp > REL_TOL);

                if (k > 100 && std::abs(p_k - p_kplus) < ABS_TOL && std::abs(tau_bar_k - tau_bar_kplus) < ABS_TOL){
                    p_k = p_kplus;
                    tau_bar_k = tau_bar_kplus;
                    r = r_tr;
                    break;
                }

                p_k = p_kplus;
                tau_bar_k = tau_bar_kplus;
                r = r_tr;
            }

            //this is the last possible case, if we got here then this is the answer
            gammap_dot_tr = (tau_bar_tr - tau_bar_k) / (G*job->dt);
            p = p_k;
            tau_bar = tau_bar_k;
            is_solved = true;
        }

        if (!std::isfinite(p) || !std::isfinite(tau_bar)){
            std::cout << "u";
        }

        //update stress
        if (p > 0 && tau_bar > 0) {
            T = tau_bar / tau_bar_tr * (T_tr - T_tr.trace() / 3.0 * MaterialTensor::Identity());
            T = T - p * MaterialTensor::Identity();
        } else {
            T.setZero();
        }

        body->points->T(i) = T;

        if (SPEC == Material::UPDATE) {
            gammap_dot(i) = gammap_dot_tr;
            gammap(i) += job->dt * gammap_dot(i);
            I_v(i) = eta(i)*gammap_dot(i)/(-T.trace()/T.rows()); //undefined
            I(i) = gammap_dot(i)*grains_d*std::sqrt(grains_rho/(-T.trace()/T.rows()));
            I_m(i) = std::sqrt(I(i)*I(i) + 2*I_v(i));
        }

    }

    return;
}
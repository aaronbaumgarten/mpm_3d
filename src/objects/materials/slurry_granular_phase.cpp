//
// Created by aaron on 5/25/18.
// slurry_mu_im.cpp
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

bool SLURRY_MU_IM_DEBUG = false;

/*----------------------------------------------------------------------------*/
//
void SlurryGranularPhase::init(Job* job, Body* body){
    if (fp64_props.size() < 13){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 13 properties defined ({E, nu, mu_1, mu_2, a, b, K_3, K_4, phi_m, grains_rho, eta_0, grains_d, fluid_rho}, {fluid_body}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        E = fp64_props[0];
        nu = fp64_props[1];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        mu_1 = fp64_props[2];
        mu_2 = fp64_props[3];
        a = fp64_props[4];
        I_0 = fp64_props[5];
        K_3 = fp64_props[6];
        K_4 = fp64_props[7];
        //K_5 = fp64_props[8];
        phi_m = fp64_props[8];
        grains_rho = fp64_props[9];
        eta_0 = fp64_props[10];
        grains_d = fp64_props[11];
        fluid_rho = fp64_props[12];

        //set body id by name
        if (str_props.size() >= 1){
            for (int b = 0; b < job->bodies.size(); b++) {
                if (str_props[0].compare(job->bodies[b]->name) == 0){
                    fluid_body_id = b;
                    break;
                }
            }
        }

        if (fluid_body_id == -1){
            std::cout << std::endl << "WARNING: No fluid body defined!" << std::endl << std::endl;
        }

        if (K_3 > mu_1 / phi_m){
            std::cout << std::endl;
            std::cout << "WARNING: K_3 = " << K_3 << " > mu_1/phi_m = " << mu_1 / phi_m << "! For dilute suspensions, this may result in negative strength!" << std::endl;
            std::cout << std::endl;
        }

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, mu_1 = %g, mu_2 = %g, a = %g, I_0 = %g K_3 = %g, K_4 = %g, phi_m = %g, grains_rho = %g, eta_0 = %g, grains_d = %g, fluid_rho = %g).\n",
               E, nu, G, K, mu_1, mu_2, a, I_0, K_3, K_4, phi_m, grains_rho, eta_0, grains_d, fluid_rho);
    }

    gammap_dot.resize(body->points->x.size());
    gammap_dot.setZero();
    gammap.resize(body->points->x.size());
    gammap.setZero();
    I_v.resize(body->points->x.size());
    I_v.setZero();
    I.resize(body->points->x.size());
    I.setZero();
    I_m.resize(body->points->x.size());
    I_m.setZero();
    phi.resize(body->points->x.size());
    phi.setZero();
    eta.resize(body->points->x.size());
    eta.setZero();

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//
double SlurryGranularPhase::getBeta(double phi, double phi_eq){
    /*if (phi > phi_m) {
        return (K_3 * (phi - phi_m) + K_4 * (phi - phi_eq));
    } else {
        return (K_4 * (phi - phi_eq));
    }*/
    return K_3 * (phi - phi_eq);
}

void SlurryGranularPhase::calcState(double &gdp, double &p, double &eta_in, double &phi_in, double &I_out, double &Iv_out, double &Im_out, double &mu_out, double &phi_eq, double &beta_out){
    I_out = gdp * grains_d * std::sqrt(grains_rho/p);
    Iv_out = gdp * eta_in / p;
    Im_out = std::sqrt(I_out*I_out + 2 * Iv_out);

    if (Im_out == 0){
        mu_out = mu_1;
    } else if (p <= 0){
        mu_out = mu_2;
    } else {
        mu_out = mu_1 + (mu_2 - mu_1)/(1 + I_0/Im_out) + 2.5 * phi_in * Iv_out/(a*Im_out);
    }

    if (p <= 0){
        phi_eq = 0;
    } else {
        phi_eq = phi_m / (1 + a * Im_out);
    }

    /*if (phi_in > phi_m) {
        beta_out = (K_3 * (phi_in - phi_m) + K_4 * (phi_in - phi_eq));
    } else {
        beta_out = (K_4 * (phi_in - phi_eq));
    }*/

    beta_out = K_3 * (phi_in - phi_eq);

    return;
}

double SlurryGranularPhase::getStep(double val, double ub, double lb){
    if (val == 0){
        return h;
    } if (val*(1+h) >= ub){
        if (val*(1-h) <= lb && SLURRY_MU_IM_DEBUG){
            std::cout << "getStep() OUT OF BOUNDS!" << std::endl;
        }
        return val*(1-h);
    } else {
        return val*(1+h);
    }
}


/*----------------------------------------------------------------------------*/
//
void SlurryGranularPhase::writeFrame(Job* job, Body* body, Serializer* serializer){
    serializer->writeScalarArray(gammap,"gammap");
    serializer->writeScalarArray(gammap_dot,"gammap_dot");
    serializer->writeScalarArray(I_v,"I_v");
    serializer->writeScalarArray(I,"I");
    serializer->writeScalarArray(I_m,"I_m");
    serializer->writeScalarArray(phi,"packing_fraction");
    serializer->writeScalarArray(eta,"eta_interp");
    return;
}

std::string SlurryGranularPhase::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int SlurryGranularPhase::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//
void SlurryGranularPhase::calculateStress(Job* job, Body* body, int SPEC){
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
    double p_max, p_min, tau_max, tau_min;
    double rp_max, rp_min, rt_max, rt_min;
    double p_zero_strength_f1;
    double p_zero_strength_f1f3;

    bool tau_too_large = false;
    bool tau_too_small = false;

    double A, B, C;

    bool is_solved;
    int k, j;

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
        tau_bar = tau_bar_tr;
        p = p_tr;

        //check admissibility ****************************************************************************************//
        if (!is_solved){
            beta = getBeta(phi(i), phi_m); //zero plastic flow limit

            /*if (p_tr >= 0 && (tau_bar_tr <= ((mu_1 + beta)*p_tr))){
                //this would've been admissible in old code
                std::cout << "[" << i << "]: would've passed!" << std::endl;
            }*/

            if ((p_tr >= 0) && (tau_bar_tr <= ((mu_1 + beta)*p_tr)) && (phi(i) >= phi_m)){
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

            //setup bisection method to find pressure
            r_p = std::max(std::abs(p_tr),std::abs(p_tr + (K*tau_bar_tr/G)*beta));
            b_p = r_p;
            tau_bar_k = 0;
            p_k = 0;
            k = 0;


            //use bisection
            beta = getBeta(phi(i), 0.0); //low pressure limit
            p_max = p_tr + (K*tau_bar_tr/G)*beta; tau_max = tau_bar_tr;
            p_min = 0.0;  tau_min = 0.0;

            //calculate gammap_dot for tau going to zero
            gammap_dot_tr = tau_bar_tr / (G*job->dt);

            //calculate state for p_max
            calcState(gammap_dot_tr,p_max,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

            //calculate residual for p_max
            rp_max = p_max - p_tr - K*job->dt*beta*gammap_dot_tr;

            //calculate state for p_min
            calcState(gammap_dot_tr,p_min,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

            //calculate residual for r_min
            rp_min = p_min - p_tr - K*job->dt*beta*gammap_dot_tr;

            if (rp_min * rp_max > 0){
                std::cout << "ERROR: f1 weak residuals in binary search have same sign! That's bad!" << std::endl;
            }

            //bisection method
            while (std::abs(r_p) > ABS_TOL and std::abs(r_p)/std::abs(b_p) > REL_TOL){
                k += 1;
                if (k > 50 && SLURRY_MU_IM_DEBUG){
                    std::cout << "f1 weak: " << p_k << ", " << r_p << std::endl;
                    break; //machine tolerance with bisection
                } else if (k > 50){
                    break;
                }

                //binary search
                p_k = 0.5*(p_max + p_min);

                //calculate state for step
                calcState(gammap_dot_tr,p_k,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

                //calculate residual
                r_p = p_k - p_tr - K*job->dt*beta*gammap_dot_tr;

                //check residual sign
                if (r_p * rp_min > 0){
                    //r_p replaces r_min
                    rp_min = r_p;
                    p_min = p_k;
                } else {
                    //r_p replaces r_max
                    rp_max = r_p;
                    p_max = p_k;
                }
            }

            //save value for next check
            p_zero_strength_f1 = p_k;

            //check that zero strength assumption is true and admissible
            gammap_dot_tr = tau_bar_tr / (G*job->dt);
            if (((mu + beta) <= 0) and ((phi(i) >= phi_m) or (p_k <= (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*(gammap_dot_tr*gammap_dot_tr*grains_d*grains_d*grains_rho + 2.0*eta(i)*gammap_dot_tr) ))){
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

            //use bisection
            beta = getBeta(phi(i), 0.0); //low pressure limit
            p_max = p_tr + (K*tau_bar_tr/G)*beta;
            p_min = 0.0;
            tau_max = tau_bar_tr;
            tau_min = 0.0;

            //tau = tau_max -> gdp = 0 -> p = p_tr
            gammap_dot_tr = 0.0;
            calcState(gammap_dot_tr, p_tr, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);
            rt_max = tau_max - (mu + beta)*p_tr;

            //tau = tau_min -> p = p_zero_strenght
            gammap_dot_tr = tau_bar_tr/(G*job->dt);
            calcState(gammap_dot_tr, p_zero_strength_f1, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);
            rt_min = tau_min - (mu + beta)*p_zero_strength_f1;

            //if on f1 only and (phi(i) < phi_m) then beta >= 0
            //as tau decreases, p increases. If both rt_min and rt_max are less than zero,
            //then we must also be on f3
            //if rt_min > 0 and rt_max > 0, then p_zero_strength_f1 -> mu+beta < 0,  must be on f3, too
            if (rt_min*rt_max > 0){ //rt_min < 0 && rt_max < 0){
                //do nothing
            } else {
                //solve for f1

                if (rt_min * rt_max > 0) {
                    std::cout << "ERROR: f1 residuals in binary search have same sign! That's bad!" << std::endl;
                    std::cout << rt_min << ", " << rt_max << ", " << phi(i) << " <? " << phi_m << std::endl;
                }

                tau_too_large = false;
                k = 0;
                //bisection method
                while (std::abs(r(0)) > std::abs(b(0)) * REL_TOL && std::abs(r(0)) > ABS_TOL) {
                    k += 1;
                    if (k > 50 && SLURRY_MU_IM_DEBUG) {
                        std::cout << "f1: " << p_k << ", " << tau_bar_k << ", " << r(0) << ", " << r(1) << std::endl;
                        std::cout << rt_min << ", " << rt_max << std::endl;
                        break; //after 50 bisections, answer is not changing
                    } else if (k > 50) {
                        break;
                    }

                    //set up initial residual
                    beta = getBeta(phi(i), phi_m); //high pressure limit
                    r_p = std::max(std::abs(p_tr), std::abs(p_tr + (K * tau_bar_tr / G) * beta)); //reference
                    b_p = r_p;

                    //set tau_bar_k
                    tau_bar_k = 0.5 * (tau_max + tau_min);

                    //find p_j
                    //use bisection
                    beta = getBeta(phi(i), 0.0); //low pressure limit
                    p_max = p_tr + (K * tau_bar_tr / G) * beta;
                    p_min = 0.0;

                    //calculate gammap_dot for tau going to zero
                    gammap_dot_tr = (tau_bar_tr - tau_bar_k) / (G * job->dt);

                    //calculate state for p_max
                    calcState(gammap_dot_tr, p_max, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                    //calculate residual for p_max
                    rp_max = p_max - p_tr - K * job->dt * beta * gammap_dot_tr;

                    //calculate state for p_min
                    calcState(gammap_dot_tr, p_min, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                    //calculate residual for r_min
                    rp_min = p_min - p_tr - K * job->dt * beta * gammap_dot_tr;

                    tau_too_large = false;
                    if (rp_min * rp_max > 0) {
                        //std::cout << "ERROR: Residuals in binary search have same sign! That's bad!" << std::endl;
                        //this likely means that tau_bar_k is too large (gdp too small, beta too small)
                        tau_max = tau_bar_k;
                        tau_too_large = true;
                        r_p = 0;
                        //std::cout << p_tr << " : " << p_max << " : " << p_min << std::endl;
                        //exit(0);
                    }

                    //bisection method
                    j = 0;
                    while (std::abs(r_p) > ABS_TOL and std::abs(r_p) / std::abs(b_p) > REL_TOL) {
                        j += 1;
                        if (j > 50 && SLURRY_MU_IM_DEBUG) {
                            std::cout << "f1 internal: " << p_k << ", " << r_p << ", " << b_p << std::endl;
                            break; //after 50 bisections, answer is not changing
                        } else if (j > 50) {
                            break;
                        }

                        //binary search
                        p_k = 0.5 * (p_max + p_min);

                        //calculate state for step
                        calcState(gammap_dot_tr, p_k, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                        //calculate residual
                        r_p = p_k - p_tr - K * job->dt * beta * gammap_dot_tr;

                        //check residual sign
                        if (r_p * rp_min > 0) {
                            //r_p replaces r_min
                            rp_min = r_p;
                            p_min = p_k;
                        } else {
                            //r_p replaces r_max
                            rp_max = r_p;
                            p_max = p_k;
                        }
                    }

                    if (!tau_too_large) {
                        //calculate equiv shear rate
                        gammap_dot_tr = (tau_bar_tr - tau_bar_k) / (G * job->dt);

                        //calculate state
                        calcState(gammap_dot_tr, p_k, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                        //calculate residual
                        r(0) = tau_bar_k - (mu + beta) * p_k;
                        r(1) = p_k - p_tr - K * job->dt * beta * gammap_dot_tr;

                        //check sign of residual
                        if (r(0) * rt_min > 0) {
                            //r(0) replaces rt_min
                            rt_min = r(0);
                            tau_min = tau_bar_k;
                        } else {
                            rt_max = r(0);
                            tau_max = tau_bar_k;
                        }
                    }
                }
                //should now have solution


                //std::cout << gammap_dot_tr << ", " << p_kplus << ", " << eta(i) << ", " << phi(i) << std::endl;
                //std::cout << mu << ", " << beta << ", " << p_k << ", " << tau_bar_k << std::endl;

                //check that solution meets criteria for f1 yield ONLY
                gammap_dot_tr = (tau_bar_tr - tau_bar_k) / (G * job->dt);
                if ((phi(i) >= phi_m) or (p_k <= (a * a * phi(i) * phi(i)) / ((phi_m - phi(i)) * (phi_m - phi(i))) *
                                                 (gammap_dot_tr * gammap_dot_tr * grains_d * grains_d * grains_rho +
                                                  2.0 * eta(i) * gammap_dot_tr))) {
                    p = p_k;
                    tau_bar = tau_bar_k;
                    is_solved = true;
                    //std::cout << "success ^" << std::endl;
                }
            }
        }

        //check zero strength f1,f3 yield condition ******************************************************************//
        if (!is_solved){
            //setup newton method to find pressure
            beta = getBeta(phi(i), phi_m); //high pressure limit
            r_p = std::max(std::abs(p_tr),std::abs(p_tr + (K*tau_bar_tr/G)*beta));
            b_p = r_p;
            tau_bar_k = 0;
            p_k = 0;

            //use bisection
            p_max = p_tr; tau_max = tau_bar_tr;
            p_min = 0.0;  tau_min = 0.0;

            //calculate gammap_dot for tau going to zero
            gammap_dot_tr = tau_bar_tr / (G*job->dt);

            //calculate state for p_max
            calcState(gammap_dot_tr,p_max,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

            //if mu + beta is positive at maximum admissible pressure,
            //then mu + beta will be positive for entire range of p
            if (mu_1 + beta > 0){
                //do nothing
                p_zero_strength_f1f3 = p_max;
            } else {

                //calculate xi_dot_2
                xi_dot_2 = (p_max - p_tr) / (K * job->dt) - beta * gammap_dot_tr;

                //calculate residual
                rp_max = p_max - (a * a * phi(i) * phi(i)) / ((phi_m - phi(i)) * (phi_m - phi(i))) *
                                 ((gammap_dot_tr - K_4 * xi_dot_2) * (gammap_dot_tr - K_4 * xi_dot_2) * grains_d *
                                  grains_d * grains_rho + 2.0 * eta(i) * (gammap_dot_tr - K_4 * xi_dot_2));

                //calculate state for p_min
                calcState(gammap_dot_tr, p_min, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                //calculate xi_dot_2
                xi_dot_2 = (p_min - p_tr) / (K * job->dt) - beta * gammap_dot_tr;

                //calculate residual
                rp_min = p_min - (a * a * phi(i) * phi(i)) / ((phi_m - phi(i)) * (phi_m - phi(i))) *
                                 ((gammap_dot_tr - K_4 * xi_dot_2) * (gammap_dot_tr - K_4 * xi_dot_2) * grains_d *
                                  grains_d * grains_rho + 2.0 * eta(i) * (gammap_dot_tr - K_4 * xi_dot_2));

                if (rp_min * rp_max > 0) {
                    std::cout << "ERROR: f1,f3 weak residuals in binary search have same sign! That's bad!"
                              << std::endl;
                    std::cout << p_min << ": " << rp_min << ", " << p_max << ": " << rp_max << std::endl;
                }

                k = 0;
                //bisection method
                while (std::abs(r_p) > ABS_TOL and std::abs(r_p) / std::abs(b_p) > REL_TOL) {
                    k += 1;
                    if (k > 50 && SLURRY_MU_IM_DEBUG) {
                        std::cout << "f1, f3: " << p_k << ", " << r_p << std::endl;
                        break; //after 50 iterations, bisection is not changing
                    } else if (k > 50) {
                        break;
                    }

                    //binary search
                    p_k = 0.5 * (p_max + p_min);

                    //calculate state for step
                    calcState(gammap_dot_tr, p_k, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                    //calculate xi_dot_2
                    xi_dot_2 = (p_k - p_tr) / (K * job->dt) - beta * gammap_dot_tr;

                    //calculate residual
                    r_p = p_k - (a * a * phi(i) * phi(i)) / ((phi_m - phi(i)) * (phi_m - phi(i))) *
                                ((gammap_dot_tr - K_4 * xi_dot_2) * (gammap_dot_tr - K_4 * xi_dot_2) * grains_d *
                                 grains_d * grains_rho + 2.0 * eta(i) * (gammap_dot_tr - K_4 * xi_dot_2));

                    //check residual sign
                    if (r_p * rp_min > 0) {
                        //r_p replaces r_min
                        rp_min = r_p;
                        p_min = p_k;
                    } else {
                        //r_p replaces r_max
                        rp_max = r_p;
                        p_max = p_k;
                    }
                }

                //save value for next check
                p_zero_strength_f1f3 = p_k;

                //check that zero strength assumption is true
                gammap_dot_tr = tau_bar_tr / (G * job->dt);
                if ((mu + beta) <= 0) {
                    p = p_k;
                    tau_bar = tau_bar_k;
                    is_solved = true;
                }
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

            //use bisection
            beta = getBeta(phi(i), 0.0); //low pressure limit
            p_max = p_tr; tau_max = tau_bar_tr;
            p_min = 0.0;  tau_min = 0.0;

            //tau = tau_max -> gdp = 0 -> p = ?
            gammap_dot_tr = 0.0;
            A = - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((K_4)*(K_4)*grains_d*grains_d*grains_rho);
            B = (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*(2.0*eta(i)*(K_4)) + K*job->dt;
            C = p_tr;
            if (K_4 < 1e-10){
                p_max = 0.0;
            } else {
                xi_dot_2 = (-B - std::sqrt(B*B - 4*A*C))/(2*A);
                p_max = p_tr - K*xi_dot_2*job->dt;
            }
            calcState(gammap_dot_tr, p_max, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);
            rt_max = tau_max - (mu + beta)*p_max;

            //tau = tau_min -> p = p_zero_strength
            gammap_dot_tr = tau_bar_tr/(G*job->dt);
            calcState(gammap_dot_tr, p_zero_strength_f1f3, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);
            rt_min = tau_min - (mu + beta)*p_zero_strength_f1f3;

            if (rt_min * rt_max > 0){
                std::cout << "ERROR: f1,f3 residuals in binary search have same sign! That's bad!" << std::endl;
                std::cout << rt_min << ", " << rt_max << std::endl;
                std::cout << (mu + beta) << ", " << p_zero_strength_f1f3 << std::endl;
                //tau_max -> gdp=0 -> p_max small -> rt_max big
                //tau_min -> max gdp -> p_max big -> rt_max small
                //if rt_min > 0, then (mu + beta) associated with p_zero_strength < 0, shouldn't happen
            }

            k = 0;
            tau_too_small = false;

            //bisection method
            while (r.norm() > b.norm() * REL_TOL && r.norm() > ABS_TOL) {
                k+=1;
                if (k>50 && SLURRY_MU_IM_DEBUG){
                    //std::cout << "f1,f3: " << p_k << ", " << tau_bar_k << ", " << r.norm() << ", " << dr.norm() << std::endl;
                    std::cout << "f1,f3: " << p_k << ", " << tau_bar_k << std::endl;
                    std::cout << "       " << r(0) << ", " << r(1) << std::endl;
                    std::cout << "       " << tau_too_small << std::endl;
                    //after 50 iterations bisection isn't changing
                    break;
                } else if (k>50){
                    break;
                }

                //set tau_k
                tau_bar_k = 0.5*(tau_max + tau_min);

                //calculate equiv shear rate
                gammap_dot_tr = (tau_bar_tr - tau_bar_k)/(G*job->dt);

                //use bisection
                p_max = p_tr;
                p_min = 0.0;

                //calculate state for p_max
                calcState(gammap_dot_tr,p_max,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

                //if beta is positive at max pressure, then beta will be positive for all admissible pressures
                tau_too_small = false;
                if (beta > 0){
                    //trial tau too small (gdp too large)
                    tau_min = tau_bar_k;
                    tau_too_small = true;
                }

                //calculate xi_dot_2
                xi_dot_2 = (p_max - p_tr)/(K*job->dt) - beta*gammap_dot_tr;

                //calculate residual
                rp_max = p_max - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((gammap_dot_tr - K_4*xi_dot_2)*(gammap_dot_tr - K_4*xi_dot_2)*grains_d*grains_d*grains_rho + 2.0*eta(i)*(gammap_dot_tr - K_4*xi_dot_2));

                //calculate state for p_min
                calcState(gammap_dot_tr,p_min,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

                //calculate xi_dot_2
                xi_dot_2 = (p_min - p_tr)/(K*job->dt) - beta*gammap_dot_tr;

                //calculate residual
                rp_min = p_min - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((gammap_dot_tr - K_4*xi_dot_2)*(gammap_dot_tr - K_4*xi_dot_2)*grains_d*grains_d*grains_rho + 2.0*eta(i)*(gammap_dot_tr - K_4*xi_dot_2));

                if (rp_min * rp_max > 0 && !tau_too_small){
                    std::cout << "ERROR: f1,f3 pressure residuals in binary search have same sign! That's bad!" << std::endl;
                }

                j = 0;
                if (!tau_too_small) {
                    r_p = rp_max;
                    b_p = r_p;
                } else {
                    p_k = p_max;
                    b_p = rp_max;
                    r_p = 0.0;
                }
                //bisection method
                while (std::abs(r_p) > ABS_TOL and std::abs(r_p)/std::abs(b_p) > REL_TOL){
                    j += 1;
                    if (j > 50 && SLURRY_MU_IM_DEBUG){
                        std::cout << "f1, f3: " << p_k << ", " << r_p << std::endl;
                        break; //bisection at machine precision
                    } else if (j > 50){
                        break;
                    }

                    //binary search
                    p_k = 0.5*(p_max + p_min);

                    //calculate state for step
                    calcState(gammap_dot_tr,p_k,eta(i),phi(i),I_tr,I_v_tr,I_m_tr,mu,phi_eq,beta);

                    //calculate xi_dot_2
                    xi_dot_2 = (p_k - p_tr)/(K*job->dt) - beta*gammap_dot_tr;

                    //calculate residual
                    r_p = p_k - (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*((gammap_dot_tr - K_4*xi_dot_2)*(gammap_dot_tr - K_4*xi_dot_2)*grains_d*grains_d*grains_rho + 2.0*eta(i)*(gammap_dot_tr - K_4*xi_dot_2));

                    //check residual sign
                    if (r_p * rp_min > 0){
                        //r_p replaces r_min
                        rp_min = r_p;
                        p_min = p_k;
                    } else {
                        //r_p replaces r_max
                        rp_max = r_p;
                        p_max = p_k;
                    }
                }

                if (!tau_too_small) {
                    //calculate state
                    calcState(gammap_dot_tr, p_k, eta(i), phi(i), I_tr, I_v_tr, I_m_tr, mu, phi_eq, beta);

                    //calculate xi_dot_2
                    xi_dot_2 = (p_k - p_tr) / (K * job->dt) - beta * gammap_dot_tr;

                    //calculate residual
                    r(0) = tau_bar_k - (mu + beta) * p_k;
                    r(1) = p_k - (a * a * phi(i) * phi(i)) / ((phi_m - phi(i)) * (phi_m - phi(i))) *
                                 ((gammap_dot_tr - K_4 * xi_dot_2) * (gammap_dot_tr - K_4 * xi_dot_2) * grains_d *
                                  grains_d * grains_rho + 2.0 * eta(i) * (gammap_dot_tr - K_4 * xi_dot_2));


                    //check sign of residual
                    if (r(0) * rt_max > 0) {
                        //r(1) replaces rt_max
                        rt_max = r(0);
                        tau_max = tau_bar_k;
                    } else {
                        rt_min = r(0);
                        tau_min = tau_bar_k;
                    }
                }
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

void SlurryGranularPhase::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    body->points->T(idIN) = stressIN;
    return;
}

void SlurryGranularPhase::assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    MaterialTensor tmp = body->points->T[idIN];
    body->points->T[idIN] = tmp - (1.0/3.0 * tmp.trace() + pressureIN)*MaterialTensor::Identity();
    return;
}

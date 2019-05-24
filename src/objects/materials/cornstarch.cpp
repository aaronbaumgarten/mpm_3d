//
// Created by aaron on 6/14/18.
// cornstarch.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <dlfcn.h>
#include <objects/grids/grids.hpp>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "materials.hpp"

bool CORNSTARCH_DEBUG = false;

/*----------------------------------------------------------------------------*/
//
void Cornstarch::init(Job* job, Body* body){
    if (fp64_props.size() < 18){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 18 properties defined ({E, nu, mu_1, a_0, a_inf, K_3, K_4, K_5, K_6, phi_j, phi_c, Delta, tau_star, alpha, grains_rho, eta_0, grains_d, fluid_rho}, {fluid_body}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        E = fp64_props[0];
        nu = fp64_props[1];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        mu_1 = fp64_props[2];
        mu_2 = mu_1;
        I_0 = 0.0;
        a_0 = fp64_props[3];
        a_inf = fp64_props[4];
        K_3 = fp64_props[5];
        K_4 = fp64_props[6];
        K_5 = fp64_props[7];
        K_6 = fp64_props[8];

        if (fp64_props.size() > 18){
            K_7 = fp64_props[18];
        } else {
            K_7 = 0;
        }

        phi_j = fp64_props[9];
        phi_c = fp64_props[10];
        Delta = fp64_props[11];
        tau_star = fp64_props[12];

        alpha = fp64_props[13];

        phi_star = phi_j + (phi_c - phi_j)*Delta;

        grains_rho = fp64_props[14];
        eta_0 = fp64_props[15];
        grains_d = fp64_props[16];
        fluid_rho = fp64_props[17];

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
            std::cout << std::endl << "WARNING: No fluid body defined! Setting eta to 0!" << std::endl << std::endl;
        }

        if (K_3 > mu_1 / phi_m){
            std::cout << std::endl;
            std::cout << "WARNING: K_3 = " << K_3 << " > mu_1/phi_m = " << mu_1 / phi_m << "! For dilute suspensions, this may result in negative strength!" << std::endl;
            std::cout << std::endl;
        }

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, mu_1 = %g, a_0 = %g, a_inf = %g, K_3 = %g, K_4 = %g, K_5 = %g, K_6 = %g, phi_j = %g, phi_c = %g, Delta = %g, tau_star = %g, alpha = %g, grains_rho = %g, eta_0 = %g, grains_d = %g, fluid_rho = %g).\n",
               E, nu, G, K, mu_1, a_0, a_inf, K_3, K_4, K_5, K_6, phi_j, phi_c, Delta, tau_star, alpha, grains_rho, eta_0, grains_d, fluid_rho);
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

    phi_m_vec.resize(body->points->x.size());
    phi_m_vec.setZero();
    c.resize(body->points->x.size());
    c.setZero();

    //if c0 given, then initialize c with value as assigned
    if (fp64_props.size() > 19){
        for (int i=0; i<c.rows(); i++){
            c(i) = fp64_props[19];
        }
    }

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//
void Cornstarch::writeFrame(Job* job, Body* body, Serializer* serializer){
    SlurryGranularPhase::writeFrame(job,body,serializer);
    serializer->writeScalarArray(phi_m_vec,"phi_m");
    serializer->writeScalarArray(c,"c");

    //post-process f field
    Eigen::VectorXd nvec(body->nodes->x.size());
    Eigen::VectorXd vi(body->nodes->x.size());
    Eigen::VectorXd pvec(body->points->x.size());
    for (int i=0; i<pvec.rows(); i++){
        pvec(i) = c(i)*body->points->v(i);
    }
    vi = body->S*body->points->v; //map volume to nodes
    nvec = body->S*pvec; //map weighted f to nodes
    for (int i=0; i<nvec.rows(); i++){
        if (nvec(i)>0){
            nvec(i) /= vi(i); //normalize by weight
        }
    }
    pvec = body->S.operate(nvec,MPMSparseMatrixBase::TRANSPOSED);
    serializer->writeScalarArray(pvec,"f");

    return;
}

/*----------------------------------------------------------------------------*/
//
double Cornstarch::getBeta(double phi, double phi_eq){
    return (K_3 * (phi - phi_eq));
}

/*----------------------------------------------------------------------------*/
//
void Cornstarch::calcState(double &gdp, double &p, double &eta_in, double &phi_in, double &I_out, double &Iv_out, double &Im_out, double &mu_out, double &phi_eq, double &beta_out){
    I_out = gdp * grains_d * std::sqrt(grains_rho/p);
    Iv_out = gdp * eta_in / p;
    Im_out = std::sqrt(I_out*I_out + 2 * Iv_out);

    if (Im_out == 0){
        mu_out = mu_1;
    } else if (p <= 0){
        mu_out = mu_1;
    } else {
        mu_out = mu_1 + 2.5 * phi_in * Iv_out/(a*Im_out);
    }

    if (p <= 0){
        phi_eq = 0;
    } else {
        phi_eq = phi_m / (1 + a * Im_out);
    }

    beta_out = (K_3 * (phi_in - phi_eq));

    return;
}

/*----------------------------------------------------------------------------*/
//
void Cornstarch::calculateStress(Job* job, Body* body, int SPEC){
    MaterialTensor T, T_tr, L, D, W;
    MaterialTensor tmpMat;

    double H, S, t_k, c_m;

    double trD, tau_bar, tau_bar_tr, p, p_tr;
    double beta, mu, phi_eq, xi_dot_1, xi_dot_2;
    double tau_bar_k, p_k, tau_bar_kplus, p_kplus;
    double I_v_tr, I_tr, I_m_tr, gammap_dot_tr, xi_dot_2_tr;
    double tmpVal;

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

        /*************************************************/
        //fix nan clumpiness
        if (!std::isfinite(c(i))){
            c(i) = 0;
        }
        //assign phi_m based off of clumpiness
        /*
        if (phi(i) > phi_c && phi(i) < phi_star) {
            phi_m = phi_j + (phi(i) - phi_j) * c(i);
        } else if (phi(i) < phi_c){
            phi_m = phi_j + (phi_c - phi_j) * c(i);
        } else {
            phi_m = phi_j + (phi_star - phi_j) * c(i);
        }
        */
        //assign c_m based off of phi_star
        if (phi(i) > phi_c && phi(i) < phi_star) {
            c_m = (phi(i) - phi_j)/(phi_c - phi_j);
        } else if (phi(i) < phi_c){
            c_m = 1;
        } else {
            c_m = (phi_star - phi_j)/(phi_c - phi_j);
        }

        //assign a based off of clumpiness
        a = a_0 + (a_inf - a_0)*c(i);
        phi_m = phi_j + (phi_c - phi_j)*c(i);
        /*************************************************/

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
                if (k > 100 && CORNSTARCH_DEBUG){
                    std::cout << "f1 weak: " << p_k << ", " << r_p << ", " << dr_dp << std::endl;
                }

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

                if (k > 5 && std::abs(p_k - p_kplus) < ABS_TOL){
                    //std::cout << "Maximum iterations exceeded; Newton scheme appears to have stagnated. Exiting loop. " << r_p_tr << std::endl;
                    p_k = p_kplus;
                    r_p = r_p_tr;
                    break;
                }

                p_k = p_kplus;
                r_p = r_p_tr;
            }

            //check that zero strength assumption is true and admissible
            gammap_dot_tr = tau_bar_tr / (G*job->dt);
            xi_dot_2_tr = 0;
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

            k = 0;
            //newton method
            while (r.norm() > b.norm() * REL_TOL && r.norm() > ABS_TOL) {
                k+=1;
                if (k>100 && CORNSTARCH_DEBUG){
                    std::cout << "f1: " << p_k << ", " << tau_bar_k << ", " << r.norm() << ", " << dr.norm() << std::endl;
                }
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

                if (k > 5 && std::abs(p_k - p_kplus) < ABS_TOL && std::abs(tau_bar_k - tau_bar_kplus) < ABS_TOL){
                    //std::cout << "Maximum iterations exceeded; Newton scheme appears to have stagnated. Exiting loop. " << r_tr.norm() << std::endl;
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
            xi_dot_2_tr = 0;
            if ((phi(i) >= phi_m) or (p_k <= (a*a*phi(i)*phi(i))/((phi_m - phi(i))*(phi_m - phi(i)))*(gammap_dot_tr*gammap_dot_tr*grains_d*grains_d*grains_rho + 2.0*eta(i)*gammap_dot_tr) )){
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
                if (k > 100 && CORNSTARCH_DEBUG){
                    std::cout << "f1,f3 weak: " << p_k << ", " << r_p << ", " << dr_dp << std::endl;
                }

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

                if (k > 5 && std::abs(p_k - p_kplus) < ABS_TOL){
                    //std::cout << "Maximum iterations exceeded; Newton scheme appears to have stagnated. Exiting loop. " << r_p_tr << std::endl;
                    p_k = p_kplus;
                    r_p = r_p_tr;
                    break;
                }

                p_k = p_kplus;
                r_p = r_p_tr;
            }

            //check that zero strength assumption is true
            gammap_dot_tr = tau_bar_tr / (G*job->dt);
            xi_dot_2_tr = 0;
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
                if (k>100 && CORNSTARCH_DEBUG){
                    //std::cout << "f1,f3: " << p_k << ", " << tau_bar_k << ", " << r.norm() << ", " << dr.norm() << std::endl;
                    std::cout << "f1,f3: " << p_k << ", " << tau_bar_k << std::endl;
                    std::cout << "       " << r(0) << ", " << r(1) << std::endl;
                    std::cout << "       " << dr(0,0) << ", " << dr(0,1) << std::endl;
                    std::cout << "       " << dr(1,0) << ", " << dr(1,1) << std::endl;
                    std::cout << "       " << i << ", " << phi(i) << ", " << ", " << p_tr << ", " << tau_bar_tr << std::endl;
                    std::cout << "       " << lambda_tmp << ", " << upd(0) << ", " << upd(1) << std::endl << std::endl;

                }
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

                if (k > 5 && std::abs(p_k - p_kplus) < ABS_TOL && std::abs(tau_bar_k - tau_bar_kplus) < ABS_TOL){
                    //std::cout << "Maximum iterations exceeded; Newton scheme appears to have stagnated. Exiting loop. " << r_tr.norm() << std::endl;
                    p_k = p_kplus;
                    tau_bar_k = tau_bar_kplus;
                    r = r_tr;
                    //std::cout << "u, " << r.transpose() << std::endl;
                    //std::cout << dr << std::endl;
                    //std::cout << xi_dot_2 << std::endl;
                    break;
                }

                p_k = p_kplus;
                tau_bar_k = tau_bar_kplus;
                r = r_tr;
            }

            //this is the last possible case, if we got here then this is the answer
            gammap_dot_tr = (tau_bar_tr - tau_bar_k) / (G*job->dt);
            xi_dot_2_tr = (p_k - p_tr)/(K*job->dt) - beta*gammap_dot_tmp;
            p = p_k;
            tau_bar = tau_bar_k;
            is_solved = true;
        }

        /*if (!std::isfinite(p) || !std::isfinite(tau_bar)){
            std::cout << "u";
        }*/

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

            //update history
            phi_m_vec(i) = phi_m;

            //K_6 has units of (Pa s)^(-1)
            H = K_5 * (tau_bar/tau_star) * std::sqrt(tau_bar/tau_star) * gammap_dot(i);
            S = K_5 * (gammap_dot(i) + (K_6 * std::pow(phi_j - phi(i), alpha) + K_7 * std::pow(phi_j - phi(i), 0.1)) * tau_bar);

            if (phi(i) > phi_j){
                c(i) = 1;
            } else {
                c(i) = (c(i) + H*c_m*job->dt)/(1 + H*job->dt + S*job->dt);
            }
        }

    }

    return;
}

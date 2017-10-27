//
// Created by aaron on 7/14/17.
// slurry_mu_iv.cpp
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <Eigen/Dense>

#include "job.hpp"
#include "serializer.hpp"

#include "body.hpp"
#include "points.hpp"
#include "material.hpp"

//residual tolerance
double ABS_TOL = 1e-3;
double REL_TOL = 1e-6;
double h = 1e-5;

//standard linear elastic terms
double E, nu, G, K;
double lambda;

//mu(Im)
//0.32, 0.7, (?), (?), (?), 0.585, (?), (?)
//0.26, 1.00
double mu_1, mu_2, K_3, K_4, K_5, phi_m, grains_rho, eta_0, grain_diam, fluid_rho;
//0.707, 0.233
double a, I_0;

//sand
Eigen::VectorXd gammap(0,1);
Eigen::VectorXd gammap_dot(0,1);
Eigen::VectorXd I_v(0,1);
Eigen::VectorXd I(0,1);
Eigen::VectorXd I_m(0,1);
Eigen::VectorXd phi(0,1);
Eigen::VectorXd eta(0,1);
int fluid_body_id = -1;


extern "C" void materialWriteFrame(Job* job, Body* body, Serializer* serializer);
extern "C" std::string materialSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
extern "C" int materialLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);

extern "C" void materialInit(Job* job, Body* body);
extern "C" void materialCalculateStress(Job* job, Body* body, int SPEC);
extern "C" void materialAssignStress(Job* job, Body* body, Eigen::MatrixXd stressIN, int idIN, int SPEC);
extern "C" void materialAssignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

/*----------------------------------------------------------------------------*/

void materialInit(Job* job, Body* body){
    if (body->material.fp64_props.size() < 14 || body->material.str_props.size() < 1){
        std::cout << body->material.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 13 properties defined ({E, nu, mu_1, mu_2, a, b, K_3, K_4, K_5, phi_m, grains_rho, eta_0, grain_diam, fluid_rho}, {fluid_body}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        mu_1 = body->material.fp64_props[2];
        mu_2 = body->material.fp64_props[3];
        a = body->material.fp64_props[4];
        I_0 = body->material.fp64_props[5];
        K_3 = body->material.fp64_props[6];
        K_4 = body->material.fp64_props[7];
        K_5 = body->material.fp64_props[8];
        phi_m = body->material.fp64_props[9];
        grains_rho = body->material.fp64_props[10];
        eta_0 = body->material.fp64_props[11];
        grain_diam = body->material.fp64_props[12];
        fluid_rho = body->material.fp64_props[13];

        //set body id by name
        if (body->material.str_props.size() >= 1){
            for (size_t b = 0; b < job->bodies.size(); b++) {
                if (body->material.str_props[0].compare(job->bodies[b].name) == 0){
                    fluid_body_id = b;
                    break;
                }
            }
        }

        if (fluid_body_id == -1){
            std::cout << std::endl << "WARNING: No fluid body defined! Setting eta to 0!" << std::endl << std::endl;
        }

        if (K_4 > mu_1 / phi_m){
            std::cout << std::endl;
            std::cout << "WARNING: K_4 = " << K_4 << " > mu_1/phi_m = " << mu_1 / phi_m << "! For dilute suspensions, this may result in negative strength!" << std::endl;
            std::cout << std::endl;
        }

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, mu_1 = %g, mu_2 = %g, a = %g, I_0 = %g K_3 = %g, K_4 = %g, K_5 = %g, phi_m = %g, grains_rho = %g, eta_0 = %g, grain_diam = %g, fluid_rho = %g).\n",
               E, nu, G, K, mu_1, mu_2, a, I_0, K_3, K_4, K_5, phi_m, grains_rho, eta_0, grain_diam, fluid_rho);
    }

    gammap_dot.resize(body->points.x.rows());
    gammap_dot.setZero();
    gammap.resize(body->points.x.rows());
    gammap.setZero();
    I_v.resize(body->points.x.rows());
    I_v.setZero();
    I.resize(body->points.x.rows());
    I.setZero();
    I_m.resize(body->points.x.rows());
    I_m.setZero();
    phi.resize(body->points.x.rows());
    phi.setZero();
    eta.resize(body->points.x.rows());
    eta.setZero();

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void materialWriteFrame(Job* job, Body* body, Serializer* serializer) {
    //nothing to report
    serializer->serializerWriteScalarArray(gammap,"gammap");
    serializer->serializerWriteScalarArray(gammap_dot,"gammap_dot");
    serializer->serializerWriteScalarArray(I_v,"I_v");
    serializer->serializerWriteScalarArray(I,"I");
    serializer->serializerWriteScalarArray(I_m,"I_m");
    serializer->serializerWriteScalarArray(phi,"packing_fraction");
    serializer->serializerWriteScalarArray(eta,"eta_interp");
    return;
}

/*----------------------------------------------------------------------------*/

std::string materialSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename
    std::stringstream s;
    s << "mpm_v2."  << body->name << "." << body->id << ".material." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    if (ffile.is_open()){
        ffile << "# mpm_v2 materials/isolin.so\n";
        ffile << E << "\n";
        ffile << nu << "\n";
        ffile << mu_1 << "\n";
        ffile << mu_2 << "\n";
        ffile << a << "\n";
        ffile << I_0 << "\n";
        ffile << K_3 << "\n";
        ffile << K_4 << "\n";
        ffile << K_5 << "\n";
        ffile << phi_m << "\n";
        ffile << grains_rho << "\n";
        ffile << eta_0 << "\n";
        ffile << grain_diam << "\n";
        ffile << fluid_rho << "\n";
        ffile << fluid_body_id << "\n";

        ffile << gammap_dot.rows() << "\n";

        ffile << "gammap_dot\n{\n";
        job->jobScalarArrayToFile(gammap_dot,ffile);
        ffile << "}\n";
        ffile << "gammap\n{\n";
        job->jobScalarArrayToFile(gammap,ffile);
        ffile << "}\n";
        ffile << "I_v\n{\n";
        job->jobScalarArrayToFile(I_v,ffile);
        ffile << "}\n";

        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filepath+filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Material Saved: [" << body->name << "]." << std::endl;

    return filename;
}

/*----------------------------------------------------------------------------*/

int materialLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    std::string line;
    std::ifstream fin(fullpath);

    if(fin.is_open()){
        std::getline(fin,line); //first line
        std::getline(fin,line); //E
        E = std::stod(line);
        std::getline(fin,line); //nu
        nu = std::stod(line);
        std::getline(fin,line); //mu_1
        mu_1 = std::stod(line);
        std::getline(fin,line); //mu_2
        mu_2 = std::stod(line);
        std::getline(fin,line); //a
        a = std::stod(line);
        std::getline(fin,line); //b
        I_0 = std::stod(line);
        std::getline(fin,line); //K_3
        K_3 = std::stod(line);
        std::getline(fin,line); //K_4
        K_4 = std::stod(line);
        std::getline(fin,line); //K_5
        K_5 = std::stod(line);
        std::getline(fin,line); //phi_m
        phi_m = std::stod(line);
        std::getline(fin,line); //grains_rho
        grains_rho = std::stod(line);
        std::getline(fin,line); //eta_0
        eta_0 = std::stod(line);
        std::getline(fin,line); //grain_diam
        grain_diam = std::stod(line);
        std::getline(fin,line); //fluid_rho
        fluid_rho = std::stod(line);
        std::getline(fin,line);//fluid body id
        fluid_body_id = std::stoi(line);

        std::getline(fin,line); //len
        int len = std::stoi(line);

        gammap.resize(len);
        gammap.setZero();
        gammap_dot.resize(len);
        gammap_dot.setZero();
        I_v.resize(len);
        I_v.setZero();
        phi.resize(len);
        phi.setZero();

        while(std::getline(fin,line)){
            if (line.compare("gammap_dot") == 0){
                std::getline(fin,line); //{
                job->jobScalarArrayFromFile(gammap_dot,fin);
                std::getline(fin,line); //}
            } else if (line.compare("gammap") == 0){
                std::getline(fin,line);
                job->jobScalarArrayFromFile(gammap,fin);
                std::getline(fin,line);
            } else if (line.compare("I_v") == 0){
                std::getline(fin,line);
                job->jobScalarArrayFromFile(I_v,fin);
                std::getline(fin,line);
            } else {
                std::cerr << "Unknown field title: " << line << std::endl;
            }
        }

        fin.close();

        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, mu_1 = %g, mu_2 = %g, a = %g, I_0 = %g K_3 = %g, K_4 = %g, K_5 = %g, phi_m = %g, grains_rho = %g, eta_0 = %g, grain_diam = %g, fluid_rho = %g).\n",
               E, nu, G, K, mu_1, mu_2, a, I_0, K_3, K_4, K_5, phi_m, grains_rho, eta_0, grain_diam, fluid_rho);
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Material Loaded: [" << body->name << "]." << std::endl;
    return 1;
}

/*----------------------------------------------------------------------------*/

void materialCalculateStress(Job* job, Body* body, int SPEC){
    Eigen::MatrixXd T = job->jobTensor<double>();
    Eigen::MatrixXd T_tr = job->jobTensor<double>();
    Eigen::MatrixXd L = job->jobTensor<double>();
    Eigen::MatrixXd D = job->jobTensor<double>();
    Eigen::MatrixXd D_0 = job->jobTensor<double>();
    Eigen::MatrixXd W = job->jobTensor<double>();

    double trD, tau_bar, tau_bar_tr, p, p_tr, p_plus, tau_bar_plus;
    double beta, mu, phi_eq, xi_dot_1, xi_dot_2;
    double mu_f;
    double I_v_tr, I_tr, I_m_tr, gammap_dot_tr;

    Eigen::MatrixXd tmpMat = job->jobTensor<double>();
    Eigen::VectorXd tmpVec;

    Eigen::Vector2d r, b, upd;
    Eigen::Matrix2d dr;

    double p_tmp, tau_bar_tmp, lambda_tmp;
    Eigen::Vector2d r_tmp;

    //calculate packing fraction
    phi = 1/grains_rho * (body->points.m.array() / body->points.v.array());

    //interpolate eta from fluid points
    if (eta_0 > 0) {
        Eigen::VectorXd pvec = job->bodies[fluid_body_id].points.m * eta_0;
        Eigen::VectorXd nvec(job->bodies[fluid_body_id].nodes.m.rows());
        job->bodies[fluid_body_id].bodyCalcNodalValues(job, nvec, pvec, Body::SET);
        for (size_t i = 0; i < nvec.rows(); i++) {
            if (job->bodies[fluid_body_id].nodes.m(i) > 0) {
                nvec(i) = nvec(i) / job->bodies[fluid_body_id].nodes.m(i);
            } else {
                nvec(i) = 0;
            }
        }
        body->bodyCalcPointValues(job, eta, nvec, Body::SET);
    } else {
        eta.setZero();
    }

    for (size_t i=0;i<body->points.x.rows();i++) {
        if (body->points.active[i] == 0) {
            continue;
        }

        tmpVec = body->points.L.row(i).transpose();
        L = job->jobTensor<double>(tmpVec.data());

        tmpVec = body->points.T.row(i).transpose();
        T = job->jobTensor<double>(tmpVec.data());

        D = 0.5 * (L + L.transpose());
        W = 0.5 * (L - L.transpose());

        trD = D.trace();
        D_0 = D - trD / D.rows() * job->jobTensor<double>(Job::IDENTITY);

        //trial stress
        T_tr = T + job->dt * ((2 * G * D) + (lambda * trD * job->jobTensor<double>(Job::IDENTITY)) + (W * T) - (T * W));
        tau_bar_tr = (T_tr - (T_tr.trace() / T_tr.rows()) * job->jobTensor<double>(Job::IDENTITY)).norm() / std::sqrt(2.0);
        p_tr = -T_tr.trace() / T_tr.rows();

        //check f2 yield surface
        if (phi(i) > phi_m) {
            beta = K_3 * (phi(i) - phi_m) + K_4 * (phi(i));
        } else {
            beta = K_4 * phi(i);
        }
        if (p_tr <= 0 && tau_bar_tr == 0){
            tau_bar = 0;
            p = 0;
            gammap_dot_tr = 0;
        } else if(p_tr + K*(tau_bar_tr/G)*beta <= 0) {
            tau_bar = 0;
            p = 0;
            gammap_dot_tr = tau_bar_tr / (G * job->dt);
        } else {

            //check for admissible stress state
            if (phi(i) > phi_m) {
                beta = (K_3 + K_4) * (phi(i) - phi_m);
            } else {
                beta = K_4 * (phi(i) - phi_m);
            }
            if (p_tr >= 0 && tau_bar_tr < (mu_1 + beta) * p_tr) {
                tau_bar = tau_bar_tr;
                p = p_tr;
                gammap_dot_tr = 0;
            } else {

                //check f1 only
                tau_bar = 0; //[0, tau_bar_tr]
                p = 0; //[0,?]
                if (phi(i) > phi_m) {
                    beta = (K_3 + K_4) * (phi(i) - phi_m);
                } else {
                    beta = K_4 * (phi(i) - phi_m);
                }
                r(0) = tau_bar_tr;
                r(1) = -p_tr - K * (tau_bar_tr / G) * beta;
                b = r;
                while (r.norm() > b.norm() * REL_TOL && r.norm() > ABS_TOL) {
                    gammap_dot_tr = (tau_bar_tr - tau_bar) / (G * job->dt);

                    if (p == 0) {
                        mu = mu_2;
                        phi_eq = 0;
                    } else {
                        I_v_tr = eta(i) * gammap_dot_tr / p;
                        I_tr = gammap_dot_tr * grain_diam * std::sqrt(grains_rho/p);
                        I_m_tr = std::sqrt(I_tr*I_tr + 2*I_v_tr);
                        mu = mu_1 + (mu_2 - mu_1) / (1 + I_0 / I_m_tr) + 5.0 / 2.0 * phi(i) * I_v_tr / (a*I_m_tr);
                        phi_eq = phi_m / (1 + a*I_m_tr);
                    }
                    if (gammap_dot_tr == 0 || eta(i) == 0) {
                        mu = mu_1;
                    }

                    if (phi(i) > phi_m) {
                        beta = K_3 * (phi(i) - phi_m) + K_4 * (phi(i) - phi_eq);
                    } else {
                        beta = K_4 * (phi(i) - phi_eq);
                    }

                    mu_f = mu + beta;
                    if (mu_f < 0){
                        mu_f = 0;
                    }

                    r(0) = tau_bar - mu_f * p;
                    r(1) = p - p_tr - K * job->dt * (beta * gammap_dot_tr);

                    if (r.norm() > b.norm() * REL_TOL && r.norm() > ABS_TOL) {
                        if (p > 0) {
                            p_plus = p * (1 + h);
                        } else {
                            p_plus = h;
                        }

                        gammap_dot_tr = (tau_bar_tr - tau_bar) / (G * job->dt);

                        if (p_plus == 0) {
                            mu = mu_2;
                            phi_eq = 0;
                        } else {
                            I_v_tr = eta(i) * gammap_dot_tr / p_plus;
                            I_tr = gammap_dot_tr * grain_diam * std::sqrt(grains_rho/p_plus);
                            I_m_tr = std::sqrt(I_tr*I_tr + 2*I_v_tr);
                            mu = mu_1 + (mu_2 - mu_1) / (1 + I_0 / I_m_tr) + 5.0 / 2.0 * phi(i) * I_v_tr / (a*I_m_tr);
                            phi_eq = phi_m / (1 + a*I_m_tr);
                        }
                        if (gammap_dot_tr == 0 || eta(i) == 0) {
                            mu = mu_1;
                        }

                        if (phi(i) > phi_m) {
                            beta = K_3 * (phi(i) - phi_m) + K_4 * (phi(i) - phi_eq);
                        } else {
                            beta = K_4 * (phi(i) - phi_eq);
                        }

                        mu_f = mu + beta;
                        if (mu_f < 0){
                            mu_f = 0;
                        }

                        dr(0, 0) = (tau_bar - mu_f * p_plus - r(0)) / (p_plus - p);
                        dr(1, 0) = (p_plus - p_tr - K * job->dt * (beta * gammap_dot_tr) - r(1)) / (p_plus - p);

                        //--------------

                        if (tau_bar > 0 && tau_bar * (1 + h) < tau_bar_tr) {
                            tau_bar_plus = tau_bar * (1 + h);
                        } else if (tau_bar * (1 + h) >= tau_bar_tr) {
                            tau_bar_plus = tau_bar * (1 - h);
                        } else {
                            tau_bar_plus = h;
                        }
                        gammap_dot_tr = (tau_bar_tr - tau_bar_plus) / (G * job->dt);

                        if (p == 0) {
                            mu = mu_2;
                            phi_eq = 0;
                        } else {
                            I_v_tr = eta(i) * gammap_dot_tr / p;
                            I_tr = gammap_dot_tr * grain_diam * std::sqrt(grains_rho/p);
                            I_m_tr = std::sqrt(I_tr*I_tr + 2*I_v_tr);
                            mu = mu_1 + (mu_2 - mu_1) / (1 + I_0 / I_m_tr) + 5.0 / 2.0 * phi(i) * I_v_tr / (a*I_m_tr);
                            phi_eq = phi_m / (1 + a*I_m_tr);
                        }
                        if (gammap_dot_tr == 0 || eta(i) == 0) {
                            mu = mu_1;
                        }

                        if (phi(i) > phi_m) {
                            beta = K_3 * (phi(i) - phi_m) + K_4 * (phi(i) - phi_eq);
                        } else {
                            beta = K_4 * (phi(i) - phi_eq);
                        }

                        mu_f = mu + beta;
                        if (mu_f < 0){
                            mu_f = 0;
                        }

                        dr(0, 1) = (tau_bar_plus - mu_f * p - r(0)) / (tau_bar_plus - tau_bar);
                        dr(1, 1) = (p - p_tr - K * job->dt * (beta * gammap_dot_tr) - r(1)) / (tau_bar_plus - tau_bar);
                    } else {
                        break;
                    }

                    upd = dr.inverse()*r;
                    lambda_tmp = 1;
                    //p = p - upd(0);
                    //tau_bar = tau_bar - upd(1);
                    do {
                        p_tmp = p - lambda_tmp*upd(0);
                        tau_bar_tmp = tau_bar - lambda_tmp*upd(1);

                        if (p_tmp < 0) {
                            p_tmp = 0;
                        }

                        if (tau_bar_tmp < 0) {
                            tau_bar_tmp = 0;
                        } else if (tau_bar_tmp > tau_bar_tr) {
                            tau_bar_tmp = tau_bar_tr;
                        }

                        gammap_dot_tr = (tau_bar_tr - tau_bar_tmp) / (G * job->dt);

                        if (p_tmp == 0) {
                            mu = mu_2;
                            phi_eq = 0;
                        } else {
                            I_v_tr = eta(i) * gammap_dot_tr / p_tmp;
                            I_tr = gammap_dot_tr * grain_diam * std::sqrt(grains_rho/p_tmp);
                            I_m_tr = std::sqrt(I_tr*I_tr + 2*I_v_tr);
                            mu = mu_1 + (mu_2 - mu_1) / (1 + I_0 / I_m_tr) + 5.0 / 2.0 * phi(i) * I_v_tr / (a*I_m_tr);
                            phi_eq = phi_m / (1 + a*I_m_tr);
                        }
                        if (gammap_dot_tr == 0 || eta(i) == 0) {
                            mu = mu_1;
                        }

                        if (phi(i) > phi_m) {
                            beta = K_3 * (phi(i) - phi_m) + K_4 * (phi(i) - phi_eq);
                        } else {
                            beta = K_4 * (phi(i) - phi_eq);
                        }

                        mu_f = mu + beta;
                        if (mu_f < 0){
                            mu_f = 0;
                        }

                        r_tmp(0) = tau_bar_tmp - mu_f * p_tmp;
                        r_tmp(1) = p_tmp - p_tr - K * job->dt * (beta * gammap_dot_tr);

                        lambda_tmp *= 0.5;
                    } while (r_tmp.norm() > r.norm() && lambda_tmp > REL_TOL);

                    p = p_tmp;
                    tau_bar = tau_bar_tmp;

                    if (p < 0) {
                        p = 0;
                    }

                    if (tau_bar < 0) {
                        tau_bar = 0;
                    } else if (tau_bar > tau_bar_tr) {
                        tau_bar = tau_bar_tr;
                    }

                    if (false) {
                        std::cout << i << " : " << r.norm() << " <? " << b.norm() << " : " << tau_bar << " : " << p << std::endl;
                        std::cout << phi_eq << ", " << phi(i) << std::endl;
                        std::cout << beta << ", " << gammap_dot_tr << std::endl;
                        std::cout << tau_bar_tr << " : " << p_tr << std::endl;
                        std::cout << tau_bar_plus << " : " << p_plus << std::endl;
                        std::cout << dr << std::endl;
                        std::cout << r << std::endl << std::endl;
                        std::cout << upd << std::endl << std::endl;
                    }
                    //std::cout << T << std::endl << L << std::endl << std::endl;


                    if (!std::isfinite(p) || !std::isfinite(tau_bar)){ // || (mu+beta) < 0){
                        std::cout << "u";
                    }
                }

                if (phi(i) > phi_m || (p - (a * phi(i) / (phi_m - phi(i)))*(a * phi(i) / (phi_m - phi(i))) * (gammap_dot_tr*gammap_dot_tr * grain_diam * grain_diam * grains_rho + 2 * eta(i) * (gammap_dot_tr))) < 0){
                    //do nothing
                } else {
                    //check f1,f3
                    if (phi(i) > phi_m) {
                        beta = K_3 * (phi(i) - phi_m) + K_4 * (phi(i));
                    } else {
                        beta = K_4 * phi(i);
                    }
                    r(0) = tau_bar_tr;
                    r(1) = -p_tr - K * (tau_bar_tr / G) * beta;
                    b = r;
                    while (r.norm() > b.norm() * REL_TOL && r.norm() > ABS_TOL) {
                        gammap_dot_tr = (tau_bar_tr - tau_bar) / (G * job->dt);
                        if (p == 0) {
                            mu = mu_2;
                            phi_eq = 0;
                        } else {
                            I_v_tr = eta(i) * gammap_dot_tr / p;
                            I_tr = gammap_dot_tr * grain_diam * std::sqrt(grains_rho/p);
                            I_m_tr = std::sqrt(I_tr*I_tr + 2*I_v_tr);
                            mu = mu_1 + (mu_2 - mu_1) / (1 + I_0 / I_m_tr) + 5.0 / 2.0 * phi(i) * I_v_tr / (a*I_m_tr);
                            phi_eq = phi_m / (1 + a*I_m_tr);
                        }
                        if (gammap_dot_tr == 0 || eta(i) == 0) {
                            mu = mu_1;
                        }

                        if (phi(i) > phi_m) {
                            beta = K_3 * (phi(i) - phi_m) + K_4 * (phi(i) - phi_eq);
                        } else {
                            beta = K_4 * (phi(i) - phi_eq);
                        }
                        xi_dot_2 = (p - p_tr) / (K * job->dt) - beta * gammap_dot_tr;

                        mu_f = mu + beta;
                        if (mu_f < 0){
                            mu_f = 0;
                        }

                        r(0) = tau_bar - mu_f * p;
                        //r(1) = p - (phi(i) / (phi_m - phi(i))) * (phi(i) / (phi_m - phi(i))) * eta(i) * (gammap_dot_tr - K_5 * xi_dot_2);
                        r(1) = p - (a * phi(i) / (phi_m - phi(i)))*(a * phi(i) / (phi_m - phi(i))) * ((gammap_dot_tr - K_5 * xi_dot_2)*(gammap_dot_tr - K_5 * xi_dot_2) * grain_diam * grain_diam * grains_rho + 2 * eta(i) * ((gammap_dot_tr - K_5 * xi_dot_2)));

                        if (r.norm() > b.norm() * REL_TOL && r.norm() > ABS_TOL) {
                            if (p > 0) {
                                p_plus = p * (1 + h);
                            } else {
                                p_plus = h;
                            }
                            gammap_dot_tr = (tau_bar_tr - tau_bar) / (G * job->dt);

                            if (p_plus == 0) {
                                mu = mu_2;
                                phi_eq = 0;
                            } else {
                                I_v_tr = eta(i) * gammap_dot_tr / p_plus;
                                I_tr = gammap_dot_tr * grain_diam * std::sqrt(grains_rho/p_plus);
                                I_m_tr = std::sqrt(I_tr*I_tr + 2*I_v_tr);
                                mu = mu_1 + (mu_2 - mu_1) / (1 + I_0 / I_m_tr) + 5.0 / 2.0 * phi(i) * I_v_tr / (a*I_m_tr);
                                phi_eq = phi_m / (1 + a*I_m_tr);
                            }
                            if (gammap_dot_tr == 0 || eta(i) == 0) {
                                mu = mu_1;
                            }

                            if (phi(i) > phi_m) {
                                beta = K_3 * (phi(i) - phi_m) + K_4 * (phi(i) - phi_eq);
                            } else {
                                beta = K_4 * (phi(i) - phi_eq);
                            }
                            xi_dot_2 = (p_plus - p_tr) / (K * job->dt) - beta * gammap_dot_tr;

                            mu_f = mu + beta;
                            if (mu_f < 0){
                                mu_f = 0;
                            }

                            dr(0, 0) = tau_bar - mu_f * p_plus;
                            dr(0, 0) = (dr(0, 0) - r(0)) / (p_plus - p);
                            //dr(1, 0) = p_plus - (phi(i) / (phi_m - phi(i))) * (phi(i) / (phi_m - phi(i))) * eta(i) * (gammap_dot_tr - K_5 * xi_dot_2);
                            dr(1, 0) = p_plus - (a * phi(i) / (phi_m - phi(i)))*(a * phi(i) / (phi_m - phi(i))) * ((gammap_dot_tr - K_5 * xi_dot_2)*(gammap_dot_tr - K_5 * xi_dot_2) * grain_diam * grain_diam * grains_rho + 2 * eta(i) * ((gammap_dot_tr - K_5 * xi_dot_2)));
                            dr(1, 0) = (dr(1, 0) - r(1)) / (p_plus - p);

                            //--------------

                            if (tau_bar > 0 && tau_bar * (1 + h) < tau_bar_tr) {
                                tau_bar_plus = tau_bar * (1 + h);
                            } else if (tau_bar * (1 + h) >= tau_bar_tr) {
                                tau_bar_plus = tau_bar * (1 - h);
                            } else {
                                tau_bar_plus = h;
                            }
                            gammap_dot_tr = (tau_bar_tr - tau_bar_plus) / (G * job->dt);

                            if (p == 0) {
                                mu = mu_2;
                                phi_eq = 0;
                            } else {
                                I_v_tr = eta(i) * gammap_dot_tr / p;
                                I_tr = gammap_dot_tr * grain_diam * std::sqrt(grains_rho/p);
                                I_m_tr = std::sqrt(I_tr*I_tr + 2*I_v_tr);
                                mu = mu_1 + (mu_2 - mu_1) / (1 + I_0 / I_m_tr) + 5.0 / 2.0 * phi(i) * I_v_tr / (a*I_m_tr);
                                phi_eq = phi_m / (1 + a*I_m_tr);
                            }
                            if (gammap_dot_tr == 0 || eta(i) == 0) {
                                mu = mu_1;
                            }

                            if (phi(i) > phi_m) {
                                beta = K_3 * (phi(i) - phi_m) + K_4 * (phi(i) - phi_eq);
                            } else {
                                beta = K_4 * (phi(i) - phi_eq);
                            }
                            xi_dot_2 = (p - p_tr) / (K * job->dt) - beta * gammap_dot_tr;


                            mu_f = mu + beta;
                            if (mu_f < 0){
                                mu_f = 0;
                            }

                            dr(0, 1) = tau_bar_plus - mu_f * p;
                            dr(0, 1) = (dr(0, 1) - r(0)) / (tau_bar_plus - tau_bar);
                            //dr(1, 1) = p - (phi(i) / (phi_m - phi(i)))*(phi(i) / (phi_m - phi(i))) * eta(i) * (gammap_dot_tr - K_5 * xi_dot_2);
                            dr(1, 1) = p - (a * phi(i) / (phi_m - phi(i)))*(a * phi(i) / (phi_m - phi(i))) * ((gammap_dot_tr - K_5 * xi_dot_2)*(gammap_dot_tr - K_5 * xi_dot_2) * grain_diam * grain_diam * grains_rho + 2 * eta(i) * ((gammap_dot_tr - K_5 * xi_dot_2)));
                            dr(1, 1) = (dr(1, 1) - r(1)) / (tau_bar_plus - tau_bar);
                        } else {
                            break;
                        }

                        upd = dr.inverse()*r;
                        lambda_tmp = 1;
                        //p = p - upd(0);
                        //tau_bar = tau_bar - upd(1);
                        do {
                            p_tmp = p - lambda_tmp*upd(0);
                            tau_bar_tmp = tau_bar - lambda_tmp*upd(1);

                            if (p_tmp < 0) {
                                p_tmp = 0;
                            }

                            if (tau_bar_tmp < 0) {
                                tau_bar_tmp = 0;
                            } else if (tau_bar_tmp > tau_bar_tr) {
                                tau_bar_tmp = tau_bar_tr;
                            }

                            gammap_dot_tr = (tau_bar_tr - tau_bar_tmp) / (G * job->dt);
                            if (p == 0) {
                                mu = mu_2;
                                phi_eq = 0;
                            } else {
                                I_v_tr = eta(i) * gammap_dot_tr / p_tmp;
                                I_tr = gammap_dot_tr * grain_diam * std::sqrt(grains_rho/p_tmp);
                                I_m_tr = std::sqrt(I_tr*I_tr + 2*I_v_tr);
                                mu = mu_1 + (mu_2 - mu_1) / (1 + I_0 / I_m_tr) + 5.0 / 2.0 * phi(i) * I_v_tr / (a*I_m_tr);
                                phi_eq = phi_m / (1 + a*I_m_tr);
                            }
                            if (gammap_dot_tr == 0 || eta(i) == 0) {
                                mu = mu_1;
                            }

                            if (phi(i) > phi_m) {
                                beta = K_3 * (phi(i) - phi_m) + K_4 * (phi(i) - phi_eq);
                            } else {
                                beta = K_4 * (phi(i) - phi_eq);
                            }
                            xi_dot_2 = (p_tmp - p_tr) / (K * job->dt) - beta * gammap_dot_tr;

                            mu_f = mu + beta;
                            if (mu_f < 0){
                                mu_f = 0;
                            }

                            r_tmp(0) = tau_bar_tmp - mu_f * p_tmp;
                            //r_tmp(1) = p_tmp - (phi(i) / (phi_m - phi(i))) * (phi(i) / (phi_m - phi(i))) * eta(i) * (gammap_dot_tr - K_5 * xi_dot_2);
                            r_tmp(1) = p_tmp - (a * phi(i) / (phi_m - phi(i)))*(a * phi(i) / (phi_m - phi(i))) * ((gammap_dot_tr - K_5 * xi_dot_2)*(gammap_dot_tr - K_5 * xi_dot_2) * grain_diam * grain_diam * grains_rho + 2 * eta(i) * ((gammap_dot_tr - K_5 * xi_dot_2)));

                            lambda_tmp *= 0.5;
                        } while (r_tmp.norm() > r.norm() && lambda_tmp > REL_TOL);

                        p = p_tmp;
                        tau_bar = tau_bar_tmp;

                        if (p < 0) {
                            p = 0;
                        }

                        if (tau_bar < 0) {
                            tau_bar = 0;
                        } else if (tau_bar > tau_bar_tr) {
                            tau_bar = tau_bar_tr;
                        }
                    }
                }
            }
        }
        //update stress
        if (p > 0 && tau_bar > 0) {
            T = tau_bar / tau_bar_tr * (T_tr - T_tr.trace() / T_tr.rows() * job->jobTensor<double>(Job::IDENTITY));
            T = T - p * job->jobTensor<double>(Job::IDENTITY);
        } else {
            T = job->jobTensor<double>(Job::ZERO);
        }

        for (size_t pos=0;pos<T.size();pos++){
            /*if (!std::isfinite(T(pos))){
                std::cout << T << std::endl;
                std::cout << tau_bar << std::endl;
                std::cout << p << std::endl;
                exit(0);
            }*/
            body->points.T(i,pos) = T(pos);
        }
        if (SPEC == Material::UPDATE) {
            gammap_dot(i) = gammap_dot_tr;
            gammap(i) += job->dt * gammap_dot(i);
            I_v(i) = eta(i)*gammap_dot(i)/(-T.trace()/T.rows()); //undefined
            I(i) = gammap_dot(i)*grain_diam*std::sqrt(grains_rho/(-T.trace()/T.rows()));
            I_m(i) = (I(i)*I(i) + 2*I_v(i));
        }

    }

    return;
}

/*----------------------------------------------------------------------------*/

void materialAssignStress(Job* job, Body* body, Eigen::MatrixXd stressIN, int idIN, int SPEC){
    for (size_t pos=0;pos<stressIN.size();pos++){
        body->points.T(idIN,pos) = stressIN(pos);
    }
    return;
}

/*----------------------------------------------------------------------------*/

void materialAssignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    Eigen::MatrixXd tmpMat;
    Eigen::VectorXd tmpVec;
    double trT;

    tmpVec = body->points.T.row(idIN).transpose();
    tmpMat = job->jobTensor<double>(tmpVec.data());

    trT = tmpMat.trace();
    tmpMat = tmpMat - (trT/tmpMat.rows() + pressureIN) * job->jobTensor<double>(Job::IDENTITY);
    for (size_t pos=0;pos<tmpMat.size();pos++){
        body->points.T(idIN,pos) = tmpMat(pos);
    }
    return;
}

/*----------------------------------------------------------------------------*/
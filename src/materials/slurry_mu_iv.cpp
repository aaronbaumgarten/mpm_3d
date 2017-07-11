//
// Created by aaron on 7/10/17.
// slurry_mu_iv.cpp
//

//
// Created by aaron on 5/24/17.
// isolin.cpp
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
double TOL = 1e-10;
double h = 1e-7;

//standard linear elastic terms
double E, nu, G, K;
double lambda;

//mu(Iv)
//0.32, 0.7, 0.005, (?), (?), (?), 0.585, (?)
double mu_1, mu_2, I_0, K_3, K_4, K_5, phi_m, grains_rho, eta_0;

//sand
Eigen::VectorXd gammap(0,1);
Eigen::VectorXd gammap_dot(0,1);
Eigen::VectorXd I_v(0,1);
Eigen::VectorXd phi(0,1);


extern "C" void materialWriteFrame(Job* job, Body* body, Serializer* serializer);
extern "C" std::string materialSaveState(Job* job, Body* body, Serializer* serializer, std::string filepath);
extern "C" int materialLoadState(Job* job, Body* body, Serializer* serializer, std::string fullpath);

extern "C" void materialInit(Job* job, Body* body);
extern "C" void materialCalculateStress(Job* job, Body* body, int SPEC);
extern "C" void materialAssignStress(Job* job, Body* body, Eigen::MatrixXd stressIN, int idIN, int SPEC);
extern "C" void materialAssignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC);

/*----------------------------------------------------------------------------*/

void materialInit(Job* job, Body* body){
    if (body->material.fp64_props.size() < 12){
        std::cout << body->material.fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 12 properties defined (E, nu, mu_1, mu_2, I_0, K_3, K_4, K_5, phi_m, grains_rho, eta_0).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        E = body->material.fp64_props[0];
        nu = body->material.fp64_props[1];
        G = E / (2.0 * (1.0 + nu));
        K = E / (3.0 * (1.0 - 2 * nu));
        lambda = K - 2.0 * G / 3.0;
        mu_1 = body->material.fp64_props[3];
        mu_2 = body->material.fp64_props[4];
        I_0 = body->material.fp64_props[5];
        K_3 = body->material.fp64_props[6];
        K_4 = body->material.fp64_props[7];
        K_5 = body->material.fp64_props[8];
        phi_m = body->material.fp64_props[9];
        grains_rho = body->material.fp64_props[10];
        eta_0 = body->material.fp64_props[11];

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, mu_1 = %g, mu_2 = %g, I_0 = %g, K_3 = %g, K_4 = %g, K_5 = %g, phi_m = %g, grains_rho = %g, eta_0 = %g).\n",
               E, nu, G, K, mu_1, mu_2, I_0, K_3, K_4, K_5, phi_m, grains_rho, eta_0);
    }

    gammap_dot.resize(body->points.x.rows());
    gammap_dot.setZero();
    gammap.resize(body->points.x.rows());
    gammap.setZero();
    I_v.resize(body->points.x.rows());
    I_v.setZero();
    phi.resize(body->points.x.rows());
    phi.setZero();

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void materialWriteFrame(Job* job, Body* body, Serializer* serializer) {
    //nothing to report
    serializer->serializerWriteScalarArray(gammap,"gammap");
    serializer->serializerWriteScalarArray(gammap_dot,"gammap_dot");
    serializer->serializerWriteScalarArray(I_v,"I_v");
    serializer->serializerWriteScalarArray(phi,"packing_fraction");
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
        ffile << I_0 << "\n";
        ffile << K_3 << "\n";
        ffile << K_4 << "\n";
        ffile << K_5 << "\n";
        ffile << phi_m << "\n";
        ffile << grains_rho << "\n";
        ffile << eta_0 << "\n";

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
        std::getline(fin,line); //I_0
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

        printf("Material properties (E = %g, nu = %g, G = %g, K = %g, mu_1 = %g, mu_2 = %g, I_0 = %g, K_3 = %g, K_4 = %g, K_5 = %g, phi_m = %g, grains_rho = %g, eta_0 = %g).\n",
               E, nu, G, K, mu_1, mu_2, I_0, K_3, K_4, K_5, phi_m, grains_rho, eta_0);
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

    double trD, p, p_tr, tau_bar, tau_bar_tr, I_v_tr;
    double xi_dot_1, xi_dot_2, gammap_dot_tr, mu, beta, phi_eq;

    Eigen::MatrixXd tmpMat = job->jobTensor<double>();
    Eigen::VectorXd tmpVec;

    //residuals
    Eigen::Vector2d r, b, upd;
    Eigen::Matrix2d dr;
    double gdp, pt;
    double r_k5, b_k5, dr_k5;

    //calculate packing fraction
    phi = 1/grains_rho * (body->points.m.array() / body->points.v.array());

    for (size_t i=0;i<body->points.x.rows();i++){
        if (body->points.active[i] == 0){
            continue;
        }

        tmpVec = body->points.L.row(i).transpose();
        L = job->jobTensor<double>(tmpVec.data());

        tmpVec = body->points.T.row(i).transpose();
        T = job->jobTensor<double>(tmpVec.data());

        D = 0.5*(L+L.transpose());
        W = 0.5*(L-L.transpose());

        trD = D.trace();
        D_0 = D - trD/3.0 * job->jobTensor<double>(Job::IDENTITY);

        //trial stress
        T_tr = T + job->dt * ((2*G*D) + (lambda*trD*job->jobTensor<double>(Job::IDENTITY)) + (W*T) - (T*W));
        tau_bar_tr = (T_tr - (T_tr.trace() / 3.0)*job->jobTensor<double>(Job::IDENTITY)).norm() / std::sqrt(2.0);
        p_tr = -T_tr.trace() / 3.0;

        //check if trial stress admissible
        if (tau_bar_tr <= mu_1*p_tr && p_tr > 0){
            if (phi(i) >= phi_m){
                for (size_t pos=0;pos<body->points.T.cols();pos++){
                    body->points.T(i,pos) = T_tr(pos);
                    if (SPEC == Material::UPDATE){
                        gammap_dot(i) = 0;
                        I_v(i) = 0;
                    }
                }
                continue;
            } else {
                //check if stress state lies on f3 only
                xi_dot_2 = -p_tr/((phi(i)/(phi_m - phi(i)))*(phi(i)/(phi_m - phi(i)))*eta_0*K_5 + K*job->dt);
                p = p_tr + K*job->dt*xi_dot_2;
                if (tau_bar_tr <= mu_1*p && p >= 0){
                    for (size_t pos=0;pos<body->points.T.cols();pos++) {
                        body->points.T(i, pos) = T_tr(pos);
                        if (SPEC == Material::UPDATE) {
                            gammap_dot(i) = 0;
                            I_v(i) = 0;
                        }
                    }
                    continue;
                }
            }
        }


        //check if stress state lies on f1,f2 yield surfaces
        if (p_tr <= 0){
            beta = K_3*(phi(i) - phi_m) + K_4*(phi(i) - 0.0);
            gammap_dot_tr = sqrt(2)*D_0.norm();
            if ((p_tr + K*job->dt*(beta*gammap_dot_tr)) <= 0){
                body->points.T.row(i).setZero();
                if (SPEC == Material::UPDATE) {
                    gammap_dot(i) = gammap_dot_tr;
                    gammap(i) += job->dt * gammap_dot(i);
                    I_v(i) = -1; //undefined
                }
                continue;
            }
        }

        //assume stress state lies on f1 yield surface only
        r << tau_bar_tr, p_tr;
        b = r;
        dr.setZero();
        p = p_tr;
        gammap_dot_tr = gammap_dot(i);

        while (r.norm() > b.norm()*TOL){
            if (p > 0){
                I_v_tr = eta_0 * gammap_dot_tr / p;
                mu = mu_1 + (mu_2 - mu_1) / (1 + I_0/I_v_tr) + 5.0/2.0 * phi(i) * sqrt(I_v_tr);
                phi_eq = phi_m / 1+sqrt(I_v_tr);
            } else {
                mu = mu_1;
                phi_eq = 0;
            }
            tau_bar = tau_bar_tr - G*gammap_dot_tr*job->dt;
            beta = K_4 * (phi(i) - phi_eq);
            if (phi(i) > phi_m){
                beta += K_3 * (phi(i) - phi_m);
            }

            r(0) = tau_bar - (mu + beta)*p;
            r(1) = p - p_tr - K*job->dt*(beta*gammap_dot_tr);
            if (r.norm() > b.norm()*TOL){
                //dr/dgamma
                if (gammap_dot_tr > 0) {
                    gdp = gammap_dot_tr * (1 + h);
                } else {
                    gdp = h;
                }
                if (p > 0){
                    I_v_tr = eta_0 * gdp / p;
                    mu = mu_1 + (mu_2 - mu_1) / (1 + I_0/I_v_tr) + 5.0/2.0 * phi(i) * sqrt(I_v_tr);
                    phi_eq = phi_m / 1+sqrt(I_v_tr);
                } else {
                    I_v_tr = -1; //undef
                    mu = mu_1;
                    phi_eq = 0;
                }
                tau_bar = tau_bar_tr - G*gdp*job->dt;
                beta = K_4 * (phi(i) - phi_eq);
                if (phi(i) > phi_m){
                    beta += K_3 * (phi(i) - phi_m);
                }

                dr(0,0) = (tau_bar - (mu + beta)*p - r(0))/(gdp - gammap_dot_tr);
                dr(1,0) = (p - p_tr - K*job->dt*(beta*gdp) -r(1))/(gdp - gammap_dot_tr);

                //dr/dp
                if (p > 0){
                    pt = p*(1+h);
                } else {
                    pt = h;
                }
                if (pt > 0){
                    I_v_tr = eta_0 * gammap_dot_tr / pt;
                    mu = mu_1 + (mu_2 - mu_1) / (1 + I_0/I_v_tr) + 5.0/2.0 * phi(i) * sqrt(I_v_tr);
                    phi_eq = phi_m / 1+sqrt(I_v_tr);
                } else {
                    I_v_tr = -1; //undef
                    mu = mu_1;
                    phi_eq = 0;
                }
                tau_bar = tau_bar_tr - G*gammap_dot_tr*job->dt;
                beta = K_4 * (phi(i) - phi_eq);
                if (phi(i) > phi_m){
                    beta += K_3 * (phi(i) - phi_m);
                }

                dr(0,0) = (tau_bar - (mu + beta)*pt - r(0))/(pt - p);
                dr(1,0) = (pt - p_tr - K*job->dt*(beta*gammap_dot_tr) -r(1))/(pt - p);
            } else {
                T = tau_bar / tau_bar_tr * (T_tr - T_tr.trace()/3.0 * job->jobTensor<double>(Job::IDENTITY)) - p*job->jobTensor<double>(Job::IDENTITY);
                break;
            }

            upd = dr.colPivHouseholderQr().solve(r);
            gammap_dot_tr = gammap_dot_tr - upd(0);
            p = p - upd(1);
        }

        //check if solution lies on f3 yield surface
        if (phi(i) >= phi_m){
            for (size_t pos=0;pos<T.size();pos++){
                body->points.T(i,pos) = T(pos);
            }
            if (SPEC == Material::UPDATE) {
                gammap_dot(i) = gammap_dot_tr;
                gammap(i) += job->dt * gammap_dot(i);
                I_v(i) = I_v_tr; //undefined
            }
            continue;
        } else if (p <= (phi(i) / (phi_m - phi(i))) * (phi(i) / (phi_m - phi(i))) * eta_0 * gammap_dot_tr ){
            for (size_t pos=0;pos<T.size();pos++){
                body->points.T(i,pos) = T(pos);
            }
            if (SPEC == Material::UPDATE) {
                gammap_dot(i) = gammap_dot_tr;
                gammap(i) += job->dt * gammap_dot(i);
                I_v(i) = I_v_tr; //undefined
            }
            continue;
        } else if (K_5 == 0){
            //cannot support pressure
            r_k5 = tau_bar_tr;
            b_k5 = r_k5;
            dr_k5 = 0;
            gammap_dot_tr = gammap_dot(i);

            while (std::abs(r_k5) > std::abs(b_k5*TOL)){
                p = (phi(i)*phi(i)) / ((phi_m - phi(i))*(phi_m - phi(i))) * eta_0 * gammap_dot_tr;
                if (p > 0){
                    I_v_tr = eta_0 * gammap_dot_tr / p;
                    mu = mu_1 + (mu_2 - mu_1) / (1 + I_0/I_v_tr) + 5.0/2.0 * phi(i) * sqrt(I_v_tr);
                    phi_eq = phi_m / 1+sqrt(I_v_tr);
                } else {
                    mu = mu_1;
                    phi_eq = 0;
                }
                tau_bar = tau_bar_tr - G*gammap_dot_tr*job->dt;
                beta = K_4 * (phi(i) - phi_eq);
                if (phi(i) > phi_m){
                    beta += K_3 * (phi(i) - phi_m);
                }

                r_k5 = tau_bar - (mu + beta)*p;

                if (std::abs(r_k5) > std::abs(b_k5*TOL)){
                    if (gammap_dot_tr > 0) {
                        gdp = gammap_dot_tr * (1 + h);
                    } else {
                        gdp = h;
                    }

                    if (p > 0){
                        I_v_tr = eta_0 * gdp / p;
                        mu = mu_1 + (mu_2 - mu_1) / (1 + I_0/I_v_tr) + 5.0/2.0 * phi(i) * sqrt(I_v_tr);
                        phi_eq = phi_m / 1+sqrt(I_v_tr);
                    } else {
                        mu = mu_1;
                        phi_eq = 0;
                    }
                    tau_bar = tau_bar_tr - G*gdp*job->dt;
                    beta = K_4 * (phi(i) - phi_eq);
                    if (phi(i) > phi_m){
                        beta += K_3 * (phi(i) - phi_m);
                    }

                    p = (phi(i)*phi(i)) / ((phi_m - phi(i))*(phi_m - phi(i))) * eta_0 * gdp;

                    dr_k5 = (tau_bar - (mu + beta)*p - r_k5)/(gdp - gammap_dot_tr);
                } else {
                    T = tau_bar / tau_bar_tr * (T_tr - T_tr.trace()/3.0 * job->jobTensor<double>(Job::IDENTITY)) - p*job->jobTensor<double>(Job::IDENTITY);
                    break;
                }

                gammap_dot_tr = gammap_dot_tr - r_k5/dr_k5;
            }

        } else {
            r << tau_bar_tr, p_tr;
            b = r;
            dr.setZero();
            p = p_tr;
            gammap_dot_tr = gammap_dot(i);

            while (r.norm() > b.norm()*TOL){
                if (p > 0){
                    I_v_tr = eta_0 * gammap_dot_tr / p;
                    mu = mu_1 + (mu_2 - mu_1) / (1 + I_0/I_v_tr) + 5.0/2.0 * phi(i) * sqrt(I_v_tr);
                    phi_eq = phi_m / 1+sqrt(I_v_tr);
                } else {
                    mu = mu_1;
                    phi_eq = 0;
                }
                tau_bar = tau_bar_tr - G*gammap_dot_tr*job->dt;
                beta = K_4 * (phi(i) - phi_eq);
                if (phi(i) > phi_m){
                    beta += K_3 * (phi(i) - phi_m);
                }
                xi_dot_2 = (gammap_dot_tr - p *((phi_m - phi(i))*(phi_m - phi(i))) / (phi(i)*phi(i)) / eta_0) / K_5;

                r(0) = tau_bar - (mu + beta)*p;
                r(1) = p - p_tr - K*job->dt*(beta*gammap_dot_tr - K_5*xi_dot_2);
                if (r.norm() > b.norm()*TOL){
                    //dr/dgamma
                    if (gammap_dot_tr > 0) {
                        gdp = gammap_dot_tr * (1 + h);
                    } else {
                        gdp = h;
                    }

                    if (p > 0){
                        I_v_tr = eta_0 * gdp / p;
                        mu = mu_1 + (mu_2 - mu_1) / (1 + I_0/I_v_tr) + 5.0/2.0 * phi(i) * sqrt(I_v_tr);
                        phi_eq = phi_m / 1+sqrt(I_v_tr);
                    } else {
                        I_v_tr = -1; //undef
                        mu = mu_1;
                        phi_eq = 0;
                    }
                    tau_bar = tau_bar_tr - G*gdp*job->dt;
                    beta = K_4 * (phi(i) - phi_eq);
                    if (phi(i) > phi_m){
                        beta += K_3 * (phi(i) - phi_m);
                    }
                    xi_dot_2 = (gdp - p *((phi_m - phi(i))*(phi_m - phi(i))) / (phi(i)*phi(i)) / eta_0) / K_5;

                    dr(0,0) = (tau_bar - (mu + beta)*p - r(0))/(gdp - gammap_dot_tr);
                    dr(1,0) = (p - p_tr - K*job->dt*(beta*gdp + xi_dot_2) -r(1))/(gdp - gammap_dot_tr);

                    //dr/dp
                    if (p > 0){
                        pt = p*(1+h);
                    } else {
                        pt = h;
                    }
                    if (pt > 0){
                        I_v_tr = eta_0 * gammap_dot_tr / pt;
                        mu = mu_1 + (mu_2 - mu_1) / (1 + I_0/I_v_tr) + 5.0/2.0 * phi(i) * sqrt(I_v_tr);
                        phi_eq = phi_m / 1+sqrt(I_v_tr);
                    } else {
                        I_v_tr = -1; //undef
                        mu = mu_1;
                        phi_eq = 0;
                    }
                    tau_bar = tau_bar_tr - G*gammap_dot_tr*job->dt;
                    beta = K_4 * (phi(i) - phi_eq);
                    if (phi(i) > phi_m){
                        beta += K_3 * (phi(i) - phi_m);
                    }
                    xi_dot_2 = (gammap_dot_tr - pt *((phi_m - phi(i))*(phi_m - phi(i))) / (phi(i)*phi(i)) / eta_0) / K_5;

                    dr(0,0) = (tau_bar - (mu + beta)*pt - r(0))/(pt - p);
                    dr(1,0) = (pt - p_tr - K*job->dt*(beta*gammap_dot_tr + xi_dot_2) -r(1))/(pt - p);
                } else {
                    T = tau_bar / tau_bar_tr * (T_tr - T_tr.trace()/3.0 * job->jobTensor<double>(Job::IDENTITY)) - p*job->jobTensor<double>(Job::IDENTITY);
                    break;
                }

                upd = dr.colPivHouseholderQr().solve(r);
                gammap_dot_tr = gammap_dot_tr - upd(0);
                p = p - upd(1);
            }
        }

        for (size_t pos=0;pos<T.size();pos++){
            body->points.T(i,pos) = T(pos);
        }
        if (SPEC == Material::UPDATE) {
            gammap_dot(i) = gammap_dot_tr;
            gammap(i) += job->dt * gammap_dot(i);
            I_v(i) = I_v_tr; //undefined
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
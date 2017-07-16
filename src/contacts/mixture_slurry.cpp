//
// Created by aaron on 6/10/17.
// mixture_slurry.cpp
//

//
// Created by aaron on 5/26/17.
// contact_huang.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "stringparser.hpp"

#include "job.hpp"

#include "serializer.hpp"
#include "contact.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "points.hpp"

int id = 0;
double grains_rho = 2500; // kg/m^3
double grains_diam = 0.001; // m
double k = 0.0055; //permeability constant
double mu_w = 8.9e-4; //viscosity of fluid

std::vector<int> bodyIDs = {-1,-1};

int solid_body_id = -1;
int liquid_body_id = -1;

Eigen::VectorXd n(0,1);
Eigen::VectorXd V(0,1);

Eigen::MatrixXd divT(0,0);

extern "C" void contactWriteFrame(Job* job, Serializer* serializer); //write frame to serializer
extern "C" std::string contactSaveState(Job* job, Serializer* serializer, std::string filepath); //save state to serializer folder with returned filename
extern "C" int contactLoadState(Job* job, Serializer* serializer, std::string fullpath); //read state from given full path

extern "C" void contactInit(Job* job, Contact* contact); //initialize contact
extern "C" void contactGenerateRules(Job* job); //generate contact rules
extern "C" void contactApplyRules(Job* job, int SPEC); //apply rules

/*----------------------------------------------------------------------------*/

void contactInit(Job* job, Contact* contact){
    //set id and initialize body id vector
    id = contact->id;
    bodyIDs = {-1,-1};

    //check that contact properties are set
    if (contact->fp64_props.size() < 4 || (contact->str_props.size() < 2 && contact->int_props.size() < 2)){
        //need to coefficient of friction and bodies
        std::cout << contact->fp64_props.size() << ", " << contact->int_props.size() << ", " << contact->str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 6 properties defined ({grain_rho, grain_diam, k, mu_w},{body_1,body_2}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //set coeff of friction
        grains_rho = contact->fp64_props[0];
        grains_diam = contact->fp64_props[1];
        k = contact->fp64_props[2];
        mu_w = contact->fp64_props[3];

        //set body ids by name
        if (contact->str_props.size() == 2){
            for (size_t i=0;i<bodyIDs.size();i++) {
                for (size_t b = 0; b < job->bodies.size(); b++) {
                    if (contact->str_props[i].compare(job->bodies[b].name) == 0){
                        bodyIDs[i] = b;
                        break;
                    }
                }
            }
        }

        // or set body ids by int
        for (size_t i=0;i<bodyIDs.size();i++) {
            if (bodyIDs[i] < 0){
                if (contact->int_props.size() == 2) {
                    bodyIDs = contact->int_props;
                } else {
                    std::cout << contact->fp64_props.size() << ", " << contact->int_props.size() << ", " << contact->str_props.size() << "\n";
                    fprintf(stderr,
                            "%s:%s: Need at least 6 properties defined ({grain_rho, grain_diam, k, mu_w},{solid,fluid}).\n",
                            __FILE__, __func__);
                    exit(0);
                }
                break;
            }
        }

        solid_body_id = bodyIDs[0];
        liquid_body_id = bodyIDs[1];

        V.resize(job->bodies[liquid_body_id].nodes.x.rows());
        n.resize(job->bodies[solid_body_id].nodes.x.rows());

        divT = job->jobVectorArray<double>(job->bodies[solid_body_id].nodes.x.rows());

        printf("Contact properties (grain_rho = %g, grain_diam = %g, k = %g, mu_w = %g, {solid: %i, fluid: %i}).\n",
               grains_rho, grains_diam, k, mu_w, solid_body_id, liquid_body_id);
    }

    std::cout << "Contact Initialized: [" << id << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void contactWriteFrame(Job* job, Serializer* serializer){
    serializer->serializerWriteScalarArray(n,("n_" + job->bodies[solid_body_id].name + "_" + job->bodies[liquid_body_id].name));
    return;
}

/*----------------------------------------------------------------------------*/

std::string contactSaveState(Job* job, Serializer* serializer, std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename
    std::stringstream s;
    s << "mpm_v2.contact." << id << "." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    filename = s.str();
    std::ofstream ffile((filepath+filename), std::ios::trunc);

    if (ffile.is_open()){
        ffile << "# mpm_v2 grids/cartesian.so\n";
        ffile << bodyIDs[0] << " " << bodyIDs[1] << "\n"; //body ids
        ffile << id << "\n"; //contact id
        ffile << grains_rho << "\n";
        ffile << grains_diam << "\n";
        ffile << k << "\n";
        ffile << mu_w << "\n";
        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filepath+filename << "\" !\n";
        return "ERR";
    }

    std::cout << "Contact Saved: [" << id << "]." << std::endl;

    return filename;
}

/*----------------------------------------------------------------------------*/

int contactLoadState(Job* job, Serializer* serializer, std::string fullpath){
    bodyIDs = {-1,-1};
    std::string line;
    std::stringstream ss;
    std::ifstream fin(fullpath);

    if(fin.is_open()){
        std::getline(fin,line); //first line

        std::getline(fin,line); //body IDs
        ss = std::stringstream(line);
        if (!(ss >> bodyIDs[0] >> bodyIDs[1])){
            std::cout << "Conact Loading Error! Unable to read body IDs in file: " << fullpath << std::endl;
            return 0;
        }

        std::getline(fin,line); //contact id
        id = std::stoi(line);

        std::getline(fin,line); //grains_rho
        grains_rho = std::stod(line);

        std::getline(fin,line); //grains_diam
        grains_diam = std::stod(line);

        std::getline(fin,line); //k
        k = std::stod(line);

        std::getline(fin,line); //mu_w
        mu_w = std::stod(line);

        fin.close();

        solid_body_id = bodyIDs[0];
        liquid_body_id = bodyIDs[1];

        V.resize(job->bodies[liquid_body_id].nodes.x.rows());
        n.resize(job->bodies[solid_body_id].nodes.x.rows());

        divT = job->jobVectorArray<double>(job->bodies[solid_body_id].nodes.x.rows());

        printf("Contact properties (grain_rho = %g, grain_diam = %g, k = %g, mu_w = %g, {solid: %i, fluid: %i}).\n",
               grains_rho, grains_diam, k, mu_w, solid_body_id, solid_body_id);
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Contact Loaded." << std::endl;
    return 1;
}

/*----------------------------------------------------------------------------*/

void contactGenerateRules(Job* job){

    //set porosity from solid mass
    for (size_t i=0;i<n.rows();i++){
        //n = 1 - phi
        n(i) = 1 - (job->bodies[solid_body_id].nodes.m(i) / (job->grid.gridNodalVolume(job,i)*grains_rho));
    }

    /*Eigen::VectorXd solid_n = job->bodies[solid_body_id].points.m.array() * job->bodies[solid_body_id].points.m.array() / job->bodies[solid_body_id].points.v.array();
    job->bodies[solid_body_id].bodyCalcNodalValues(job,n,solid_n,Body::SET);
    for (size_t i = 0; i < n.rows(); i++) {
        if (job->bodies[solid_body_id].nodes.m(i) > 0){
            n(i) = 1 - n(i)/(job->bodies[solid_body_id].nodes.m(i)*grains_rho);
        } else {
            n(i) = 1;
        }
    }*/

    //approximate V as integrated liquid volume
    job->bodies[liquid_body_id].bodyCalcNodalValues<Eigen::VectorXd,Eigen::VectorXd>(job,V,job->bodies[liquid_body_id].points.v,Body::SET);

    return;
}

/*----------------------------------------------------------------------------*/

void contactApplyRules(Job* job, int SPEC){
    double m1, m2;
    Eigen::VectorXd fsfi = job->jobVector<double>();
    Eigen::VectorXd mv1i = job->jobVector<double>();
    Eigen::VectorXd mv2i = job->jobVector<double>();
    Eigen::VectorXd vCMi = job->jobVector<double>();

    double C, k_eff;

    //calculate divergence of liquid stress
    job->bodies[liquid_body_id].bodyCalcNodalDivergence(job,divT,job->bodies[liquid_body_id].points.T,Body::SET);

    //calculate nodal momentum exchange
    for (size_t i = 0; i < n.rows(); i++) {
        //test every node for contact
        if (job->bodies[solid_body_id].nodes.m[i] > 0 && job->bodies[liquid_body_id].nodes.m[i] > 0) {
            //determine 'center of mass' velocity
            m1 = job->bodies[solid_body_id].nodes.m[i];
            m2 = job->bodies[liquid_body_id].nodes.m[i];
            if (SPEC == Contact::EXPLICIT) {
                mv1i << (job->bodies[solid_body_id].nodes.mx_t.row(i) + job->dt * job->bodies[solid_body_id].nodes.f.row(i)).transpose();
                mv2i << (job->bodies[liquid_body_id].nodes.mx_t.row(i) + job->dt * job->bodies[liquid_body_id].nodes.f.row(i)).transpose();
            } else if (SPEC == Contact::IMPLICIT){
                mv1i << job->bodies[solid_body_id].nodes.mx_t.row(i).transpose();
                mv2i << job->bodies[liquid_body_id].nodes.mx_t.row(i).transpose();
            } else {
                std::cerr << "ERROR: Unknown SPEC in mixture_slurry.so: " << SPEC << "!" << std::endl;
                return;
            }
            vCMi  = (mv1i + mv2i) / (m1 + m2);

            //permeability
            //k = 0.0055 * n^3 * d^2 / (1-n)^2
            k_eff = k*n(i)*n(i)*n(i)*grains_diam*grains_diam / ((1-n(i))*(1-n(i)));

            //f_d = (v_s - v_f)*C
            //C = n^2 * mu_w/k_eff
            //(v_s - f_d/m_s*dt - v_f - f_d/m_f*dt)/(v_s - v_f) > 0
            //1 - C*dt*(1/m_s + 1/m_f) > 0
            //C*dt*(1/m_s + 1/m_f) < 1
            C = n(i) * V(i) * mu_w / k_eff;
            if (!std::isfinite(k_eff)){
                C = 0;
            } else if ((C*job->dt*(1.0/m1 + 1.0/m2)) > 1){
                C = 1.0/(job->dt*(1.0/m1 + 1.0/m2)); //enforce damping
            }

            double one_minus_n = 1 - job->bodies[liquid_body_id].nodes.m(i) / (job->grid.gridNodalVolume(job,i) * 1000);

            fsfi = (mv1i / m1 - mv2i / m2)*C - one_minus_n*divT.row(i).transpose();//(1-n(i))*divT.row(i).transpose();

            job->bodies[solid_body_id].nodes.f.row(i) -= fsfi.transpose();
            job->bodies[liquid_body_id].nodes.f.row(i) += fsfi.transpose();
        }
    }

    return;
}


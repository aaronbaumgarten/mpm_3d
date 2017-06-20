//
// Created by aaron on 6/20/17.
// contact_rigid.cpp
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
double mu_f = 0.4;
std::vector<int> bodyIDs = {-1,-1};
Eigen::MatrixXd contact_normal(0,0);
Eigen::MatrixXd rigid_force(0,0);
Eigen::MatrixXd f_sum_write(0,0);
Eigen::MatrixXd mv1_k;
Eigen::MatrixXd mv2_k;
Eigen::MatrixXd f1_k;
Eigen::MatrixXd f2_k;

extern "C" void contactWriteFrame(Job* job, Serializer* serializer); //write frame to serializer
extern "C" std::string contactSaveState(Job* job, Serializer* serializer, std::string filepath); //save state to serializer folder with returned filename
extern "C" int contactLoadState(Job* job, Serializer* serializer, std::string fullpath); //read state from given full path

extern "C" void contactInit(Job* job, Contact* contact); //initialize contact
extern "C" void contactGenerateRules(Job* job); //generate contact rules
extern "C" void contactApplyRules(Job* job, int); //apply rules

/*----------------------------------------------------------------------------*/

void contactWriteFrame(Job* job, Serializer* serializer){
    for (size_t pos=0; pos<rigid_force.cols(); pos++){
        f_sum_write.col(0) = Eigen::VectorXd::Ones(f_sum_write.rows()) * rigid_force.col(pos).sum();
    }
    serializer->serializerWriteVectorArray(contact_normal,("rigid_normal_" + job->bodies[bodyIDs[0]].name + "_" + job->bodies[bodyIDs[1]].name));
    serializer->serializerWriteVectorArray(rigid_force,("rigid_force_on_" + job->bodies[bodyIDs[0]].name));
    serializer->serializerWriteVectorArray(f_sum_write,("rigid_total_force_" + job->bodies[bodyIDs[0]].name));
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
        ffile << id << "\n" << mu_f << "\n"; //contact if and mu_f
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

        std::getline(fin,line); //mu_f
        mu_f = std::stod(line);


        printf("Contact properties (mu_f = %g).\n",
               mu_f);

        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    contact_normal = job->jobVectorArray<double>(job->bodies[bodyIDs[0]].nodes.x.rows());
    rigid_force = job->jobVectorArray<double>(job->bodies[bodyIDs[0]].nodes.x.rows());
    f_sum_write = job->jobVectorArray<double>(job->bodies[bodyIDs[0]].points.x.rows());

    std::cout << "Contact Loaded." << std::endl;
    return 1;
}

/*----------------------------------------------------------------------------*/

void contactInit(Job* job, Contact* contact){
    //set id and initialize body id vector
    id = contact->id;
    bodyIDs = {-1,-1};

    //check that contact properties are set
    if (contact->fp64_props.size() < 1 || (contact->str_props.size() < 2 && contact->int_props.size() < 2)){
        //need to coefficient of friction and bodies
        std::cout << contact->fp64_props.size() << ", " << contact->int_props.size() << ", " << contact->str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 3 properties defined ({mu_f},{body_1,body_2}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //set coeff of friction
        mu_f = contact->fp64_props[0];

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
                            "%s:%s: Need at least 3 properties defined ({mu_f},{body_1,body_2}).\n",
                            __FILE__, __func__);
                    exit(0);
                }
                break;
            }
        }

        contact_normal = job->jobVectorArray<double>(job->bodies[bodyIDs[0]].nodes.x.rows());
        rigid_force = job->jobVectorArray<double>(job->bodies[bodyIDs[0]].nodes.x.rows());
        f_sum_write = job->jobVectorArray<double>(job->bodies[bodyIDs[0]].points.x.rows());

        printf("Contact properties (mu_f = %g, {%i, %i}).\n",
               mu_f, bodyIDs[0], bodyIDs[1]);
    }

    std::cout << "Contact Initialized: [" << id << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void contactGenerateRules(Job* job){
    //set normal for problem
    //use normal from body 1
    job->bodies[bodyIDs[0]].bodyCalcNodalGradient<Eigen::MatrixXd,Eigen::VectorXd>(job,contact_normal,job->bodies[bodyIDs[0]].points.m,Body::SET);
    for (size_t i=0;i<contact_normal.rows();i++){
        //normalize
        contact_normal.row(i) *= -1.0/contact_normal.row(i).norm();
    }

    rigid_force.setZero();

    //store initial velocity for implicit update
    mv1_k = job->bodies[bodyIDs[0]].nodes.mx_t;
    mv2_k = job->bodies[bodyIDs[1]].nodes.mx_t;

    f1_k = job->bodies[bodyIDs[0]].nodes.f;
    f2_k = job->bodies[bodyIDs[1]].nodes.f;

    return;
}

/*----------------------------------------------------------------------------*/

void contactApplyRules(Job* job, int SPEC){
    Eigen::VectorXd normal = job->jobVector<double>();
    double m1, m2;
    double fn1i, ft1i;
    Eigen::VectorXd mv1i = job->jobVector<double>();
    Eigen::VectorXd mv2i = job->jobVector<double>();
    Eigen::VectorXd vCMi = job->jobVector<double>();
    Eigen::VectorXd fcti = job->jobVector<double>();
    Eigen::VectorXd s1i = job->jobVector<double>();

    Eigen::VectorXd tmpVec = job->jobVector<double>();

    size_t b1 = bodyIDs[0];
    size_t b2 = bodyIDs[1];

    //look for contacts if there are two bodies
    for (size_t i = 0; i < contact_normal.rows(); i++) {
        //test every node for contact
        if (job->bodies[b1].nodes.m[i] > 0 && job->bodies[b2].nodes.m[i] > 0) {

            normal << -contact_normal.row(i).transpose();

            //determine 'center of mass' velocity
            m1 = job->bodies[b1].nodes.m[i];
            m2 = job->bodies[b2].nodes.m[i];

            if (SPEC == Contact::EXPLICIT) {
                mv1i << (job->bodies[b1].nodes.mx_t.row(i)).transpose();
                mv2i << (job->bodies[b2].nodes.mx_t.row(i) + job->dt * job->bodies[b2].nodes.f.row(i)).transpose();
            } else if (SPEC == Contact::IMPLICIT){
                //mv1i << (mv1_k.row(i) + job->dt * job->bodies[b1].nodes.f.row(i)).transpose();
                mv2i << (mv2_k.row(i) + job->dt * job->bodies[b2].nodes.f.row(i)).transpose();
                mv1i << job->bodies[b1].nodes.mx_t.row(i).transpose();
                //mv2i << job->bodies[b2].nodes.mx_t.row(i).transpose();
            } else {
                std::cerr << "ERROR: Unknown SPEC in contact_huang.so: " << SPEC << "!" << std::endl;
                return;
            }

            //check if converging
            if ((mv2i / m2 - mv1i / m1).dot(normal) > 0) {
                //determine normal force
                fn1i = m2 / job->dt * (mv1i / m1 - mv2i / m2).dot(normal);

                //determine shear force and shear vector
                s1i = m2 / job->dt * (mv1i / m1 - mv2i / m2) - fn1i * normal;
                ft1i = sqrt(s1i.dot(s1i));
                s1i /= ft1i;

                //add forces
                fcti = std::min(0.0, fn1i) * normal + std::min(mu_f * std::abs(fn1i), std::abs(ft1i)) * s1i;

                //set contact forces
                rigid_force.row(i) = -fcti.transpose();
                job->bodies[b2].nodes.f.row(i) += fcti.transpose();
            }
        }
    }

    return;
}


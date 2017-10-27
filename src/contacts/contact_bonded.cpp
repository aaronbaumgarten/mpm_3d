//
// Created by aaron on 7/19/17.
// contact_traction.cpp
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
std::vector<int> bodyIDs = {-1,-1};
Eigen::MatrixXd interaction_force(0,0);
Eigen::MatrixXd mv1_k;
Eigen::MatrixXd mv2_k;

extern "C" void contactWriteFrame(Job* job, Serializer* serializer); //write frame to serializer
extern "C" std::string contactSaveState(Job* job, Serializer* serializer, std::string filepath); //save state to serializer folder with returned filename
extern "C" int contactLoadState(Job* job, Serializer* serializer, std::string fullpath); //read state from given full path

extern "C" void contactInit(Job* job, Contact* contact); //initialize contact
extern "C" void contactGenerateRules(Job* job); //generate contact rules
extern "C" void contactApplyRules(Job* job, int); //apply rules

/*----------------------------------------------------------------------------*/

void contactWriteFrame(Job* job, Serializer* serializer){
    serializer->serializerWriteVectorArray(interaction_force,("interaction_force_" + job->bodies[bodyIDs[0]].name + "_" + job->bodies[bodyIDs[1]].name));
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
        ffile << "# mpm_v2 contacts/contact_trction.so\n";
        ffile << bodyIDs[0] << " " << bodyIDs[1] << "\n"; //body ids
        ffile << id << "\n"; //contact id
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

        printf("Contact properties ({%i, %i}).\n", bodyIDs[0], bodyIDs[1]);

        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    interaction_force = job->jobVectorArray<double>(job->bodies[bodyIDs[0]].nodes.x.rows());

    std::cout << "Contact Loaded." << std::endl;
    return 1;
}

/*----------------------------------------------------------------------------*/

void contactInit(Job* job, Contact* contact){
    //set id and initialize body id vector
    id = contact->id;
    bodyIDs = {-1,-1};

    //check that contact properties are set
    if ((contact->str_props.size() < 2 && contact->int_props.size() < 2)){
        //need to coefficient of friction and bodies
        std::cout << contact->fp64_props.size() << ", " << contact->int_props.size() << ", " << contact->str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined ({body_1,body_2}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
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
                            "%s:%s: Need at least 2 properties defined ({body_1,body_2}).\n",
                            __FILE__, __func__);
                    exit(0);
                }
                break;
            }
        }

        interaction_force = job->jobVectorArray<double>(job->bodies[bodyIDs[0]].nodes.x.rows());

        printf("Contact properties ({%i, %i}).\n", bodyIDs[0], bodyIDs[1]);
    }

    std::cout << "Contact Initialized: [" << id << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/

void contactGenerateRules(Job* job){
    //set normal for problem
    //use normal from body 1
    //store initial velocity for implicit update
    mv1_k = job->bodies[bodyIDs[0]].nodes.mx_t;
    mv2_k = job->bodies[bodyIDs[1]].nodes.mx_t;

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

    Eigen::VectorXd tmpVec = job->jobVector<double>();

    size_t b1 = bodyIDs[0];
    size_t b2 = bodyIDs[1];

    //look for contacts if there are two bodies
    for (size_t i = 0; i < interaction_force.rows(); i++) {
        //test every node for contact
        if (job->bodies[b1].nodes.m[i] > 0 && job->bodies[b2].nodes.m[i] > 0) {
            //determine 'center of mass' velocity
            m1 = job->bodies[b1].nodes.m[i];
            m2 = job->bodies[b2].nodes.m[i];

            mv1i << (mv1_k.row(i) + job->dt * job->bodies[b1].nodes.f.row(i)).transpose();
            mv2i << (mv2_k.row(i) + job->dt * job->bodies[b2].nodes.f.row(i)).transpose();
            //mv1i << job->bodies[b1].nodes.mx_t.row(i).transpose();
            //mv2i << job->bodies[b2].nodes.mx_t.row(i).transpose();

            vCMi = (mv1i + mv2i) / (m1 + m2);

            fcti = m1 / job->dt * (vCMi - mv1i / m1);
            interaction_force.row(i) = fcti.transpose();
            job->bodies[b1].nodes.f.row(i) += fcti.transpose();
            job->bodies[b2].nodes.f.row(i) -= fcti.transpose();
        } else {
            interaction_force.row(i).setZero();
        }
    }

    return;
}


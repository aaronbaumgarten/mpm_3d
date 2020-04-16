//
// Created by aaron on 4/15/20.
// contact_linked.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "contacts.hpp"

/*----------------------------------------------------------------------------*/
//initialize contact, requires that job and bodies be setup first
void ContactLinked::init(Job* job){
    //check that contact properties are set
    if ((str_props.size() < 2 && int_props.size() < 2)){
        //need bodies
        std::cout << int_props.size() << ", " << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined ({body_1,body_2}).\n",
                __FILE__, __func__);
        exit(0);
    } else {

        //set body ids by name
        if (str_props.size() >= 2){
            for (int i=0;i<bodyIDs.size();i++) {
                for (int b = 0; b < job->bodies.size(); b++) {
                    if (str_props[i].compare(job->bodies[b]->name) == 0){
                        bodyIDs[i] = b;
                        break;
                    }
                }
            }
        }

        // or set body ids by int
        for (int i=0;i<bodyIDs.size();i++) {
            if (bodyIDs[i] < 0){
                if (int_props.size() == 2) {
                    bodyIDs = int_props;
                } else {
                    std::cout << fp64_props.size() << ", " << int_props.size() << ", " << str_props.size() << "\n";
                    fprintf(stderr,
                            "%s:%s: Need at least 2 properties defined ({body_1,body_2}).\n",
                            __FILE__, __func__);
                    exit(0);
                }
                break;
            }
        }
        contact_force = KinematicVectorArray(job->bodies[bodyIDs[0]]->nodes->x.size(),job->JOB_TYPE);
        contact_normal = KinematicVectorArray(job->bodies[bodyIDs[0]]->nodes->x.size(),job->JOB_TYPE);

        printf("Contact properties ({%i, %i}).\n", bodyIDs[0], bodyIDs[1]);
    }

    //loop over str-props and assign relevant flags
    std::vector<std::string> options = {"SLIP_BOUNDARY"};
    for (int i=0; i<str_props.size(); i++){
        switch (Parser::findStringID(options, str_props[i])){
            case 0:
                //USE_REDUCED_QUADRATURE
                SLIP = true;
                std::cout << "ContactLinked using slip boundary." << std::endl;
                break;
            default:
                //do nothing
                break;
        }
    }

    std::cout << "Contact Initialized: [" << id << "]." << std::endl;

    return;
}


/*----------------------------------------------------------------------------*/
//generate contact rules
void ContactLinked::generateRules(Job* job){
    //set normal for problem
    if (SLIP) {
        //use normal from body 1
        if (job->JOB_TYPE == job->JOB_AXISYM) {
            Eigen::VectorXd pval = Eigen::VectorXd(job->bodies[bodyIDs[0]]->points->x.size());
            for (int i = 0; i < pval.rows(); i++) {
                pval(i) = job->bodies[bodyIDs[0]]->points->m(i) / job->bodies[bodyIDs[0]]->points->x(i, 0);
            }
            //calculate normal using 2D integral of density
            contact_normal = job->bodies[bodyIDs[0]]->gradS * pval;
        } else {
            contact_normal = job->bodies[bodyIDs[0]]->gradS * job->bodies[bodyIDs[0]]->points->m;
        }
        for (int i = 0; i < contact_normal.size(); i++) {
            //normalize
            contact_normal(i) /= contact_normal(i).norm();
        }
    }
    return;
}


/*----------------------------------------------------------------------------*/
//apply contact rules
void ContactLinked::applyRules(Job* job, int SPEC){
    KinematicVector normal(job->JOB_TYPE);
    double m1, m2;
    double fn1i, ft1i;
    KinematicVector mv1i(job->JOB_TYPE);
    KinematicVector mv2i(job->JOB_TYPE);
    KinematicVector vCMi(job->JOB_TYPE);
    KinematicVector fcti(job->JOB_TYPE);
    KinematicVector s1i(job->JOB_TYPE);

    KinematicVector tmpVec(job->JOB_TYPE);

    int b1 = bodyIDs[0];
    int b2 = bodyIDs[1];

    //look for contacts if there are two bodies
    for (int i = 0; i < contact_force.size(); i++) {
        //test every node for contact
        if (job->bodies[b1]->nodes->m[i] > 0 && job->bodies[b2]->nodes->m[i] > 0) {

            fcti.setZero();
            normal = contact_normal(i);

            //determine 'center of mass' velocity
            m1 = job->bodies[b1]->nodes->m[i];
            m2 = job->bodies[b2]->nodes->m[i];

            if (SPEC == Contact::IMPLICIT) {
                mv1i = (job->bodies[b1]->nodes->mx_t(i) + job->dt * job->bodies[b1]->nodes->f(i));
                mv2i = (job->bodies[b2]->nodes->mx_t(i) + job->dt * job->bodies[b2]->nodes->f(i));
            } else if (SPEC == Contact::EXPLICIT){
                mv1i = job->bodies[b1]->nodes->mx_t(i);
                mv2i = job->bodies[b2]->nodes->mx_t(i);
            } else {
                std::cerr << "ERROR: Unknown SPEC in contact_huang.so: " << SPEC << "!" << std::endl;
                return;
            }

            vCMi = (mv1i + mv2i) / (m1 + m2);

            if (!SLIP) {
                //two bodies move together
                tmpVec = (m1 * vCMi - mv1i) / job->dt; //force on body 1
                job->bodies[b1]->nodes->f(i) += tmpVec;
                job->bodies[b2]->nodes->f(i) -= tmpVec;
            } else {
                //determine normal force
                fn1i = (m1 * vCMi.dot(normal) - mv1i.dot(normal)) / job->dt;

                //set contact forces
                job->bodies[b1]->nodes->f(i) += fn1i * normal;
                job->bodies[b2]->nodes->f(i) -= fn1i * normal;
            }

            contact_force[i] = tmpVec;
        } else {
            //zero reported force
            contact_force[i].setZero();
        }
    }

    return;
}


/*----------------------------------------------------------------------------*/
//write to serializer and save/load
void ContactLinked::writeFrame(Job* job, Serializer* serializer){
    serializer->writeVectorArray(contact_force,("contact_force_by_" + job->bodies[bodyIDs[1]]->name + "_on_" + job->bodies[bodyIDs[0]]->name));
    return;
}

std::string ContactLinked::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

int ContactLinked::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}
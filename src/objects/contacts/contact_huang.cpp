//
// Created by aaron on 5/15/18.
// contact_huang.cpp
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
void ContactHuang::init(Job* job){
    //check that contact properties are set
    if (fp64_props.size() < 1 || (str_props.size() < 2 && int_props.size() < 2)){
        //need to coefficient of friction and bodies
        std::cout << fp64_props.size() << ", " << int_props.size() << ", " << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 3 properties defined ({mu_f},{body_1,body_2}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //set coeff of friction
        mu_f = fp64_props[0];

        //set body ids by name
        if (str_props.size() == 2){
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
                            "%s:%s: Need at least 3 properties defined ({mu_f},{body_1,body_2}).\n",
                            __FILE__, __func__);
                    exit(0);
                }
                break;
            }
        }

        contact_normal = KinematicVectorArray(job->bodies[bodyIDs[0]]->nodes->x.size(),job->JOB_TYPE);
        contact_force = KinematicVectorArray(job->bodies[bodyIDs[0]]->nodes->x.size(),job->JOB_TYPE);

        printf("Contact properties (mu_f = %g, {%i, %i}).\n",
               mu_f, bodyIDs[0], bodyIDs[1]);
    }

    std::cout << "Contact Initialized: [" << id << "]." << std::endl;

    return;
}


/*----------------------------------------------------------------------------*/
//generate contact rules
void ContactHuang::generateRules(Job* job){
    //set normal for problem
    //use normal from body 1
    contact_normal = job->bodies[bodyIDs[0]]->gradS * job->bodies[bodyIDs[0]]->points->m;
    for (int i=0;i<contact_normal.size();i++){
        //normalize
        contact_normal(i) /= contact_normal(i).norm();
    }

    return;
}


/*----------------------------------------------------------------------------*/
//apply contact rules
void ContactHuang::applyRules(Job* job, int SPEC){
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
    for (int i = 0; i < contact_normal.size(); i++) {
        //test every node for contact
        if (job->bodies[b1]->nodes->m[i] > 0 && job->bodies[b2]->nodes->m[i] > 0) {

            normal = contact_normal(i);

            //determine 'center of mass' velocity
            m1 = job->bodies[b1]->nodes->m[i];
            m2 = job->bodies[b2]->nodes->m[i];

            if (SPEC == Contact::IMPLICIT) {
                mv1i = (job->bodies[b1]->nodes->mx_t(i) + job->dt * job->bodies[b1]->nodes->f(i));
                mv2i = (job->bodies[b2]->nodes->mx_t(i) + job->dt * job->bodies[b2]->nodes->f(i));

                vCMi = (mv1i + mv2i) / (m1 + m2);

                //check if converging
                if ((mv1i / m1 - vCMi).dot(normal) > 0) {
                    //determine normal force
                    //fn1i = m1 * m2 / (job->dt * (m1 + m2)) * (mv2i.dot(n1i) / m2 - mv1i.dot(n1i) / m1);
                    fn1i = m1 / job->dt * (vCMi.dot(normal) - mv1i.dot(normal) / m1);

                    //determine shear force and shear vector
                    s1i = m1 / job->dt * (vCMi - mv1i / m1) - fn1i * normal;
                    ft1i = sqrt(s1i.dot(s1i));
                    s1i /= ft1i;

                    //add forces
                    fcti = std::min(0.0, fn1i) * normal + std::min(mu_f * std::abs(fn1i), std::abs(ft1i)) * s1i;

                    //set contact forces
                    job->bodies[b1]->nodes->f(i) += fcti;
                    job->bodies[b2]->nodes->f(i) -= fcti;
                }
            } else if (SPEC == Contact::EXPLICIT){
                mv1i = job->bodies[b1]->nodes->mx_t(i);
                mv2i = job->bodies[b2]->nodes->mx_t(i);
                //mv1i << job->bodies[b1].nodes.mx_t.row(i).transpose();
                //mv2i << job->bodies[b2].nodes.mx_t.row(i).transpose();

                vCMi = (mv1i + mv2i) / (m1 + m2);

                //check if converging
                if ((mv1i / m1 - vCMi).dot(normal) > 0) {
                    //determine normal force
                    //fn1i = m1 * m2 / (job->dt * (m1 + m2)) * (mv2i.dot(n1i) / m2 - mv1i.dot(n1i) / m1);
                    fn1i = m1 / job->dt * (vCMi.dot(normal) - mv1i.dot(normal) / m1);

                    //determine shear force and shear vector
                    s1i = m1 / job->dt * (vCMi - mv1i / m1) - fn1i * normal;
                    ft1i = sqrt(s1i.dot(s1i));
                    s1i /= ft1i;

                    //add forces
                    fcti = std::min(0.0, fn1i) * normal + std::min(mu_f * std::abs(fn1i), std::abs(ft1i)) * s1i;

                    //set contact forces
                    job->bodies[b1]->nodes->f(i) += fcti;
                    job->bodies[b2]->nodes->f(i) -= fcti;
                }
            } else {
                std::cerr << "ERROR: Unknown SPEC in contact_huang.so: " << SPEC << "!" << std::endl;
                return;
            }
            contact_force[i] = fcti;
        } else {
            //zero reported force
            contact_force[i].setZero();
        }
    }

    return;
}


/*----------------------------------------------------------------------------*/
//write to serializer and save/load
void ContactHuang::writeFrame(Job* job, Serializer* serializer){
    serializer->writeVectorArray(contact_normal,("contact_normal_" + job->bodies[bodyIDs[0]]->name + "_" + job->bodies[bodyIDs[1]]->name));
    serializer->writeVectorArray(contact_force,("contact_force_by_" + job->bodies[bodyIDs[1]]->name + "_on_" + job->bodies[bodyIDs[0]]->name));
    return;
}

std::string ContactHuang::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

int ContactHuang::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}
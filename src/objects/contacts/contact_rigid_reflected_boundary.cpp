//
// Created by aaron on 11/29/18.
// contact_rigid_reflected_boundary.cpp
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
void ContactRigid_ReflectedBoundary::init(Job* job){
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

        Lx = KinematicVector(job->JOB_TYPE);
        Lx.setZero();
        for (int i=0; i < job->bodies[bodyIDs[0]]->nodes->x.size(); i++){
            for (int pos=0; pos < job->bodies[bodyIDs[0]]->nodes->x.DIM; pos++){
                if (job->bodies[bodyIDs[0]]->nodes->x(i,pos) > Lx(pos)){
                    Lx(pos) = job->bodies[bodyIDs[0]]->nodes->x(i,pos);
                }
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
void ContactRigid_ReflectedBoundary::generateRules(Job* job){
    //set normal for problem
    //use normal from body 1
    if (job->JOB_TYPE == job->JOB_AXISYM){
        Eigen::VectorXd pval = Eigen::VectorXd(job->bodies[bodyIDs[0]]->points->x.size());
        for (int i=0; i<pval.rows(); i++){
            pval(i) = job->bodies[bodyIDs[0]]->points->m(i)/job->bodies[bodyIDs[0]]->points->x(i,0);
        }
        //calculate normal using 2D integral of density
        contact_normal = job->bodies[bodyIDs[0]]->gradS * pval;
    } else {
        contact_normal = job->bodies[bodyIDs[0]]->gradS * job->bodies[bodyIDs[0]]->points->m;
    }

    for (int i=0;i<contact_normal.size();i++){
        if (job->bodies[bodyIDs[0]]->nodes->m[i] > 0) {
            for (int pos = 0; pos < job->DIM; pos++) {
                if (job->bodies[bodyIDs[0]]->nodes->x(i, pos) == 0 ||
                    job->bodies[bodyIDs[0]]->nodes->x(i, pos) == Lx(pos)) {
                    contact_normal(i, pos) = 0;
                }
            }
            //normalize
            contact_normal(i) /= contact_normal(i).norm();
        }
    }

    return;
}

/*----------------------------------------------------------------------------*/
//apply contact rules
void ContactRigid_ReflectedBoundary::applyRules(Job* job, int SPEC){
    KinematicVector normal(job->JOB_TYPE);
    double m1, m2;
    double fn2i, ft2i;
    KinematicVector mv1i(job->JOB_TYPE);
    KinematicVector mv2i(job->JOB_TYPE);
    KinematicVector vCMi(job->JOB_TYPE);
    KinematicVector fcti(job->JOB_TYPE);
    KinematicVector s2i(job->JOB_TYPE);

    KinematicVector tmpVec(job->JOB_TYPE);

    int b1 = bodyIDs[0];
    int b2 = bodyIDs[1];

    //look for contacts if there are two bodies
    for (int i = 0; i < contact_normal.size(); i++) {
        //test every node for contact
        if (job->bodies[b1]->nodes->m[i] > 0 && job->bodies[b2]->nodes->m[i] > 0) {

            normal = -contact_normal(i);

            //determine velocity
            m1 = job->bodies[b1]->nodes->m[i];
            m2 = job->bodies[b2]->nodes->m[i];

            mv1i = job->bodies[b1]->nodes->mx_t(i);
            if (SPEC == Contact::EXPLICIT){
                mv2i = job->bodies[b2]->nodes->mx_t(i);
            } else if (SPEC == Contact::IMPLICIT) {
                mv2i = (job->bodies[b2]->nodes->mx_t(i) + job->dt * job->bodies[b2]->nodes->f(i));
            } else {
                std::cerr << "ERROR: Unknown SPEC in contact_huang.so: " << SPEC << "!" << std::endl;
                return;
            }
            vCMi = mv1i/m1;

            //check if converging
            if ((mv2i / m2 - vCMi).dot(normal) > 0) {
                //determine normal force
                //fn1i = m1 * m2 / (job->dt * (m1 + m2)) * (mv2i.dot(n1i) / m2 - mv1i.dot(n1i) / m1);
                fn2i = m2 / job->dt * (vCMi.dot(normal) - mv2i.dot(normal) / m2);

                //determine shear force and shear vector
                s2i = m2 / job->dt * (vCMi - mv2i / m2) - fn2i * normal;
                ft2i = s2i.norm();
                if (ft2i != 0) {
                    s2i /= ft2i;
                }

                //add forces
                fcti = std::min(0.0, fn2i) * normal + std::min(mu_f * std::abs(fn2i), std::abs(ft2i)) * s2i;

                //set contact forces
                job->bodies[b2]->nodes->f(i) += fcti;
                job->bodies[b1]->nodes->f(i) -= fcti;

                contact_force[i] = -fcti;
            } else {
                contact_force[i] = 0;
            }
        } else {
            //zero reported force
            contact_force[i].setZero();
        }
    }

    return;
}
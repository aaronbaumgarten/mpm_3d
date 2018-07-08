//
// Created by aaron on 6/22/18.
// contact_slurry.cpp
//

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
void SlurryContact::init(Job* job){
    //check that contact properties are set
    if (fp64_props.size() < 1 || (str_props.size() < 3 && int_props.size() < 3)){
        //need to coefficient of friction and bodies
        std::cout << fp64_props.size() << ", " << int_props.size() << ", " << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 4 properties defined ({mu_f},{solid_body,granular_body,fluid_body}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //set coeff of friction
        mu_f = fp64_props[0];

        //set body ids by name
        if (str_props.size() == 3){
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
                if (int_props.size() == 3) {
                    bodyIDs = int_props;
                } else {
                    std::cout << fp64_props.size() << ", " << int_props.size() << ", " << str_props.size() << "\n";
                    fprintf(stderr,
                            "%s:%s: Need at least 4 properties defined ({mu_f},{solid_body,granular_body,fluid_body}).\n",
                            __FILE__, __func__);
                    exit(0);
                }
                break;
            }
        }

        contact_normal = KinematicVectorArray(job->bodies[bodyIDs[0]]->nodes->x.size(),job->JOB_TYPE);
        contact_force = KinematicVectorArray(job->bodies[bodyIDs[0]]->nodes->x.size(),job->JOB_TYPE);

        printf("Contact properties (mu_f = %g, {%i, %i, %i}).\n",
               mu_f, bodyIDs[0], bodyIDs[1], bodyIDs[2]);
    }

    std::cout << "Contact Initialized: [" << id << "]." << std::endl;

    return;
}


/*----------------------------------------------------------------------------*/
//generate contact rules
void SlurryContact::generateRules(Job* job){
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
void SlurryContact::applyRules(Job* job, int SPEC){
    KinematicVector normal(job->JOB_TYPE);
    double m1, m2, m3;
    double fn1i, ft1i;
    double fn2i, ft2i;
    KinematicVector mv1i(job->JOB_TYPE);
    KinematicVector mv2i(job->JOB_TYPE);
    KinematicVector mv3i(job->JOB_TYPE);
    KinematicVector f1i(job->JOB_TYPE);
    KinematicVector f2i(job->JOB_TYPE);
    KinematicVector s1i(job->JOB_TYPE);
    KinematicVector s2i(job->JOB_TYPE);

    KinematicVector tmpVec(job->JOB_TYPE);

    int b1 = bodyIDs[0];
    int b2 = bodyIDs[1];
    int b3 = bodyIDs[3];

    //look for contacts if there are two bodies
    for (int i = 0; i < contact_normal.size(); i++) {
        //test every node for contact
        if (job->bodies[b1]->nodes->m[i] > 0 && (job->bodies[b2]->nodes->m[i] > 0 || job->bodies[b3]->nodes->m[i] > 0)){
            //contact must involve body 1
            normal = contact_normal(i);
            m1 = job->bodies[b1]->nodes->m[i];
            m2 = job->bodies[b2]->nodes->m[i];
            m3 = job->bodies[b3]->nodes->m[i];

            if (SPEC == Contact::IMPLICIT){
                mv1i = (job->bodies[b1]->nodes->mx_t(i) + job->dt * job->bodies[b1]->nodes->f(i));
                mv2i = (job->bodies[b2]->nodes->mx_t(i) + job->dt * job->bodies[b2]->nodes->f(i));
                mv3i = (job->bodies[b3]->nodes->mx_t(i) + job->dt * job->bodies[b3]->nodes->f(i));
            } else if (SPEC == Contact::EXPLICIT){
                mv1i = job->bodies[b1]->nodes->mx_t(i);
                mv2i = job->bodies[b2]->nodes->mx_t(i);
                mv3i = job->bodies[b3]->nodes->mx_t(i);
            } else {
                std::cerr << "ERROR: Unknown SPEC in contact_huang.so: " << SPEC << "!" << std::endl;
                return;
            }

            if ((m2 > 0 && (mv1i/m1 - mv2i/m2).dot(normal) > 0) && (m3 > 0 && (mv1i/m1 - mv3i/m3).dot(normal) > 0)){
                //both contact forces are non-zero
                fn1i = (m1*m2*(mv1i/m1 - mv2i/m2).dot(normal) + m2*m3*(mv3i/m3 - mv2i/m2).dot(normal))/(job->dt * (m1 + m2 + m3));
                fn2i = (m1*m3*(mv1i/m1 - mv3i/m3).dot(normal) + m3*m2*(mv2i/m2 - mv3i/m3).dot(normal))/(job->dt * (m1 + m2 + m3));
            } else if (m2 > 0 && (mv1i/m1 - mv2i/m2).dot(normal) > 0) {
                //only body 2 is in contact
                fn1i = (m1*m2*(mv1i/m1 - mv2i/m2).dot(normal))/(job->dt * (m1 + m2));
                fn2i = 0;
            } else if (m3 > 0 && (mv1i/m1 - mv3i/m3).dot(normal) > 0) {
                //only body 3 is in contact
                fn1i = 0;
                fn2i = (m1*m3*(mv1i/m1 - mv3i/m3).dot(normal))/(job->dt * (m1 + m3));
            } else {
                //no contact
                fn1i = 0;
                fn2i = 0;
            }

            //granular tangential force
            if (m2 > 0) {
                s1i = (mv1i - (fn1i + fn2i) * normal * job->dt) / m1 -
                      (mv2i + fn1i * normal * job->dt) / m2; //new tangent velocity
                s1i = s1i / s1i.norm(); //normalized negative tangent direction
                ft1i = mu_f * fn1i; //tangential force

                //store resulting relative velocity
                tmpVec = (mv1i - (fn1i + fn2i)*normal*job->dt - ft1i*s1i*job->dt)/m1 -
                         (mv2i + fn1i*normal*job->dt + ft1i*s1i*job->dt)/m2;

                if (tmpVec.dot(s1i) < 0){
                    //then they should stick
                    ft1i = (m2*(mv1i - (fn1i + fn2i)*normal*job->dt) - m1*(mv2i + fn1i*normal*job->dt)).norm() / (m1 + m2);
                }
            } else {
                s1i.setZero();
                ft1i = 0;
            }

            //fluid tangential force
            /*if (m3 > 0) {
                s2i = (mv3i + fn2i * normal * job->dt) / m3 -
                      (mv1i - (fn1i + fn2i) * normal * job->dt) / m1; //new tangent velocity
                s2i = -s2i / s2i.norm(); //normalized negative tangent direction
            } else {
                s2i.setZero();
                ft2i = 0;
            }*/
            //zero for now
            s2i.setZero();
            ft2i = 0;

            if (!std::isfinite(s1i.norm())){
                s1i.setZero();
                ft1i = 0;
            }

            f1i = fn1i*normal + ft1i*s1i;
            f2i = fn1i*normal + ft2i*s2i;

            job->bodies[b1]->nodes->f[i] -= f1i + f2i;
            job->bodies[b2]->nodes->f[i] += f1i;
            job->bodies[b3]->nodes->f[i] += f2i;
            contact_force[i] = -f1i - f2i;
        } else {
            //zero reported force
            contact_force[i].setZero();
        }
    }

    return;
}


/*----------------------------------------------------------------------------*/
//write to serializer and save/load
void SlurryContact::writeFrame(Job* job, Serializer* serializer){
    serializer->writeVectorArray(contact_normal,("contact_normal_" + job->bodies[bodyIDs[0]]->name + "_" + job->bodies[bodyIDs[1]]->name));
    serializer->writeVectorArray(contact_force,("contact_force_by_" + job->bodies[bodyIDs[1]]->name + "_and_" + job->bodies[bodyIDs[2]]->name +"_on_" + job->bodies[bodyIDs[0]]->name));
    return;
}

std::string SlurryContact::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

int SlurryContact::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}
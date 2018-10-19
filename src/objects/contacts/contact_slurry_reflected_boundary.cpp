//
// Created by aaron on 7/8/18.
// contact_slurry_reflected_boundary.cpp
//

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
void SlurryContact_ReflectedBoundary::init(Job* job){
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

        Lx = KinematicVector(job->JOB_TYPE);
        Lx.setZero();
        for (int i=0; i < job->bodies[bodyIDs[0]]->nodes->x.size(); i++){
            for (int pos=0; pos < job->bodies[bodyIDs[0]]->nodes->x.DIM; pos++){
                if (job->bodies[bodyIDs[0]]->nodes->x(i,pos) > Lx(pos)){
                    Lx(pos) = job->bodies[bodyIDs[0]]->nodes->x(i,pos);
                }
            }
        }

        printf("Contact properties (mu_f = %g, {%i, %i, %i}).\n",
               mu_f, bodyIDs[0], bodyIDs[1], bodyIDs[2]);
    }

    std::cout << "Contact Initialized: [" << id << "]." << std::endl;

    return;
}


/*----------------------------------------------------------------------------*/
//generate contact rules
void SlurryContact_ReflectedBoundary::generateRules(Job* job){
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
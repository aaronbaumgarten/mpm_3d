//
// Created by aaron on 12/23/16.
// huang2d.cpp
//

//contact law derived in p huang 2011

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <Eigen/Core>
#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "process.hpp"

#define CONTACT_VERSION_STRING "1.0" __DATE__ " " __TIME__

/* Contact Constants (set by configuration file). */
double mu_f;

extern "C" void contact_init(job_t *job, size_t id);

extern "C" void resolve_contact(job_t *job, size_t id);

/*----------------------------------------------------------------------------*/
void contact_init(job_t *job, size_t id) {

    if (job->contacts[id].num_fp64_props < 1) {
        std::cout << job->contacts[id].num_fp64_props << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 1 properties defined (mu_f).\n",
                __FILE__, __func__);
        exit(0); //replace error code handler from sachith's work with 0
    } else {
        mu_f = job->contacts[id].fp64_props[0];
        /*printf("%s:%s: properties (mu_f = %g).\n",
               __FILE__, __func__, mu_f);*/
        printf("Contact properties (mu_f = %g).\n", mu_f);
    }

    /*printf("%s:%s: (contact version %s) done initializing contact.\n",
           __FILE__, __func__, CONTACT_VERSION_STRING);*/
    std::cout << "Done initializing contact {" << job->contacts[id].bodyIDs[0]  << ", " << job->contacts[id].bodyIDs[1] << "}.\n";
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void resolve_contact(job_t *job, size_t id) {
    size_t b1 = job->contacts[id].bodyIDs[0];
    size_t b2 = job->contacts[id].bodyIDs[1];
    //look for contacts if there are two bodies
    for (size_t i = 0; i < job->num_nodes; i++) {
        //test every node for contact
        if (job->bodies[b1].nodes.m[i] > TOL && job->bodies[b2].nodes.m[i] > TOL) {
            //use normal from body 1
            Eigen::Vector3d n1i;
            n1i << job->bodies[b1].nodes.contact_normal_x[i],
                    job->bodies[b1].nodes.contact_normal_y[i],
                    job->bodies[b1].nodes.contact_normal_z[i];
            //enforce unit length
            n1i /= sqrt(n1i.dot(n1i));

            //determine 'center of mass' velocity
            double m1 = job->bodies[b1].nodes.m[i];
            double m2 = job->bodies[b2].nodes.m[i];
            Eigen::Vector3d mv1i;
            Eigen::Vector3d mv2i;
            Eigen::Vector3d vCMi;
            mv1i << job->bodies[b1].nodes.contact_mx_t[i],
                    job->bodies[b1].nodes.contact_my_t[i],
                    job->bodies[b1].nodes.contact_mz_t[i];
            mv2i << job->bodies[b2].nodes.contact_mx_t[i],
                    job->bodies[b2].nodes.contact_my_t[i],
                    job->bodies[b2].nodes.contact_mz_t[i];
            vCMi = (mv1i + mv2i) / (m1 + m2);

            //check if converging
            if ((mv1i/m1 - vCMi).dot(n1i) > 0 ) {
                //determine normal force
                double fn1i;
                fn1i = m1 * m2 / (job->dt * (m1 + m2)) * (mv2i.dot(n1i) / m2 - mv1i.dot(n1i) / m1);

                //determine shear force and shear vector
                double ft1i;
                Eigen::Vector3d s1i;
                s1i = m1 / job->dt * (vCMi - mv1i / m1) - fn1i * n1i;
                ft1i = sqrt(s1i.dot(s1i));
                s1i /= ft1i;

                //add forces
                Eigen::Vector3d fcti;
                fcti = std::min(0.0, fn1i) * n1i + std::min(mu_f * std::abs(fn1i), std::abs(ft1i)) * s1i;

                //set contact forces
                job->bodies[b1].nodes.contact_fx[i] = fcti[0];
                job->bodies[b1].nodes.contact_fy[i] = fcti[1];
                job->bodies[b1].nodes.contact_fz[i] = fcti[2];

                job->bodies[b2].nodes.contact_fx[i] = -fcti[0];
                job->bodies[b2].nodes.contact_fy[i] = -fcti[1];
                job->bodies[b2].nodes.contact_fz[i] = -fcti[2];

                if (job->use_implicit == 0) {
                    //adjust nodal velocities for non-penetration
                    mv1i = mv1i - n1i.dot(mv1i - m1 * vCMi) * n1i;
                    mv2i = mv2i - n1i.dot(mv2i - m2 * vCMi) * n1i;

                    job->bodies[b1].nodes.contact_mx_t[i] = mv1i[0];
                    job->bodies[b1].nodes.contact_my_t[i] = mv1i[1];
                    job->bodies[b1].nodes.contact_mz_t[i] = mv1i[2];

                    job->bodies[b1].nodes.contact_x_t[i] = mv1i[0] / m1;
                    job->bodies[b1].nodes.contact_y_t[i] = mv1i[1] / m1;
                    job->bodies[b1].nodes.contact_z_t[i] = mv1i[2] / m1;

                    job->bodies[b2].nodes.contact_mx_t[i] = mv2i[0];
                    job->bodies[b2].nodes.contact_my_t[i] = mv2i[1];
                    job->bodies[b2].nodes.contact_mz_t[i] = mv2i[2];

                    job->bodies[b2].nodes.contact_x_t[i] = mv2i[0] / m2;
                    job->bodies[b2].nodes.contact_y_t[i] = mv2i[1] / m2;
                    job->bodies[b2].nodes.contact_z_t[i] = mv2i[2] / m2;
                }
            }
        }
    }
    return;
}
/*----------------------------------------------------------------------------*/
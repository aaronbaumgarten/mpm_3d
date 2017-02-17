//
// Created by aaron on 2/10/17.
// slurry.cpp
//

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
double drag, grains_rho;

extern "C" void contact_init(job_t *job, size_t id);

extern "C" void resolve_contact(job_t *job, size_t id);

/*----------------------------------------------------------------------------*/
void contact_init(job_t *job, size_t id) {

    if (job->contacts[id].num_fp64_props < 2) {
        std::cout << job->contacts[id].num_fp64_props << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 properties defined (drag, grains_rho).\n",
                __FILE__, __func__);
        exit(0); //replace error code handler from sachith's work with 0
    } else {
        drag = job->contacts[id].fp64_props[0];
        grains_rho = job->contacts[id].fp64_props[1];
        /*printf("%s:%s: properties (mu_f = %g).\n",
               __FILE__, __func__, mu_f);*/
        printf("Contact properties (drag = %g, grains_rho = %g).\n", drag, grains_rho);
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
    //assume that volume fraction given by solid portion (first body) and cast to nodes
    Eigen::VectorXd pvec(job->bodies[b1].p);
    pvec = job->bodies[b1].particles.m.array()/job->bodies[b1].particles.v.array();
    Eigen::VectorXd packing_fraction(job->num_nodes);
    packing_fraction = job->bodies[b1].Phi * pvec;
    packing_fraction *= 1.0/grains_rho; //packing fraction

    //remove body force from initial nodal fx (second body)
    Eigen::VectorXd fluid_fx(job->num_nodes);
    Eigen::VectorXd fluid_fy(job->num_nodes);
    Eigen::VectorXd fluid_fz(job->num_nodes);
    
    fluid_fx.setZero();
    fluid_fy.setZero();
    fluid_fz.setZero();

    pvec = job->bodies[b2].particles.v.array() * job->bodies[b2].particles.T.col(XX).array();
    fluid_fx += job->bodies[b2].gradPhiX*pvec;

    pvec = job->bodies[b2].particles.v.array() * job->bodies[b2].particles.T.col(XY).array();
    fluid_fx += job->bodies[b2].gradPhiY*pvec;
    fluid_fy += job->bodies[b2].gradPhiX*pvec;

    pvec = job->bodies[b2].particles.v.array() * job->bodies[b2].particles.T.col(XZ).array();
    fluid_fx += job->bodies[b2].gradPhiZ*pvec;
    fluid_fz += job->bodies[b2].gradPhiX*pvec;

    pvec = job->bodies[b2].particles.v.array() * job->bodies[b2].particles.T.col(YY).array();
    fluid_fy += job->bodies[b2].gradPhiY*pvec;

    pvec = job->bodies[b2].particles.v.array() * job->bodies[b2].particles.T.col(YZ).array();
    fluid_fy += job->bodies[b2].gradPhiZ*pvec;
    fluid_fz += job->bodies[b2].gradPhiY*pvec;

    pvec = job->bodies[b2].particles.v.array() * job->bodies[b2].particles.T.col(ZZ).array();
    fluid_fz += job->bodies[b2].gradPhiZ*pvec;


    Eigen::Vector3d fcti;
    for (size_t i = 0; i < job->num_nodes; i++) {
        //test every node for contact
        if (job->bodies[b1].nodes.m[i] > TOL && job->bodies[b2].nodes.m[i] > TOL) {
            //determine 'center of mass' velocity
            double m1 = job->bodies[b1].nodes.m[i];
            double m2 = job->bodies[b2].nodes.m[i];
            Eigen::Vector3d mv1i;
            Eigen::Vector3d mv2i;
            Eigen::Vector3d vCMi;
            if (job->use_implicit == 1){
                mv1i << job->bodies[b1].nodes.mx_t_k[i] + job->dt * job->bodies[b1].nodes.contact_fx[i],
                        job->bodies[b1].nodes.my_t_k[i] + job->dt * job->bodies[b1].nodes.contact_fy[i],
                        job->bodies[b1].nodes.mz_t_k[i] + job->dt * job->bodies[b1].nodes.contact_fz[i];
                mv2i << job->bodies[b2].nodes.mx_t_k[i] + job->dt * job->bodies[b2].nodes.contact_fx[i],
                        job->bodies[b2].nodes.my_t_k[i] + job->dt * job->bodies[b2].nodes.contact_fy[i],
                        job->bodies[b2].nodes.mz_t_k[i] + job->dt * job->bodies[b2].nodes.contact_fz[i];
            } else {
                mv1i << job->bodies[b1].nodes.contact_mx_t[i] + job->dt * job->bodies[b1].nodes.contact_fx[i],
                        job->bodies[b1].nodes.contact_my_t[i] + job->dt * job->bodies[b1].nodes.contact_fy[i],
                        job->bodies[b1].nodes.contact_mz_t[i] + job->dt * job->bodies[b1].nodes.contact_fz[i];
                mv2i << job->bodies[b2].nodes.contact_mx_t[i] + job->dt * job->bodies[b2].nodes.contact_fx[i],
                        job->bodies[b2].nodes.contact_my_t[i] + job->dt * job->bodies[b2].nodes.contact_fy[i],
                        job->bodies[b2].nodes.contact_mz_t[i] + job->dt * job->bodies[b2].nodes.contact_fz[i];
            }
            vCMi = (mv1i + mv2i) / (m1 + m2);

            fcti = (1-packing_fraction[i])*(1-packing_fraction[i])*drag*(mv2i/m2 - mv1i/m1);
            fcti[0] += packing_fraction[i]*fluid_fx[i];
            fcti[1] += packing_fraction[i]*fluid_fy[i];
            fcti[2] += packing_fraction[i]*fluid_fz[i];
            
            //set contact forces
            job->bodies[b1].nodes.contact_fx[i] += fcti[0];
            job->bodies[b1].nodes.contact_fy[i] += fcti[1];
            job->bodies[b1].nodes.contact_fz[i] += fcti[2];

            job->bodies[b2].nodes.contact_fx[i] -= fcti[0];
            job->bodies[b2].nodes.contact_fy[i] -= fcti[1];
            job->bodies[b2].nodes.contact_fz[i] -= fcti[2];
        }
    }
    return;
}
/*----------------------------------------------------------------------------*/
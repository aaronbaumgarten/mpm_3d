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
double drag, grains_rho, k, mu_w;

extern "C" void contact_init(job_t *job, size_t id);

extern "C" void resolve_contact(job_t *job, size_t id);

/*----------------------------------------------------------------------------*/
void contact_init(job_t *job, size_t id) {

    if (job->contacts[id].num_fp64_props < 4) {
        std::cout << job->contacts[id].num_fp64_props << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 4 properties defined (drag, grains_rho, k, mu_w).\n",
                __FILE__, __func__);
        exit(0); //replace error code handler from sachith's work with 0
    } else {
        drag = job->contacts[id].fp64_props[0];
        grains_rho = job->contacts[id].fp64_props[1];
        k = job->contacts[id].fp64_props[2];
        mu_w = job->contacts[id].fp64_props[3];
        /*printf("%s:%s: properties (mu_f = %g).\n",
               __FILE__, __func__, mu_f);*/
        printf("Contact properties (drag = %g, grains_rho = %g, k = %g, mu_w = %g).\n", drag, grains_rho, k, mu_w);
    }

    /*printf("%s:%s: (contact version %s) done initializing contact.\n",
           __FILE__, __func__, CONTACT_VERSION_STRING);*/
    std::cout << "Done initializing contact {" << job->contacts[id].bodyIDs[0]  << ", " << job->contacts[id].bodyIDs[1] << "}.\n";
    return;
}

/*----------------------------------------------------------------------------*/
void resolve_contact(job_t *job, size_t id) {
    //solid body 1
    //liquid body 2
    size_t solid_body_id = job->contacts[id].bodyIDs[0];
    size_t liquid_body_id = job->contacts[id].bodyIDs[1];

    //for giggles, try assigning contact here
    job->boundary.generate_dirichlet_bcs(job);

    //assume that porosity is given by solid skeleton and cast to nodes
    Eigen::VectorXd pvec(job->bodies[solid_body_id].p);
    Eigen::VectorXd n(job->num_nodes);
    //Eigen::VectorXd gradnx(job->num_nodes);
    //Eigen::VectorXd gradny(job->num_nodes);
    //Eigen::VectorXd gradnz(job->num_nodes);

    //find presure gradient from fluid
    Eigen::VectorXd gradpx(job->num_nodes);
    Eigen::VectorXd gradpy(job->num_nodes);
    Eigen::VectorXd gradpz(job->num_nodes);

    pvec = -job->bodies[solid_body_id].particles.m.array() / job->bodies[solid_body_id].particles.v.array();
    pvec *= 1.0 / grains_rho;
    pvec += Eigen::VectorXd::Ones(job->bodies[solid_body_id].p); //n = 1-phi
    n = job->bodies[solid_body_id].Phi * pvec;
    Eigen::VectorXd one = job->bodies[solid_body_id].Phi * Eigen::VectorXd::Ones(job->bodies[solid_body_id].p);
    for (size_t i=0;i<job->num_nodes;i++){
        if (job->bodies[solid_body_id].nodes.m[i] != 0){
            n[i] /= one[i];//job->bodies[solid_body_id].nodes.m[i];
        } else {
            n[i] = 1.0;
        }
    }
    //gradnx = job->bodies[solid_body_id].gradPhiX * pvec;
    //gradny = job->bodies[solid_body_id].gradPhiY * pvec;
    //gradnz = job->bodies[solid_body_id].gradPhiZ * pvec;

    pvec = job->bodies[solid_body_id].particles.v.array() * (job->bodies[solid_body_id].particles.T.col(XX) + job->bodies[solid_body_id].particles.T.col(YY) +
           job->bodies[solid_body_id].particles.T.col(ZZ)).array();
    pvec *= -1.0 / 3.0;
    gradpx = job->bodies[solid_body_id].gradPhiX * pvec;
    gradpy = job->bodies[solid_body_id].gradPhiY * pvec;
    gradpz = job->bodies[solid_body_id].gradPhiZ * pvec;

    //nodal contact
    Eigen::Vector3d fsfi;
    Eigen::VectorXd nvecx(job->num_nodes);
    nvecx.setZero();
    Eigen::VectorXd nvecy(job->num_nodes);
    nvecy.setZero();
    Eigen::VectorXd nvecz(job->num_nodes);
    nvecz.setZero();
    for (size_t i = 0; i < job->num_nodes; i++) {
        //test every node for contact
        if (job->bodies[solid_body_id].nodes.m[i] != 0 && job->bodies[liquid_body_id].nodes.m[i] != 0) {
            //determine 'center of mass' velocity
            double m1 = job->bodies[solid_body_id].nodes.m[i];
            double m2 = job->bodies[liquid_body_id].nodes.m[i];
            Eigen::Vector3d mv1i;
            Eigen::Vector3d mv2i;
            Eigen::Vector3d vCMi;
            if (job->use_implicit == 1) {
                /*mv1i << job->bodies[solid_body_id].nodes.mx_t_k[i] + job->dt * job->bodies[solid_body_id].nodes.contact_fx[i],
                        job->bodies[solid_body_id].nodes.my_t_k[i] + job->dt * job->bodies[solid_body_id].nodes.contact_fy[i],
                        job->bodies[solid_body_id].nodes.mz_t_k[i] + job->dt * job->bodies[solid_body_id].nodes.contact_fz[i];
                mv2i << job->bodies[liquid_body_id].nodes.mx_t_k[i] + job->dt * job->bodies[liquid_body_id].nodes.contact_fx[i],
                        job->bodies[liquid_body_id].nodes.my_t_k[i] + job->dt * job->bodies[liquid_body_id].nodes.contact_fy[i],
                        job->bodies[liquid_body_id].nodes.mz_t_k[i] + job->dt * job->bodies[liquid_body_id].nodes.contact_fz[i];*/
                mv1i << job->bodies[solid_body_id].nodes.contact_mx_t[i],
                        job->bodies[solid_body_id].nodes.contact_my_t[i],
                        job->bodies[solid_body_id].nodes.contact_mz_t[i];
                mv2i << job->bodies[liquid_body_id].nodes.contact_mx_t[i],
                        job->bodies[liquid_body_id].nodes.contact_my_t[i],
                        job->bodies[liquid_body_id].nodes.contact_mz_t[i];
            } else {
                mv1i << job->bodies[solid_body_id].nodes.contact_mx_t[i] +
                        job->dt * job->bodies[solid_body_id].nodes.contact_fx[i],
                        job->bodies[solid_body_id].nodes.contact_my_t[i] +
                        job->dt * job->bodies[solid_body_id].nodes.contact_fy[i],
                        job->bodies[solid_body_id].nodes.contact_mz_t[i] +
                        job->dt * job->bodies[solid_body_id].nodes.contact_fz[i];
                mv2i << job->bodies[liquid_body_id].nodes.contact_mx_t[i] +
                        job->dt * job->bodies[liquid_body_id].nodes.contact_fx[i],
                        job->bodies[liquid_body_id].nodes.contact_my_t[i] +
                        job->dt * job->bodies[liquid_body_id].nodes.contact_fy[i],
                        job->bodies[liquid_body_id].nodes.contact_mz_t[i] +
                        job->dt * job->bodies[liquid_body_id].nodes.contact_fz[i];
            }
            vCMi = (mv1i + mv2i) / (m1 + m2);

            //force by solid on fluid
            fsfi = (mv1i / m2 - mv2i / m1) * n[i] * n[i] * mu_w / k;
            fsfi[0] += (1 - n[i]) * gradpx[i];
            fsfi[1] += (1 - n[i]) * gradpy[i];
            fsfi[2] += (1 - n[i]) * gradpz[i];

            //set contact forces
            job->bodies[solid_body_id].nodes.contact_fx[i] -= fsfi[0];
            job->bodies[solid_body_id].nodes.contact_fy[i] -= fsfi[1];
            job->bodies[solid_body_id].nodes.contact_fz[i] -= fsfi[2];

            job->bodies[liquid_body_id].nodes.contact_fx[i] += fsfi[0];
            job->bodies[liquid_body_id].nodes.contact_fy[i] += fsfi[1];
            job->bodies[liquid_body_id].nodes.contact_fz[i] += fsfi[2];

            //determine strainrate variable on nodes
            nvecx[i] = (1 - n[i]) * mv1i[0] / m1 + n[i] * mv2i[0] / m2;
            nvecy[i] = (1 - n[i]) * mv1i[1] / m1 + n[i] * mv2i[1] / m2;
            nvecz[i] = (1 - n[i]) * mv1i[2] / m1 + n[i] * mv2i[2] / m2;
        }
    }

    //determine strainrate on fluid
    Eigen::VectorXd pvec2(job->bodies[liquid_body_id].p);
    job->bodies[liquid_body_id].particles.state.col(10) << (job->bodies[liquid_body_id].gradPhiX.transpose() * nvecx +
                                                           job->bodies[liquid_body_id].gradPhiY.transpose() * nvecy +
                                                           job->bodies[liquid_body_id].gradPhiY.transpose() * nvecz);
    pvec2 = job->bodies[liquid_body_id].Phi.transpose() * n;

    for (size_t p = 0; p < job->bodies[liquid_body_id].p; p++) {
        if (pvec2[p] != 0) {
            job->bodies[liquid_body_id].particles.state(p, 9) = 1;
            job->bodies[liquid_body_id].particles.state(p, 10) /= pvec2[p];
        } else {
            job->bodies[liquid_body_id].particles.state(p, 9) = 0;
        }

        if (!std::isfinite(job->bodies[liquid_body_id].particles.state(p, 10))) {
            std::cout << "SHOOT HER! " << job->bodies[liquid_body_id].particles.state(p, 10) << std::endl;
            exit(0);
        }
    }

    std::cout << "[" << n.maxCoeff() << "," << n.minCoeff() << "]" << std::endl;

    return;
}
/*----------------------------------------------------------------------------*/
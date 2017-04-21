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
double grains_rho, grains_diam, k, mu_w, t_start;

extern "C" void contact_init(job_t *job, size_t id);

extern "C" void resolve_contact(job_t *job, size_t id);

/*----------------------------------------------------------------------------*/
void contact_init(job_t *job, size_t id) {

    if (job->contacts[id].num_fp64_props < 4) {
        std::cout << job->contacts[id].num_fp64_props << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 5 properties defined (drag, grains_rho, k, mu_w, t_start).\n",
                __FILE__, __func__);
        exit(0); //replace error code handler from sachith's work with 0
    } else {
        grains_rho = job->contacts[id].fp64_props[0];
        grains_diam = job->contacts[id].fp64_props[1];
        k = job->contacts[id].fp64_props[2];
        mu_w = job->contacts[id].fp64_props[3];
        t_start = job->contacts[id].fp64_props[4];
        /*printf("%s:%s: properties (mu_f = %g).\n",
               __FILE__, __func__, mu_f);*/
        printf("Contact properties (grains_rho = %g, grains_diam = %g, k = %g, mu_w = %g, t_start = %g).\n", grains_rho, grains_diam, k, mu_w, t_start);
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

    //damping
    if (job->t < t_start){
        job->bodies[solid_body_id].particles.x_t *= 0.99;
        job->bodies[solid_body_id].particles.y_t *= 0.99;
        job->bodies[solid_body_id].particles.z_t *= 0.99;
        job->bodies[liquid_body_id].particles.x_t *= 0.99;
        job->bodies[liquid_body_id].particles.y_t *= 0.99;
        job->bodies[liquid_body_id].particles.z_t *= 0.99;
        //return;
    }

    //assume that porosity is given by solid skeleton and cast to nodes
    Eigen::VectorXd pvec(job->bodies[liquid_body_id].p);
    Eigen::VectorXd n(job->num_nodes);
    Eigen::VectorXd V(job->num_nodes);
    //Eigen::VectorXd gradnx(job->num_nodes);
    //Eigen::VectorXd gradny(job->num_nodes);
    //Eigen::VectorXd gradnz(job->num_nodes);

    //find presure gradient from fluid
    Eigen::VectorXd gradpx(job->num_nodes);
    Eigen::VectorXd gradpy(job->num_nodes);
    Eigen::VectorXd gradpz(job->num_nodes);

    /*
    pvec = job->bodies[solid_body_id].particles.m.array() / job->bodies[solid_body_id].particles.v.array();
    pvec *= 1.0 / grains_rho;
    //pvec += Eigen::VectorXd::Ones(job->bodies[solid_body_id].p); //n = 1-phi
    n = job->bodies[solid_body_id].Phi * pvec;
    */
    n = Eigen::VectorXd::Ones(job->num_nodes);
    V = job->bodies[liquid_body_id].Phi * job->bodies[liquid_body_id].particles.v;

    /*
    //counted average of n on nodes
    Eigen::VectorXd one = job->bodies[solid_body_id].Phi * Eigen::VectorXd::Ones(job->bodies[solid_body_id].p);
    */

    //element volume average (empty element contribution)
    //Eigen::VectorXd filled(job->num_elements);
    //Eigen::VectorXd v_filled(job->num_nodes);
    Eigen::VectorXd v_total(job->num_nodes);
    //filled.setZero();
    //v_filled.setZero();
    v_total.setZero();
    //would be nice to have 'isfilled' on elements
    /*for (size_t p=0;p<job->bodies[solid_body_id].p;p++){
        for (size_t c=0;c<8;c++) {
            filled[job->bodies[solid_body_id].particles.corner_elements(p, c)] = 1;
        }
    }*/
    for (size_t e=0;e<job->num_elements;e++){
        for (size_t c = 0; c < 8; c++) {
            /*if (filled[e] != 0) {
                v_filled[job->elements.nodeID(e, c)] += 1;
            }*/
            v_total[job->elements.nodeID(e,c)] += job->hx*job->hy*job->hz/8;//1;
        }
    }

    /*for (size_t i=0;i<job->num_nodes;i++){
        if (job->bodies[solid_body_id].nodes.m[i] != 0){
            n[i] /= one[i];//job->bodies[solid_body_id].nodes.m[i];
            n[i] *= v_filled[i]/v_total[i];
            n[i] = 1.0 - n[i];
        } else {
            n[i] = 1.0;
        }
    }*/

    n = Eigen::VectorXd::Ones(job->num_nodes).array() - (job->bodies[solid_body_id].nodes.m.array()/(grains_rho * v_total).array());

    //gradnx = job->bodies[solid_body_id].gradPhiX * pvec;
    //gradny = job->bodies[solid_body_id].gradPhiY * pvec;
    //gradnz = job->bodies[solid_body_id].gradPhiZ * pvec;

    /*pvec = job->bodies[liquid_body_id].particles.v.array() * (job->bodies[liquid_body_id].particles.T.col(XX) + job->bodies[liquid_body_id].particles.T.col(YY) +
           job->bodies[liquid_body_id].particles.T.col(ZZ)).array();
    pvec *= -1.0 / 3.0;
    Eigen::VectorXd p_w(job->num_nodes);
    p_w = job->bodies[liquid_body_id].Phi * pvec;
    Eigen::VectorXd p_ws(job->bodies[solid_body_id].p);
    p_ws = job->bodies[solid_body_id].Phi.transpose() * p_w;
    p_ws = job->bodies[solid_body_id].particles.v.array() * p_ws.array();
    gradpx = job->bodies[solid_body_id].gradPhiX * p_ws;
    gradpy = job->bodies[solid_body_id].gradPhiY * p_ws;
    gradpz = job->bodies[solid_body_id].gradPhiZ * p_ws;*/

    Eigen::VectorXd sigma_ijw_n(job->num_nodes);
    Eigen::VectorXd sigma_ijw_p(job->bodies[solid_body_id].p);

    gradpx.setZero();
    gradpy.setZero();
    gradpz.setZero();

    /**************************************************************************/
    //find mass averaged stress on node
    pvec = job->bodies[liquid_body_id].particles.m.array() *
           job->bodies[liquid_body_id].particles.T.col(XX).array();
    sigma_ijw_n = job->bodies[liquid_body_id].Phi * pvec;
    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->bodies[liquid_body_id].nodes.m[i] > TOL) {
            sigma_ijw_n[i] /= job->bodies[liquid_body_id].nodes.m[i];
        } else {
            sigma_ijw_n[i] = 0;
        }
    }

    //interpolate stress onto solid
    sigma_ijw_p = job->bodies[solid_body_id].Phi.transpose() * sigma_ijw_n;
    //form volumetric stress vector
    sigma_ijw_p = job->bodies[solid_body_id].particles.m.array() * sigma_ijw_p.array() / grains_rho;
    //add to gradient
    //gradpx -= job->bodies[solid_body_id].gradPhiX * sigma_ijw_p;

    pvec = job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(XX).array();
    gradpx -= job->bodies[liquid_body_id].gradPhiX * pvec;

    /**************************************************************************/

    /**************************************************************************/
    //find mass averaged stress on node
    pvec = job->bodies[liquid_body_id].particles.m.array() *
           job->bodies[liquid_body_id].particles.T.col(XY).array();
    sigma_ijw_n = job->bodies[liquid_body_id].Phi * pvec;
    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->bodies[liquid_body_id].nodes.m[i] > TOL) {
            sigma_ijw_n[i] /= job->bodies[liquid_body_id].nodes.m[i];
        } else {
            sigma_ijw_n[i] = 0;
        }
    }
    //interpolate stress onto solid
    sigma_ijw_p = job->bodies[solid_body_id].Phi.transpose() * sigma_ijw_n;
    //form volumetric stress vector
    sigma_ijw_p = job->bodies[solid_body_id].particles.m.array() * sigma_ijw_p.array() / grains_rho;
    //add to gradient
    //gradpx -= job->bodies[solid_body_id].gradPhiY * sigma_ijw_p;
    //gradpy -= job->bodies[solid_body_id].gradPhiX * sigma_ijw_p;

    pvec = job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(XY).array();
    gradpx -= job->bodies[liquid_body_id].gradPhiY * pvec;
    gradpy -= job->bodies[liquid_body_id].gradPhiX * pvec;

    /**************************************************************************/

    /**************************************************************************/
    //find mass averaged stress on node
    pvec = job->bodies[liquid_body_id].particles.m.array() *
           job->bodies[liquid_body_id].particles.T.col(XZ).array();
    sigma_ijw_n = job->bodies[liquid_body_id].Phi * pvec;
    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->bodies[liquid_body_id].nodes.m[i] > TOL) {
            sigma_ijw_n[i] /= job->bodies[liquid_body_id].nodes.m[i];
        } else {
            sigma_ijw_n[i] = 0;
        }
    }
    //interpolate stress onto solid
    sigma_ijw_p = job->bodies[solid_body_id].Phi.transpose() * sigma_ijw_n;
    //form volumetric stress vector
    sigma_ijw_p = job->bodies[solid_body_id].particles.m.array() * sigma_ijw_p.array() / grains_rho;
    //add to gradient
    //gradpx -= job->bodies[solid_body_id].gradPhiZ * sigma_ijw_p;
    //gradpz -= job->bodies[solid_body_id].gradPhiX * sigma_ijw_p;

    pvec = job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(XZ).array();
    gradpx -= job->bodies[liquid_body_id].gradPhiZ * pvec;
    gradpz -= job->bodies[liquid_body_id].gradPhiX * pvec;

    /**************************************************************************/

    /**************************************************************************/
    //find mass averaged stress on node
    pvec = job->bodies[liquid_body_id].particles.m.array() *
           job->bodies[liquid_body_id].particles.T.col(YY).array();
    sigma_ijw_n = job->bodies[liquid_body_id].Phi * pvec;
    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->bodies[liquid_body_id].nodes.m[i] > TOL) {
            sigma_ijw_n[i] /= job->bodies[liquid_body_id].nodes.m[i];
        } else {
            sigma_ijw_n[i] = 0;
        }
    }
    //interpolate stress onto solid
    sigma_ijw_p = job->bodies[solid_body_id].Phi.transpose() * sigma_ijw_n;
    //form volumetric stress vector
    sigma_ijw_p = job->bodies[solid_body_id].particles.m.array() * sigma_ijw_p.array() / grains_rho;
    //add to gradient
    //gradpy -= job->bodies[solid_body_id].gradPhiY * sigma_ijw_p;

    pvec = job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(YY).array();
    gradpy -= job->bodies[liquid_body_id].gradPhiY * pvec;

    /**************************************************************************/

    /**************************************************************************/
    //find mass averaged stress on node
    pvec = job->bodies[liquid_body_id].particles.m.array() *
           job->bodies[liquid_body_id].particles.T.col(YZ).array();
    sigma_ijw_n = job->bodies[liquid_body_id].Phi * pvec;
    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->bodies[liquid_body_id].nodes.m[i] > TOL) {
            sigma_ijw_n[i] /= job->bodies[liquid_body_id].nodes.m[i];
        } else {
            sigma_ijw_n[i] = 0;
        }
    }
    //interpolate stress onto solid
    sigma_ijw_p = job->bodies[solid_body_id].Phi.transpose() * sigma_ijw_n;
    //form volumetric stress vector
    sigma_ijw_p = job->bodies[solid_body_id].particles.m.array() * sigma_ijw_p.array() / grains_rho;
    //add to gradient
    //gradpy -= job->bodies[solid_body_id].gradPhiZ * sigma_ijw_p;
    //gradpz -= job->bodies[solid_body_id].gradPhiY * sigma_ijw_p;

    pvec = job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(YZ).array();
    gradpy -= job->bodies[liquid_body_id].gradPhiZ * pvec;
    gradpz -= job->bodies[liquid_body_id].gradPhiY * pvec;

    /**************************************************************************/

    /**************************************************************************/
    //find mass averaged stress on node
    pvec = job->bodies[liquid_body_id].particles.m.array() *
           job->bodies[liquid_body_id].particles.T.col(ZZ).array();
    sigma_ijw_n = job->bodies[liquid_body_id].Phi * pvec;
    for (size_t i = 0; i < job->num_nodes; i++) {
        if (job->bodies[liquid_body_id].nodes.m[i] > TOL) {
            sigma_ijw_n[i] /= job->bodies[liquid_body_id].nodes.m[i];
        } else {
            sigma_ijw_n[i] = 0;
        }
    }
    //interpolate stress onto solid
    sigma_ijw_p = job->bodies[solid_body_id].Phi.transpose() * sigma_ijw_n;
    //form volumetric stress vector
    sigma_ijw_p = job->bodies[solid_body_id].particles.m.array() * sigma_ijw_p.array() / grains_rho;
    //add to gradient
    //gradpz -= job->bodies[solid_body_id].gradPhiZ * sigma_ijw_p;

    pvec = job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(ZZ).array();
    gradpz -= job->bodies[liquid_body_id].gradPhiZ * pvec;

    /**************************************************************************/



    Eigen::VectorXd sgradnx(job->num_nodes);
    Eigen::VectorXd sgradny(job->num_nodes);
    Eigen::VectorXd sgradnz(job->num_nodes);
    sgradnx.setZero();
    sgradny.setZero();
    sgradnz.setZero();

    //calculate sigma grad n term
    pvec = job->bodies[liquid_body_id].gradPhiX.transpose() * n;
    pvec = pvec.array() *
           job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(XX).array();
    sgradnx += job->bodies[liquid_body_id].Phi * pvec;


    pvec = job->bodies[liquid_body_id].gradPhiY.transpose() * n;
    pvec = pvec.array() *
           job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(XY).array();
    sgradnx += job->bodies[liquid_body_id].Phi * pvec;

    pvec = job->bodies[liquid_body_id].gradPhiZ.transpose() * n;
    pvec = pvec.array() *
           job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(XZ).array();
    sgradnx += job->bodies[liquid_body_id].Phi * pvec;

    //sgradnx = sgradnx.array() * v_total.array() / (job->bodies[liquid_body_id].Phi * job->bodies[liquid_body_id].particles.v).array();


    pvec = job->bodies[liquid_body_id].gradPhiX.transpose() * n;
    pvec = pvec.array() *
           job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(XY).array();
    sgradny += job->bodies[liquid_body_id].Phi * pvec;

    pvec = job->bodies[liquid_body_id].gradPhiY.transpose() * n;
    pvec = pvec.array() *
           job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(YY).array();
    sgradny += job->bodies[liquid_body_id].Phi * pvec;

    pvec = job->bodies[liquid_body_id].gradPhiZ.transpose() * n;
    pvec = pvec.array() *
           job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(YZ).array();
    sgradny += job->bodies[liquid_body_id].Phi * pvec;

    //sgradny = sgradny.array() * v_total.array() / (job->bodies[liquid_body_id].Phi * job->bodies[liquid_body_id].particles.v).array();


    pvec = job->bodies[liquid_body_id].gradPhiX.transpose() * n;
    pvec = pvec.array() *
           job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(XZ).array();
    sgradnz += job->bodies[liquid_body_id].Phi * pvec;

    pvec = job->bodies[liquid_body_id].gradPhiY.transpose() * n;
    pvec = pvec.array() *
           job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(YZ).array();
    sgradnz += job->bodies[liquid_body_id].Phi * pvec;

    pvec = job->bodies[liquid_body_id].gradPhiZ.transpose() * n;
    pvec = pvec.array() *
           job->bodies[liquid_body_id].particles.v.array() *
           job->bodies[liquid_body_id].particles.T.col(ZZ).array();
    sgradnz += job->bodies[liquid_body_id].Phi * pvec;

    //sgradnz = sgradnz.array() * v_total.array() / (job->bodies[liquid_body_id].Phi * job->bodies[liquid_body_id].particles.v).array();

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
            //damping enforcement
            double C;

            //k = 0.0055 [calliope.dem.uniud.it/CLASS/DES-IND-PLA/CR-FlowThroughBed.pdf]
            double k_eff = k*n[i]*n[i]*n[i]*grains_diam*grains_diam/((1-n[i])*(1-n[i]));

            C = n[i]*V[i]*mu_w/k_eff; //n[i]*n[i]*V[i]*mu_w/k;
            if (!std::isfinite(k_eff)){
                C = 0;
            }

            if ((C*job->dt*(1.0/m1 + 1.0/m2)) > 1){
                //std::cout << "timestep too large for slurry model" << std::endl;
                C = 1.0/(job->dt*(1.0/m1 + 1.0/m2));
            }
            //fsfi = (mv1i / m1 - mv2i / m2) * n[i] * nV[i] * mu_w / k;
            fsfi = (mv1i / m1 - mv2i / m2)*C;
            //fsfi[0] -= sgradnx[i];
            //fsfi[1] -= sgradny[i];
            //fsfi[2] -= sgradnz[i];

            //set contact forces
            if (job->t >= t_start) {
                job->bodies[solid_body_id].nodes.contact_fx[i] += -fsfi[0] + (1-n[i])/n[i] * (gradpx[i] - sgradnx[i]);
                job->bodies[solid_body_id].nodes.contact_fy[i] += -fsfi[1] + (1-n[i])/n[i] * (gradpy[i] - sgradny[i]);
                job->bodies[solid_body_id].nodes.contact_fz[i] += -fsfi[2] + (1-n[i])/n[i] * (gradpz[i] - sgradnz[i]);
            }

            job->bodies[liquid_body_id].nodes.contact_fx[i] += fsfi[0] - sgradnx[i];
            job->bodies[liquid_body_id].nodes.contact_fy[i] += fsfi[1] - sgradny[i];
            job->bodies[liquid_body_id].nodes.contact_fz[i] += fsfi[2] - sgradnz[i];

            //std::cout << gradpx[i] << " , " << job->bodies[liquid_body_id].nodes.contact_fx[i] << std::endl;

            //determine strainrate variable on nodes with new nodal velocity
            /*mv1i -= job->dt * fsfi;
            mv2i += job->dt * fsfi;
            nvecx[i] = (1 - n[i]) * mv1i[0] / m1 + n[i] * mv2i[0] / m2;
            nvecy[i] = (1 - n[i]) * mv1i[1] / m1 + n[i] * mv2i[1] / m2;
            nvecz[i] = (1 - n[i]) * mv1i[2] / m1 + n[i] * mv2i[2] / m2;*/
        } /*else if(job->bodies[liquid_body_id].nodes.m[i] != 0){
            double m2 = job->bodies[liquid_body_id].nodes.m[i];
            Eigen::Vector3d mv2i;
            if (job->use_implicit == 1) {
                mv2i << job->bodies[liquid_body_id].nodes.contact_mx_t[i],
                        job->bodies[liquid_body_id].nodes.contact_my_t[i],
                        job->bodies[liquid_body_id].nodes.contact_mz_t[i];
            } else {
                mv2i << job->bodies[liquid_body_id].nodes.contact_mx_t[i] +
                        job->dt * job->bodies[liquid_body_id].nodes.contact_fx[i],
                        job->bodies[liquid_body_id].nodes.contact_my_t[i] +
                        job->dt * job->bodies[liquid_body_id].nodes.contact_fy[i],
                        job->bodies[liquid_body_id].nodes.contact_mz_t[i] +
                        job->dt * job->bodies[liquid_body_id].nodes.contact_fz[i];
            }
            nvecx[i] = n[i] * mv2i[0] / m2;
            nvecy[i] = n[i] * mv2i[1] / m2;
            nvecz[i] = n[i] * mv2i[2] / m2;
        }*/
    }

    //zero what needs to be zeroed
    job->addBoundaryConditions();
    if (job->use_3d == 0){
        job->bodies[solid_body_id].nodes.contact_fz.setZero();
        job->bodies[liquid_body_id].nodes.contact_fz.setZero();
    }

    for (size_t i = 0; i < job->num_nodes; i++) {
        double m1 = job->bodies[solid_body_id].nodes.m[i];
        double m2 = job->bodies[liquid_body_id].nodes.m[i];

        if (m1 != 0 && m2 != 0){
            Eigen::Vector3d mv1i;
            Eigen::Vector3d mv2i;
            Eigen::Vector3d vCMi;
            if (job->use_implicit == 1) {
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

            nvecx[i] = (1 - n[i]) * mv1i[0] / m1 + n[i] * mv2i[0] / m2;
            nvecy[i] = (1 - n[i]) * mv1i[1] / m1 + n[i] * mv2i[1] / m2;
            nvecz[i] = (1 - n[i]) * mv1i[2] / m1 + n[i] * mv2i[2] / m2;
        } else if(m2 != 0){
            Eigen::Vector3d mv2i;
            if (job->use_implicit == 1) {
                mv2i << job->bodies[liquid_body_id].nodes.contact_mx_t[i],
                        job->bodies[liquid_body_id].nodes.contact_my_t[i],
                        job->bodies[liquid_body_id].nodes.contact_mz_t[i];
            } else {
                mv2i << job->bodies[liquid_body_id].nodes.contact_mx_t[i] +
                        job->dt * job->bodies[liquid_body_id].nodes.contact_fx[i],
                        job->bodies[liquid_body_id].nodes.contact_my_t[i] +
                        job->dt * job->bodies[liquid_body_id].nodes.contact_fy[i],
                        job->bodies[liquid_body_id].nodes.contact_mz_t[i] +
                        job->dt * job->bodies[liquid_body_id].nodes.contact_fz[i];
            }
            nvecx[i] = n[i] * mv2i[0] / m2;
            nvecy[i] = n[i] * mv2i[1] / m2;
            nvecz[i] = n[i] * mv2i[2] / m2;
        }
    }

    //determine strainrate on fluid
    Eigen::VectorXd pvec2(job->bodies[liquid_body_id].p);
    job->bodies[liquid_body_id].particles.state.col(10) << (job->bodies[liquid_body_id].gradPhiX.transpose() * nvecx +
                                                           job->bodies[liquid_body_id].gradPhiY.transpose() * nvecy +
                                                           job->bodies[liquid_body_id].gradPhiZ.transpose() * nvecz);
    pvec2 = job->bodies[liquid_body_id].Phi.transpose() * n;

    job->bodies[liquid_body_id].particles.updateElementIDs(job);
    for (size_t p = 0; p < job->bodies[liquid_body_id].p; p++) {
        if (pvec2[p] != 0 && pvec2[p] < 1 && job->bodies[liquid_body_id].particles.active[p] == 1) {
            job->bodies[liquid_body_id].particles.state(p, 9) = 1;
            job->bodies[liquid_body_id].particles.state(p, 10) /= pvec2[p];

            //check that the element has internal solid mass
            /*
            for (size_t n=0; n<8; n++){
                if (job->bodies[solid_body_id].nodes.m[job->elements.nodeID(job->bodies[liquid_body_id].particles.elementIDs[p],n)] <= TOL){
                    job->bodies[liquid_body_id].particles.state(p, 9) = 0;
                    break;
                }
            }*/
        } else {
            job->bodies[liquid_body_id].particles.state(p, 9) = 0;
        }

        if (!std::isfinite(job->bodies[liquid_body_id].particles.state(p, 10))) {
            std::cout << "SHOOT HER! " << job->bodies[liquid_body_id].particles.state(p, 10) << std::endl;
            exit(0);
        }
    }

    //std::cout << "[" << n.maxCoeff() << "," << n.minCoeff() << "]" << std::endl;
    //std::cout << nvecy[14*5 + 2] << "," << nvecy[15*5 + 2] << "," << nvecy[16*5 + 2] << "," << nvecy[17*5 + 2] << std::endl;
    //std::cout << gradpy.maxCoeff() << " , " << gradpy.minCoeff() << std::endl;

    return;
}
/*----------------------------------------------------------------------------*/
//
// Created by aaron on 5/23/18.
// mixture_slurry.cpp
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

void SlurryMixture::init(Job* job) {
    //check that contact properties are set
    if (fp64_props.size() < 4 || (str_props.size() < 2 && int_props.size() < 2)) {
        //need mixture properties and bodies
        std::cout << fp64_props.size() << ", " << int_props.size() << ", "
                  << str_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 4 properties defined ({grain_rho, grain_diam, eta_0, fluid_rho},{body_1,body_2}).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //set coeff of friction
        grains_rho = fp64_props[0];
        grains_d = fp64_props[1];
        eta_0 = fp64_props[2];
        fluid_rho = fp64_props[3];

        if (int_props.size() == 1) {
            spec_override = int_props[0];
        } else if (int_props.size() == 3) {
            spec_override = int_props[2];
        }

        //set body ids by name
        if (str_props.size() >= 2) {
            for (int i = 0; i < bodyIDs.size(); i++) {
                for (int b = 0; b < job->bodies.size(); b++) {
                    if (str_props[i].compare(job->bodies[b]->name) == 0) {
                        bodyIDs[i] = b;
                        break;
                    }
                }
            }
        }

        // or set body ids by int
        for (int i = 0; i < bodyIDs.size(); i++) {
            if (bodyIDs[i] < 0) {
                if (int_props.size() == 2) {
                    bodyIDs = int_props;
                } else {
                    std::cout << fp64_props.size() << ", " << int_props.size() << ", "
                              << str_props.size() << "\n";
                    fprintf(stderr,
                            "%s:%s: Need at least 4 properties defined ({grain_rho, grain_diam, eta_0, fluid_rho},{solid,fluid}).\n",
                            __FILE__, __func__);
                    exit(0);
                }
                break;
            }
        }

        solid_body_id = bodyIDs[0];
        fluid_body_id = bodyIDs[1];

        // Parse Input Strings for Simulation Flags
        for (int i = 0; i < str_props.size(); i++) {
            std::vector<std::string> options = {"USE_CARMAN_KOZENY"};
            switch (Parser::findStringID(options, str_props[i])) {
                case 0:
                    //USE_CARMAN_KOZENY
                    USE_CARMAN_KOZENY = true;
                    std::cout << "MixtureSlurry using Carman-Kozeny drag formula." << std::endl;
                    break;
                default:
                    //do nothing
                    break;
            }
        }

        V.resize(job->bodies[fluid_body_id]->nodes->x.size());
        n.resize(job->bodies[solid_body_id]->nodes->x.size());

        divT = KinematicVectorArray(job->bodies[solid_body_id]->nodes->x.size(),job->JOB_TYPE);

        printf("Contact properties (grains_rho = %g, grains_d = %g, eta_0 = %g, fluid_rho = %g, {solid: %i, fluid: %i}).\n",
               grains_rho, grains_d, eta_0, fluid_rho, solid_body_id, fluid_body_id);
    }

    std::cout << "Contact Initialized: [" << id << "]." << std::endl;

    return;
}

void SlurryMixture::generateRules(Job* job){
    //set porosity from solid mass
    if (job->JOB_TYPE == job->JOB_AXISYM){
        Eigen::MatrixXd pval = Eigen::VectorXd(job->bodies[solid_body_id]->points->x.size());
        Eigen::MatrixXd nval = Eigen::VectorXd(job->bodies[solid_body_id]->nodes->x.size());

        //adjust mass to 2D integral of density
        for (int i=0;i<job->bodies[solid_body_id]->points->x.size();i++){
            pval(i) = job->bodies[solid_body_id]->points->m(i)/job->bodies[solid_body_id]->points->x(i,0); // A*rho = v/r * m/v
        }
        nval = job->bodies[solid_body_id]->S * pval;

        for (int i = 0; i < n.rows(); i++) {
            //n = 1 - phi
            n(i) = 1 - (nval(i) / (job->grid->nodeVolume(job, i) * grains_rho));

            /*if (n(i) < 0.2) {
                n(i) = 0.2; //keep packing from overestimates...
            }*/
            if (n(i) < 1e-2){
                n(i) = 1e-2; //keep packing from overestimates...
            }
        }

        //approximate volume as liquid volume
        pval = Eigen::VectorXd(job->bodies[fluid_body_id]->points->x.size());
        //adjust volume to 2D integral of area
        for (int i=0;i<job->bodies[fluid_body_id]->points->x.size();i++){
            pval(i) = job->bodies[fluid_body_id]->points->v(i)/job->bodies[fluid_body_id]->points->x(i,0); // A = v/r
        }
        V = job->bodies[fluid_body_id]->S * pval;

    } else {
        for (int i = 0; i < n.rows(); i++) {
            //n = 1 - phi
            n(i) = 1 - (job->bodies[solid_body_id]->nodes->m(i) / (job->grid->nodeVolume(job, i) * grains_rho));

            /*if (n(i) < 0.2) {
                n(i) = 0.2; //keep packing from overestimates...
            }*/
            if (n(i) < 1e-2){
                n(i) = 1e-2; //keep packing from overestimates...
            }
        }

        //approximate V as integrated liquid volume
        V = job->bodies[fluid_body_id]->S * job->bodies[fluid_body_id]->points->v;
    }
    return;
}

void SlurryMixture::applyRules(Job* job, int SPEC){
    //check for override
    if (spec_override != -1){
        SPEC = spec_override;
    }

    double m1, m2;
    KinematicVector fsfi = KinematicVector(job->JOB_TYPE);
    KinematicVector mv1i = KinematicVector(job->JOB_TYPE);
    KinematicVector mv2i = KinematicVector(job->JOB_TYPE);
    KinematicVector vCMi = KinematicVector(job->JOB_TYPE);

    double C, k_eff, Re;

    //calculate gradient of fluid pressure
    int len = job->bodies[fluid_body_id]->points->x.size();
    Eigen::VectorXd fluid_pressure = Eigen::VectorXd(len);
    for (int i=0;i<len;i++){
        fluid_pressure(i) = -job->bodies[fluid_body_id]->points->T(i).trace()/3.0 * job->bodies[fluid_body_id]->points->v(i);
    }

    //job->bodies[liquid_body_id].bodyCalcNodalGradient(job,divT,fluid_pressure,Body::SET);
    divT = job->bodies[fluid_body_id]->gradS * fluid_pressure;

    //calculate nodal momentum exchange
    for (int i = 0; i < n.rows(); i++) {
        //test every node for contact
        if (job->bodies[solid_body_id]->nodes->m[i] > 0 && job->bodies[fluid_body_id]->nodes->m[i] > 0) {
            //distribute solid stess
            job->bodies[solid_body_id]->nodes->f(i) += (1-n(i))*divT(i);
            job->bodies[fluid_body_id]->nodes->f(i) -= (1-n(i))*divT(i);

            //determine 'center of mass' velocity
            m1 = job->bodies[solid_body_id]->nodes->m[i];
            m2 = job->bodies[fluid_body_id]->nodes->m[i];
            if (SPEC == Contact::IMPLICIT) {
                mv1i = (job->bodies[solid_body_id]->nodes->mx_t(i) + job->dt * job->bodies[solid_body_id]->nodes->f(i));
                mv2i = (job->bodies[fluid_body_id]->nodes->mx_t(i) + job->dt * job->bodies[fluid_body_id]->nodes->f(i));
            } else if (SPEC == Contact::EXPLICIT){
                mv1i = job->bodies[solid_body_id]->nodes->mx_t(i);
                mv2i = job->bodies[fluid_body_id]->nodes->mx_t(i);
            } else {
                std::cerr << "ERROR: Unknown SPEC given to SlurryMixture: " << SPEC << "!" << std::endl;
                return;
            }
            vCMi  = (mv1i + mv2i) / (m1 + m2);

            //permeability
            Re = (mv1i/m1 - mv2i/m2).norm() * n(i) * fluid_rho * grains_d / eta_0;
            if (USE_CARMAN_KOZENY){
                C = V(i) * 180.0 * (1 - n(i)) * (1 - n(i)) * eta_0 / (n(i) * grains_d * grains_d);
            } else if (Re == 0){
                //Beetstra
                C = V(i) * 18.0 * (1 - n(i)) * eta_0 / (grains_d * grains_d) *
                    (10.0 * (1 - n(i)) / n(i) + n(i) * n(i) * n(i) * (1.0 + 1.5 * std::sqrt(1 - n(i))));
            } else {
                //Beetstra
                C = V(i) * 18.0 * (1 - n(i)) * eta_0 / (grains_d * grains_d) *
                    (10.0 * (1 - n(i)) / n(i) + n(i) * n(i) * n(i) * (1.0 + 1.5 * std::sqrt(1 - n(i))) +
                     0.413 * Re / (24.0 * n(i)) *
                     (1 / n(i) + 3 * n(i) * (1 - n(i)) + 8.4 * std::pow(Re, -0.343)) /
                     (1 + std::pow(10.0, 3 * (1 - n(i))) * std::pow(Re, -0.5 * (1 + 4 * (1 - n(i))))));
            }

            //fsfi = (mv1i / m1 - mv2i / m2)*C;
            fsfi = C/(1 + job->dt*C*(1/m1 + 1/m2)) * (mv1i/m1 - mv2i/m2);

            job->bodies[solid_body_id]->nodes->f(i) -= fsfi;
            job->bodies[fluid_body_id]->nodes->f(i) += fsfi;
        }
    }

    return;
}

/*----------------------------------------------------------------------------*/

void SlurryMixture::writeFrame(Job* job, Serializer* serializer){
    serializer->writeScalarArray(n,("n_" + job->bodies[solid_body_id]->name + "_" + job->bodies[fluid_body_id]->name));
    serializer->writeVectorArray(divT,"divT");
    return;
}

std::string SlurryMixture::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

int SlurryMixture::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

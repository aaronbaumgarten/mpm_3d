//
// Created by aaron on 5/16/18.
// explicit_usl.cpp
//

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "solvers.hpp"
#include "objects/bodies/bodies.hpp"
#include "objects/materials/materials.hpp"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

/*----------------------------------------------------------------------------*/
//
void ExplicitMixtureSolver::init(Job* job){

    //initialize parent class
    ExplicitUSL::init(job);

    //identify bodies
    if (str_props.size() < 3 || fp64_props.size() < 2){
        std::cerr << str_props.size() << " ?= 3\n";
        std::cerr << fp64_props.size() << " ?= 2\n";
        fprintf(stderr, "%s:%s:", __FILE__, __func__);

        std::cerr << "ExplicitMixtureSolver needs at least 3 bodies defined by name.\n";
        std::cerr << "    granular_body:\n";
        std::cerr << "    fluid_body:\n";
        std::cerr << "    solid_body:\n";
        std::cerr << "ExplicitMixtureSolver needs at least 2 material properties defined.\n";
        std::cerr << "    grains_d:\n";
        std::cerr << "    rhof_0:\n";
        exit(0);
    } else {

        //assign material property
        grains_d    = fp64_props[0];
        rhof_0      = fp64_props[1];

        //initialization boolean
        bool initilization_fail = false;

        //set bodyIDs by name
        std::vector<int> bodyIDs = {-1, -1, -1};
        for (int i = 0; i < bodyIDs.size(); i++) {
            for (int b = 0; b < job->bodies.size(); b++) {
                if (str_props[i].compare(job->bodies[b]->name) == 0) {
                    bodyIDs[i] = b;
                    break;
                }
            }
            if (bodyIDs[i] < 0){
                std::cerr << "Uh oh! [" << i << "]" << std::endl;
                initilization_fail = true;
                break;
            }
        }

        if (initilization_fail){
            std::cerr << "ExplicitMixtureSolver needs at least 3 bodies defined by name.\n";
            std::cerr << "    granular_body: " << str_props[0] << "\n";
            std::cerr << "    fluid_body: " << str_props[1] << "\n";
            std::cerr << "    solid_body: " << str_props[2] << "\n";

            std::cerr << "Names: ";
            for (int b = 0; b<job->bodies.size(); b++){
                std::cerr << job->bodies[b]->name << ", ";
            }
            std::cerr << "\n";

            exit(0);
        }

        //assign pointers
        granular_body = job->bodies[bodyIDs[0]].get();
        fluid_body    = job->bodies[bodyIDs[1]].get();
        solid_body    = job->bodies[bodyIDs[2]].get();

        //tell console about findings!
        std::cout << "ExplicitMixtureSolver has identified 3 bodies!\n";
        std::cout << "    granular_body: " << granular_body->name << "\n";
        std::cout << "    fluid_body: "    << fluid_body->name << "\n";
        std::cout << "    solid_body: "    << solid_body->name << "\n";
        std::cout << "ExplicitMixtureSolver has 2 material properties!\n";
        std::cout << "    grains_d: "       << grains_d << "\n";
        std::cout << "    rhof_0: "         << rhof_0 << "\n";

        //attempt to identify granular material model
        if (granular_body->material->object_name.compare("CompressibleBreakageMechanicsSand") == 0){
            //success!

            //attempt high risk of pointer to granular material model
            compressible_breakage_mechanics_sand_model = dynamic_cast<CompressibleBreakageMechanicsSand*>(granular_body->material.get());

            //assume success and deal with the consequences
            std::cout << "ExplicitMixtureSolver has identified granular material model!\n";
            std::cout << "   " << compressible_breakage_mechanics_sand_model->object_name << " =? CompressibleBreakageMechanicsSand\n";
        } else {
            std::cerr << "ExplicitMixtureSolver has failed to identify granular material model!\n";
            std::cerr << "   " << granular_body->material->object_name << " =? CompressibleBreakageMechanicsSand\n";
            exit(0);
        }

        //attempt to identify fluid material model
        if (fluid_body->material->object_name.compare("BarotropicViscousFluid") == 0){
            //success!
            fluid_model = 0;

            //attempt high risk of pointer to granular material model
            barotropic_viscous_fluid_model = dynamic_cast<BarotropicViscousFluid*>(fluid_body->material.get());

            //assume success and deal with the consequences
            std::cout << "ExplicitMixtureSolver has identified fluid material model!\n";
            std::cout << "   " << barotropic_viscous_fluid_model->object_name << " =? BarotropicViscousFluid\n";
        } else if (fluid_body->material->object_name.compare("TillotsonEOSFluid") == 0){
            //success!
            fluid_model = 1;

            //attempt high risk of pointer to granular material model
            tillotson_eos_fluid_model = dynamic_cast<TillotsonEOSFluid*>(fluid_body->material.get());

            //assume success and deal with the consequences
            std::cout << "ExplicitMixtureSolver has identified fluid material model!\n";
            std::cout << "   " << tillotson_eos_fluid_model->object_name << " =? TillotsonEOSFluid\n";
        } else {
            std::cerr << "ExplicitMixtureSolver has failed to identify granular material model!\n";
            std::cerr << "   " << fluid_body->material->object_name << " =? BarotropicViscousFluid\n";
            std::cerr << "   " << fluid_body->material->object_name << " =? TillotsonEOSFluid\n";
            exit(0);
        }


        //check for string flags
        for (int s=0; s<str_props.size(); s++){
            if (str_props[s].compare("ISOTHERMAL") == 0){
                //isothermal simulation
                is_adiabatic = false;
            } else if (str_props[s].compare("INCOMPRESSIBLE") == 0){
                //incompressible model
                is_compressible = false;
            } else if (str_props[s].compare("REFLECTED_BOUNDARY") == 0){
                //use reflected boundary for contact computation
                use_reflected_boundary = true;
            }
        }


        if (is_compressible){
            std::cout << "  Using Compressible Model.\n";
        } else {
            std::cout << "  Using Incompressible Model.\n";
        }

        if (is_adiabatic){
            std::cout << "  Using Adiabatic Model.\n";
        } else {
            std::cout << "  Using Isothermal Model.\n";
        }


        //initialize solver data structures
        //need to resize after simulation start
        n = Eigen::VectorXd(0);
        rho_f = Eigen::VectorXd(0);
        drhos_dt = Eigen::VectorXd(0);
        drhof_dt = Eigen::VectorXd(0);

        //initialization complete!
        std::cout << "ExplicitMixtureSolver Initialized!" << std::endl;

    }

    return;
}

/*----------------------------------------------------------------------------*/
//
void ExplicitMixtureSolver::step(Job* job){
    //create map
    createMappings(job);

    //map particles to grid
    mapPointsToNodes(job);

    //add interaction forces
    addInteractionForces(job);

    //enforce boundary conditions
    generateBoundaryConditions(job);
    addBoundaryConditions(job);

    //move grid
    moveGrid(job);

    //move particles
    movePoints(job);

    //calculate strainrate
    calculateStrainRate(job);

    //update density
    updateDensity(job);

    //add body forces
    job->driver->generateGravity(job);
    job->driver->applyGravity(job);

    //update stress
    updateStress(job);

    return;
}


void ExplicitMixtureSolver::addInteractionForces(Job *job) {

    //-----------------------------------------------
    // 0 -- Initialize Uninitialized Data Structures
    //-----------------------------------------------

    // Reflected Boundary Lengths
    if (!reflected_boundary_initialized) {

        // Assume Box Domain w/ 0,0,0 as Lower-Left Boundary
        Lx = KinematicVector(job->JOB_TYPE);
        Lx.setZero();
        for (int i = 0; i < granular_body->nodes->x.size(); i++) {
            for (int pos = 0; pos < granular_body->nodes->x.DIM; pos++) {
                if (granular_body->nodes->x(i, pos) > Lx(pos)) {
                    Lx(pos) = granular_body->nodes->x(i, pos);
                }
            }
        }

        // Change Initialization Flag
        reflected_boundary_initialized = true;
    }

    // Porosity, Density, and Dr/Dt
    if (n.size() != fluid_body->nodes->x.size()){
        n.resize(fluid_body->nodes->x.size());
        n.setOnes();
    }
    if (rho_f.size() != fluid_body->points->x.size()){
        rho_f.resize(fluid_body->points->x.size());
        rho_f.setConstant(rhof_0);
    }
    if (drhos_dt.size() != granular_body->points->x.size()){
        drhos_dt.resize(granular_body->points->x.size());
        drhos_dt.setZero();
    }
    if (drhof_dt.size() != fluid_body->points->x.size()){
        drhof_dt.resize(fluid_body->points->x.size());
        drhof_dt.setZero();
    }



    //----------------------------------
    // 1 -- Compute Interaction Normals
    //----------------------------------

    // Create Normals
    KinematicVectorArray contact_normal = KinematicVectorArray(fluid_body->nodes->x.size(), job->JOB_TYPE);

    // Compute Normals
    if (job->JOB_TYPE == job->JOB_AXISYM){
        Eigen::VectorXd pval = Eigen::VectorXd(solid_body->points->x.size());
        for (int i=0; i<pval.rows(); i++){
            pval(i) = solid_body->points->m(i)/solid_body->points->x(i,0);
        }
        //calculate normal using 2D integral of density
        contact_normal = solid_body->gradS * pval;
    } else {
        contact_normal = solid_body->gradS * solid_body->points->m;
    }

    for (int i=0;i<contact_normal.size();i++){
        if (solid_body->nodes->m[i] > 0) {
            if (use_reflected_boundary) {
                for (int pos = 0; pos < job->DIM; pos++) {
                    if (solid_body->nodes->x(i, pos) == 0 ||
                        solid_body->nodes->x(i, pos) == Lx(pos)) {
                        contact_normal(i, pos) = 0;
                    }
                }
            }
            //normalize
            contact_normal(i) /= contact_normal(i).norm();
        }
    }



    //----------------------------------------------------
    // 2 -- Compute Porosity and Fluid Integration Volume
    //----------------------------------------------------

    // fluid integration volumes
    Eigen::VectorXd V = Eigen::VectorXd(fluid_body->nodes->x.size());

    // set porosity from solid volume fraction
    Eigen::VectorXd phiS = compressible_breakage_mechanics_sand_model->phiS;

    // intermediate vectors
    Eigen::MatrixXd pval = Eigen::VectorXd(granular_body->points->x.size());
    Eigen::MatrixXd nval = Eigen::VectorXd(granular_body->nodes->x.size());

    // compute nodal porosity
    if (job->JOB_TYPE == job->JOB_AXISYM){
        // adjust point volume for 2D integration
        for (int i=0; i<granular_body->points->x.size(); i++){
            pval(i) = phiS(i) * granular_body->points->v(i) / granular_body->points->x(i,0); // A = v/r
        }
        nval = granular_body->S * pval;

        // divide integrated solid volume by 2D nodal volume
        for (int i = 0; i < n.rows(); i++) {
            // n = 1 - phiS
            n(i) = 1.0 - (nval(i) / job->grid->nodeVolume(job, i));

            // integration errors need to be handled
            if (n(i) < 1e-2){
                n(i) = 1e-2; //keep packing from overestimates...
            }
        }

        // approximate volume as liquid volume
        pval = Eigen::VectorXd(fluid_body->points->x.size());

        // adjust volume to 2D integral of area
        for (int i=0; i<fluid_body->points->x.size(); i++){
            pval(i) = fluid_body->points->v(i)/fluid_body->points->x(i,0); // A = v/r
        }
        V = fluid_body->S * pval;

    } else {
        for (int i=0; i<granular_body->points->x.size(); i++){
            pval(i) = phiS(i) * granular_body->points->v(i);
        }
        nval = granular_body->S * pval;

        // divide integrated solid volume by nodal volume
        for (int i = 0; i < n.rows(); i++) {
            // n = 1 - phiS
            n(i) = 1.0 - (nval(i) / job->grid->nodeVolume(job, i));

            // integration errors need to be handled
            if (n(i) < 1e-2){
                n(i) = 1e-2; //keep packing from overestimates...
            }
        }

        //approximate V as integrated liquid volume
        V = fluid_body->S * fluid_body->points->v;
    }



    //---------------------------------------------------------
    // 3 -- Compute Pressure Contribution to Interaction Force
    //---------------------------------------------------------

    //calculate gradient of fluid pressure
    int len = fluid_body->points->x.size();
    Eigen::VectorXd p_f = Eigen::VectorXd(len);
    for (int i=0;i<len;i++){
        p_f(i) = -fluid_body->points->T(i).trace()/3.0 * fluid_body->points->v(i);
    }
    KinematicVectorArray gradP = fluid_body->gradS * p_f;

    //calculate force contribution
    // SUBTRACT from solid
    // ADD to fluid
    KinematicVectorArray fb_i = KinematicVectorArray(fluid_body->nodes->x.size(), job->JOB_TYPE);
    for (int i=0; i<fluid_body->nodes->x.size(); i++){
        fb_i[i] = -(1.0 - n(i)) * gradP(i);

        granular_body->nodes->f(i)  -= fb_i[i];
        fluid_body->nodes->f(i)     += fb_i[i];
    }



    //-----------------------------------------------------
    // 4 -- Compute Drag Contribution to Interaction Force
    //-----------------------------------------------------

    double m1, m2, eta_0, C;
    KinematicVector mv1i, mv2i, vCMi;
    KinematicVectorArray fd_i = KinematicVectorArray(fluid_body->nodes->x.size(), job->JOB_TYPE);

    //get eta
    if (fluid_model == 0){
        //BarotropicViscousFluid
        eta_0 = barotropic_viscous_fluid_model->eta;
    } else if (fluid_model == 1){
        //TillotsonEOSFluid
        eta_0 = tillotson_eos_fluid_model->eta0;
    } else {
        std::cerr << "Something has gone terribly wrong in ExplicitMixtureSolver!\n";
        std::cerr << "    fluid_model = " << fluid_model << std::endl;
    }

    //calculate nodal momentum exchange
    for (int i = 0; i < n.rows(); i++) {
        //test every node for contact
        if (granular_body->nodes->m[i] > 0 && fluid_body->nodes->m[i] > 0) {

            //get nodal masses
            m1 = granular_body->nodes->m[i];
            m2 = fluid_body->nodes->m[i];

            //get nodal momenta
            mv1i = granular_body->nodes->mx_t(i) + job->dt * granular_body->nodes->f(i);
            mv2i = fluid_body->nodes->mx_t(i)    + job->dt * fluid_body->nodes->f(i);

            //determine 'center of mass' velocity
            vCMi  = (mv1i + mv2i) / (m1 + m2);

            //permeability
            C = V(i) * 180.0 * (1 - n(i)) * (1 - n(i)) * eta_0 / (n(i) * grains_d * grains_d);

            //fsfi = (mv1i / m1 - mv2i / m2)*C;
            fd_i[i] = C/(1 + job->dt*C*(1/m1 + 1/m2)) * (mv1i/m1 - mv2i/m2);

            granular_body->nodes->f(i)  -= fd_i[i];
            fluid_body->nodes->f(i)     += fd_i[i];
        }
    }



    //--------------------------------------------------------
    // 5 -- Compute Contact Contribution to Interaction Force
    //--------------------------------------------------------

    //contact force vector
    KinematicVectorArray fc_i = KinematicVectorArray(fluid_body->nodes->x.size(), job->JOB_TYPE);

    //set friction coefficient to granular friction angle
    double mu_f = compressible_breakage_mechanics_sand_model->M_0 / std::sqrt(3.0);

    //compoutation variables
    KinematicVector normal(job->JOB_TYPE);
    double m3;
    double fn1i, ft1i;
    double fn2i, ft2i;
    KinematicVector mv3i(job->JOB_TYPE);
    KinematicVector f1i(job->JOB_TYPE);
    KinematicVector f2i(job->JOB_TYPE);
    KinematicVector s1i(job->JOB_TYPE);
    KinematicVector s2i(job->JOB_TYPE);

    KinematicVector tmpVec(job->JOB_TYPE);

    //look for contacts if there are two bodies
    for (int i = 0; i < contact_normal.size(); i++) {
        //test every node for contact
        if (solid_body->nodes->m[i] > 0 && (granular_body->nodes->m[i] > 0 || fluid_body->nodes->m[i] > 0)){
            f1i.setZero();
            f2i.setZero();

            //contact must involve body 1
            normal = contact_normal(i);
            m1 = solid_body->nodes->m[i];
            m2 = granular_body->nodes->m[i];
            m3 = fluid_body->nodes->m[i];

            mv1i = (solid_body->nodes->mx_t(i)      + job->dt * solid_body->nodes->f(i));
            mv2i = (granular_body->nodes->mx_t(i)   + job->dt * granular_body->nodes->f(i));
            mv3i = (fluid_body->nodes->mx_t(i)      + job->dt * fluid_body->nodes->f(i));

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

            //zero fluid tangential force
            s2i.setZero();
            ft2i = 0;

            if (!std::isfinite(s1i.norm())){
                s1i.setZero();
                ft1i = 0;
            }

            f1i = fn1i*normal + ft1i*s1i;
            f2i = fn2i*normal + ft2i*s2i;

            solid_body->nodes->f[i] -= f1i + f2i;
            granular_body->nodes->f[i] += f1i;
            fluid_body->nodes->f[i] += f2i;
            fc_i[i] = -f1i - f2i;
        } else {
            //zero reported force
            fc_i[i].setZero();
        }
    }


    return;
}


void ExplicitMixtureSolver::updateDensity(Job* job){

    // update integration and material point volumes
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (int i=0;i<job->bodies[b]->points->v.rows();i++) {
            job->bodies[b]->points->v(i) *= std::exp(job->dt * job->bodies[b]->points->L(i).trace());
        }

        //this is new, but maybe useful
        job->bodies[b]->points->updateIntegrators(job,job->bodies[b].get());
    }



    //----------------------------------------------
    // 1 -- Compute Density Rate of Change of Solid
    //----------------------------------------------

    // intermediate variables
    double b = compressible_breakage_mechanics_sand_model->b;
    double a, dadp, Je, phiS, rho_s;
    double trD, trDe;

    if (!is_compressible){
        // Incompressible Solid Constituent
        drhos_dt.setZero();
    } else {
        // Compressible Solid Constituent
        for (int i=0; i<granular_body->points->x.size(); i++){

            // Effective Strain-Rates
            trD     = granular_body->points->L[i].trace();
            trDe    = trD - compressible_breakage_mechanics_sand_model->evDot(i);

            // Elastic Volume Ratio
            Je      = std::sqrt(compressible_breakage_mechanics_sand_model->Be[i].det());

            // Solid Volume Fraction and a(phi_s)
            phiS    = compressible_breakage_mechanics_sand_model->phiS(i);
            a       = std::pow(phiS, b);
            dadp    = b * a / phiS;

            // Solid Constituent Density
            rho_s   = compressible_breakage_mechanics_sand_model->rho_0 * (1.0 + a * (1.0/Je - 1.0));

            // Solid Constituent Density Rate of Change
            drhos_dt(i) = -rho_s * (phiS * dadp * (1.0 - Je) * trD + a * trDe) /
                    (Je + a * (1.0 - Je) + phiS * dadp * (1.0 - Je));
        }
    }



    //-----------------------------------------------
    // 2 -- Compute Density Rate of Change of Fluid
    //-----------------------------------------------

    double m1, m2;
    KinematicVector mv1i, mv2i;
    KinematicVectorArray nMat = KinematicVectorArray(n.rows(), job->JOB_TYPE);

    for (int i=0;i<n.rows();i++){
        // nodal masses
        m1 = granular_body->nodes->m[i];
        m2 = fluid_body->nodes->m[i];

        // nodal momenta
        mv1i = granular_body->nodes->mx_t(i);
        mv2i = fluid_body->nodes->mx_t(i);

        // determine nodal field ((1-n)*vs + n*vw)
        if (m1 > 0 && m2 > 0) {
            mv1i = granular_body->nodes->mx_t(i);
            mv2i = fluid_body->nodes->mx_t(i);

            nMat(i) = (1-n(i))/m1 * mv1i + n(i)/m2 * mv2i;
        } else if (m2 > 0) {
            mv2i = fluid_body->nodes->mx_t(i);
            nMat(i) = 1.0/m2 * mv2i;
        } else {
            nMat(i).setZero();
        }
    }

    // interpolate n to fluid points
    Eigen::VectorXd n_p = fluid_body->S.operate(n, MPMSparseMatrixBase::TRANSPOSED);

    // calculate divergence of ((1-n)*vs + n*vw)
    Eigen::VectorXd pvec = fluid_body->gradS.dot(nMat, MPMSparseMatrixBase::TRANSPOSED);

    // adjust divergence for axisym case
    if (job->JOB_TYPE == job->JOB_AXISYM){
        KinematicVectorArray tmpVec(fluid_body->points->x.size(),fluid_body->points->x.VECTOR_TYPE);
        tmpVec = fluid_body->S.operate(nMat, MPMSparseMatrixBase::TRANSPOSED);
        for (int p=0; p<pvec.rows(); p++){
            pvec(p) += tmpVec(p,0)/fluid_body->points->x(p,0);
        }
    }

    // project drhos_dt onto fluid points
    Eigen::VectorXd nodal_drhos_dt = granular_body->S * drhos_dt;
    Eigen::VectorXd fluid_drhos_dt = fluid_body->S.operate(nodal_drhos_dt, MPMSparseMatrixBase::TRANSPOSED);

    // compute drhof_dt
    for (int i=0; i<fluid_body->points->x.size(); i++){
        drhof_dt(i) = -rho_f(i) / n_p(i) * (pvec(i) + fluid_drhos_dt(i) * (1.0 - n_p(i)));
    }


    //--------------------------------------------
    // 3 -- Apply Equivalent Strain-rate to Fluid
    //--------------------------------------------

    // update fluid density
    for (int i=0; i<fluid_body->points->x.size(); i++){
        if (fluid_body->points->active[i] == 0){
            continue;
        }

        rho_f(i) += job->dt * drhof_dt(i);
    }

    return;
}

void ExplicitMixtureSolver::updateStress(Job* job){

    // Compute Stresses for KNOWN Materials

    //--------------------------
    // 1 -- Update Solid Stress
    //--------------------------

    solid_body->material->calculateStress(job, solid_body, Material::UPDATE);

    //-----------------------------
    // 2 -- Update Granular Stress
    //-----------------------------

    granular_body->material->calculateStress(job, granular_body, Material::UPDATE);

    //--------------------------
    // 3 -- Update Fluid Stress
    //--------------------------

    // compute stresses for the fluid directly!
    if (fluid_model == 0){
        // BarotropicViscousFluid
        double trD;
        MaterialTensor D;

        // compute volume ratio
        Eigen::VectorXd J = Eigen::VectorXd(fluid_body->points->x.size());
        for (int i=0;i<J.size();i++){
            J(i) = rhof_0 / rho_f(i);
        }

        // interpolate n to fluid points
        Eigen::VectorXd n_p = fluid_body->S.operate(n, MPMSparseMatrixBase::TRANSPOSED);

        // compute stresses
        for (int i=0;i<fluid_body->points->x.size();i++){
            if (fluid_body->points->active[i] == 0){
                continue;
            }

            D = 0.5 * (fluid_body->points->L[i] + fluid_body->points->L[i].transpose());
            trD = D.trace();

            //T = 2*mu*D_0 + K*log(J)*I
            //mu = 1 + 5/2 * phiS
            fluid_body->points->T[i] = 2 * (1.0 + 5.0/2.0 * (1.0 - n_p(i))) * barotropic_viscous_fluid_model->eta
                    * (D - (trD/3.0)*MaterialTensor::Identity())
                    + barotropic_viscous_fluid_model->K*std::log(J(i))*MaterialTensor::Identity();
        }

    } else if (fluid_model == 1){
        //TillotsonEOSFluid

        // interpolate n to fluid points
        Eigen::VectorXd n_p = fluid_body->S.operate(n, MPMSparseMatrixBase::TRANSPOSED);

        //true density rate and effective density
        double trDp, rhoBar;

        //material state for internal variables
        TillotsonEOSFluid::MaterialState mat_state;

        //tensors for matrix calculation
        MaterialTensor S, Sn, L, D, D0, W;

        //artificial viscosity
        double C0 = std::sqrt(tillotson_eos_fluid_model->A/tillotson_eos_fluid_model->r0);

        for (int i=0;i<fluid_body->points->x.size();i++){
            if (fluid_body->points->active[i] == 0){
                continue;
            }

            // Initialize Material State Structure
            mat_state = TillotsonEOSFluid::MaterialState();
            mat_state.rho  = rho_f(i);
            rhoBar         = fluid_body->points->m(i) / fluid_body->points->v(i);
            mat_state.Tf   = tillotson_eos_fluid_model->Tf(i);
            mat_state.Ef   = tillotson_eos_fluid_model->Ef(i);
            mat_state.S    = fluid_body->points->T[i];
            mat_state.p    = -fluid_body->points->T[i].trace()/3.0;

            // Strain-Rate
            L = fluid_body->points->L[i];

            // Stress Tensor
            S = fluid_body->points->T[i];
            Sn = S;

            // Deformation and Spin Rate
            D = 0.5*(L+L.transpose());
            W = 0.5*(L-L.transpose());

            // Microscopic Volumetric Strain-Rate
            trDp = -drhof_dt(i) / rho_f(i);

            // Estimate Pressure
            mat_state.p = tillotson_eos_fluid_model->PressureFromMaterialState(mat_state);

            // Estimate Material Stress
            D0 = D - D.trace()/3.0 * MaterialTensor::Identity();
            mat_state.S = 2.0 * (1.0 + 5.0/2.0 * (1.0 - n_p(i))) * tillotson_eos_fluid_model->eta0 * D0;
            mat_state.S -= mat_state.p * MaterialTensor::Identity();

            // Check for Cavitation
            if (mat_state.p >= tillotson_eos_fluid_model->PCav) {

                // Add Artificial Viscosity in Compression
                if (tillotson_eos_fluid_model->use_artificial_viscosity && trDp < 0) {
                    mat_state.S += tillotson_eos_fluid_model->r0 * C0 * tillotson_eos_fluid_model->h_i * trDp * MaterialTensor::Identity();
                }

                // Compute Temperature/Energy Increment
                // Shear Contribution to Energy Rate
                mat_state.Ef += 0.5 * (Sn.dot(D0) + mat_state.S.dot(D0)) * job->dt / rhoBar;
                // Pressure Contribution to Energy Rate
                mat_state.Ef += 0.5 * (Sn.trace() + mat_state.S.trace()) / 3.0 * trDp * job->dt / mat_state.rho;

                //mat_state.Ef += 0.5 * (Sn.dot(D) + mat_state.S.dot(D)) * job->dt / mat_state.rho;
                mat_state.Tf = tillotson_eos_fluid_model->TemperatureFromMaterialState(mat_state);

                // Final Stress Update
                mat_state.p = tillotson_eos_fluid_model->PressureFromMaterialState(mat_state);

                // Estimate Material Stress
                fluid_body->points->T[i] = 2.0 * (1.0 + 5.0/2.0 * (1.0 - n_p(i))) * tillotson_eos_fluid_model->eta0 * D0;
                fluid_body->points->T[i] -= mat_state.p * MaterialTensor::Identity();

                // Final State Updtae
                tillotson_eos_fluid_model->Tf(i) = mat_state.Tf;
                tillotson_eos_fluid_model->Ef(i) = mat_state.Ef;
                tillotson_eos_fluid_model->rho(i) = rho_f(i);
                tillotson_eos_fluid_model->J(i)   = tillotson_eos_fluid_model->r0 / mat_state.rho;

                // Add Artificial Viscosity in Compression
                if (tillotson_eos_fluid_model->use_artificial_viscosity && trDp < 0) {
                    fluid_body->points->T[i] += tillotson_eos_fluid_model->r0 * C0 * tillotson_eos_fluid_model->h_i * trDp * MaterialTensor::Identity();
                }

            } else {

                // Fluid has Cavitated

                // Set Stress to Cavitation Pressure
                fluid_body->points->T[i] = -tillotson_eos_fluid_model->PCav * MaterialTensor::Identity();

                // Set Temperature to T0
                tillotson_eos_fluid_model->Tf(i) = mat_state.Tf;
                tillotson_eos_fluid_model->Ef(i) = 0.0;
                tillotson_eos_fluid_model->rho(i) = mat_state.rho;
                tillotson_eos_fluid_model->J(i) = tillotson_eos_fluid_model->r0 / mat_state.rho;
            }
        }

    } else {
        std::cerr << "Something has gone terribly wrong in ExplicitMixtureSolver!\n";
        std::cerr << "    fluid_model = " << fluid_model << std::endl;
    }

    return;
}
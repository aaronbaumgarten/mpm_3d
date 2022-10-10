//
// Created by aaron on 5/15/18.
// isolin.cpp
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <dlfcn.h>

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "materials.hpp"

/*----------------------------------------------------------------------------*/
//initialize assuming that general properties have been assigned correctly
//fp64_props etc. have been filled by configuration object
void BreakageMechanicsSand::init(Job* job, Body* body){
    if (fp64_props.size() < 2){
        std::cout << fp64_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 13 properties defined (pr, K, G, Ec, M_0, gamma, phi_l, phi_u, l, u, theta, rho_0, B_0, phi_0).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        pr      = fp64_props[0];
        K       = fp64_props[1];
        G       = fp64_props[2];
        Ec      = fp64_props[3];
        M_0     = fp64_props[4];
        g       = fp64_props[5];
        phi_l   = fp64_props[6];
        phi_u   = fp64_props[7];
        l       = fp64_props[8];
        u       = fp64_props[9];
        theta   = fp64_props[10];
        rho_0   = fp64_props[11];

        double B_0 = fp64_props[12];
        double phi_0 = fp64_props[13];

        std::cout << "Material properties ("
                  << "pr = " << pr << " Pa, "
                  << "K = " << K << ", "
                  << "G = " << G << ", "
                  << "Ec = " << Ec << " Pa, "
                  << "M_0 = " << M_0 << ", "
                  << "gamma = " << g << ", "
                  << "phi_l = " << phi_l << ", "
                  << "phi_u = " << phi_u << ", "
                  << "l = " << l << ", "
                  << "u = " << u << ", "
                  << "theta = " << theta << ", "
                  << "rho_0 = " << rho_0 << " kg/m^3, "
                  << "B_0 = " << B_0 << ", "
                  << "phi_0 = " << phi_0
                  << ").\n";

        //initialize B, phi, Be
        B = B_0 * Eigen::VectorXd::Ones(body->points->x.size());
        phi = phi_0 * Eigen::VectorXd::Ones(body->points->x.size());
        Be = MaterialTensorArray(body->points->x.size());

        //initialize scalar outputs
        evDot = Eigen::VectorXd(body->points->x.size());
        esDot = Eigen::VectorXd(body->points->x.size());
        BDot = Eigen::VectorXd(body->points->x.size());

        //debug variables
        kVec = Eigen::VectorXd(body->points->x.size());
        yVec = Eigen::VectorXd(body->points->x.size());

        //fill Be with initial deformation
        for (int i=0; i<body->points->x.size(); i++){
            Be(i) = MaterialTensor::Identity();
        }
    }

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//calculate stress state given prior state and strain-rate
void BreakageMechanicsSand::calculateStress(Job* job, Body* body, int SPEC){
    //convergence criteria
    double AbsTOL = 1e-7;
    int MaxIter = 50;

    //tensors for matrix calculation
    MaterialTensor T, T0, L, D, W;
    MaterialTensor Ben, BeA, BeB, BeC;
    MaterialTensor Ee, Ee0, Be0, BeS, BeRate, BeRate0;

    //state for matrix calculations
    MaterialState mat_state, mat_stateA, mat_stateB, mat_stateC;

    //scalars for matrix calculation
    std::vector<double> pq, Rates;
    double p, q, EB, evRate, esRate, BRate, phi_max;
    double lambda, lambdaA, lambdaB, lambdaC, lambdaMAX, lambdaMAXEST;
    double y, y0, yA, yB, yC;

    //iteration variables
    int k;
    bool tensileYield, compactionYield;

    for (int i=0;i<body->points->x.size();i++){
        if (body->points->active[i] == 0){
            continue;
        }
        
        // Set Plasticity Flags
        tensileYield = false;
        compactionYield = false;

        // Initialize Material State Structure
        mat_state = MaterialState();
        mat_state.B = B(i);
        mat_state.phi = phi(i);
        mat_state.rho = body->points->m(i) / body->points->v(i);

        // Strain-Rate
        L = body->points->L[i];

        // Stress Tensor
        T = body->points->T[i];

        // Deformation and Spin Rate
        D = 0.5*(L+L.transpose());
        W = 0.5*(L-L.transpose());

        // Deformation at n
        Ben = Be[i];

        // Estimated Deformation w/o Plasticity
        BeS = Be[i] + job->dt * (L * Ben + Ben * L.transpose());
        Ee = 0.5 * (BeS - MaterialTensor::Identity());
        Ee0 = Ee - Ee.trace()/3.0 * MaterialTensor::Identity();
        Be0 = BeS - BeS.trace()/3.0 * MaterialTensor::Identity();

        // Compute Strain Invariants
        mat_state.ev = Ee.trace();
        mat_state.es = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

        // Check for Yielding in Tension
        if (mat_state.ev >= 0){
            // No Tensile Elastic Strains
            tensileYield = true;
        }

        // Check for Yielding in Compaction
        phi_max = phi_u * std::pow(1.0 - mat_state.B, u);
        if (mat_state.phi > phi_max){
            // No Stresses in Under-Compacted State
            compactionYield = true;
        }

        // Update Material State
        if (tensileYield){
            // Material has Yielded in Tension
            // - no stress
            body->points->T[i].setZero();
            // - zero elastic strain
            Be[i].setIdentity();
            // - porosity evolution continues
            phi(i) += job->dt * D.trace() * (1 - phi(i));

            // Set Scalar Outputs and History Variables
            evDot(i) = D.trace();
            esDot(i) = std::sqrt(2.0 / 3.0 * D.deviator().dot(D.deviator()));
            BDot(i) = 0;
            yVec(i) = 0;
            kVec(i) = 0;

        } else if (compactionYield) {
            // Material has Yielded in Under-Compacted State
            // - no stress
            body->points->T[i].setZero();
            // - zero elastic strain
            Be[i].setIdentity();
            // - porosity evolution continues
            phi(i) += job->dt * D.trace() * (1 - phi(i));

            // Set Scalar Outputs and History Variables
            evDot(i) = D.trace();
            esDot(i) = std::sqrt(2.0 / 3.0 * D.deviator().dot(D.deviator()));
            BDot(i) = 0;
            yVec(i) = 0;
            kVec(i) = 0;

        } else {
            // Material is NOT in Tension or Under-Compacted

            // Estimate Stress State
            T = CauchyStressFromMaterialStateandDeformation(mat_state, BeS);
            T0 = T - T.trace() / 3.0 * MaterialTensor::Identity();

            // Estimate True P, Q
            p = -T.trace() / 3.0;
            q = std::sqrt(3.0 / 2.0 * T0.dot(T0));

            // Estimate Breakage Energy
            EB = EBFromMaterialState(mat_state);

            // Get Relative Deformation Rates
            Rates = RelativePlasticityRatesFromMaterialStateandDeformation(mat_state, BeS);
            evRate = Rates[0];
            esRate = Rates[1];
            BRate = Rates[2];

            // Estimate Maximum Plastic Step Size
            BeRate = MaterialTensor();
            if (q > 0) {
                BeRate = 3.0 / 2.0 * esRate / q * (T0 * Ben + Ben * T0) + 2.0 / 3.0 * evRate * Ben;
            } else {
                BeRate = 2.0 / 3.0 * evRate * Ben;
            }

            lambdaMAX = 1e10; //a ludicrous number
            lambdaMAXEST = lambdaMAX;

            // Check for Plastic Compaction/Dilation
            if (BeRate.trace() < 0) {
                lambdaMAX = (BeS.trace() - 3) / BeRate.trace();
            } else {
                lambdaMAX = BeS.trace() / BeRate.trace();
            }

            // Check for Shearing Reversal
            BeRate0 = BeRate - BeRate.trace() / 3.0 * MaterialTensor::Identity();
            if (Be0.dot(BeRate0) > 0){
                lambdaMAXEST = (Be0.dot(Be0)) / (Be0.dot(BeRate0));
            }
            if (lambdaMAXEST < lambdaMAX){
                lambdaMAX = lambdaMAXEST;
            }

            // Check for Over-Breakage
            if (BRate > 0){
                lambdaMAXEST = (1 - mat_state.B) / BRate;
            }
            if (lambdaMAXEST < lambdaMAX){
                lambdaMAX = lambdaMAXEST;
            }
            
            // Check for Valid Limits
            /*
            if (lambdaMAX < 0){
                std::cout << "Hmmm... lambdaMAX = " << lambdaMAX
                          << ", evRate = " << evRate
                          << ", esRate = " << esRate
                          << ", BRate = " << BRate
                          << ", i = " << i << std::endl;
            }
            */

            // Compute Yield Function Value
            y = YieldFunctionFromMaterialStateandDeformation(mat_state, BeS);
            y0 = y;

            // Initial Guess for Plastic Deformation
            lambda = 0;

            // Initialize Counter
            k = 0;

            // Check if Material has Yielded
            if (y > 0){

                // Use Secant Method
                lambdaA = lambda;               //Plastic Step Lengths
                lambdaB = 1e-3 * lambdaMAX;
                lambdaC = lambdaB;
                yA = 1;                         //Yield Function Values
                yB = 10;
                yC = 100;

                // Compute Material State and Yield Function for A
                mat_stateA      = mat_state;
                mat_stateA.B    = mat_stateA.B      + lambdaA * BRate;
                mat_stateA.phi  = mat_stateA.phi    + lambdaA * evRate * (1.0 - mat_stateA.phi);
                BeA             = BeS               - lambdaA * BeRate;
                Ee              = 0.5 * (BeA - MaterialTensor::Identity());
                Ee0             = Ee - Ee.trace()/3.0 * MaterialTensor::Identity();
                mat_stateA.ev   = Ee.trace();
                mat_stateA.es   = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

                yA = YieldFunctionFromMaterialStateandDeformation(mat_stateA, BeA);

                // Compute Material State and Yield Function for B
                mat_stateB      = mat_state;
                mat_stateB.B    = mat_stateB.B      + lambdaB * BRate;
                mat_stateB.phi  = mat_stateB.phi    + lambdaB * evRate * (1.0 - mat_stateA.phi);
                BeB             = BeS               - lambdaB * BeRate;
                Ee              = 0.5 * (BeB - MaterialTensor::Identity());
                Ee0             = Ee - Ee.trace()/3.0 * MaterialTensor::Identity();
                mat_stateB.ev   = Ee.trace();
                mat_stateB.es   = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

                yB = YieldFunctionFromMaterialStateandDeformation(mat_stateB, BeB);

                // Ensure Yield Function Values Differ
                while (std::abs(yA) > AbsTOL && std::abs(yB) > AbsTOL && std::abs(yA - yB) < AbsTOL && k < MaxIter){
                    // Increment Counter
                    k++;

                    // Double lambdaB
                    lambdaB *= 2.0;

                    // Recompute Material State and Yield Function for B
                    mat_stateB      = mat_state;
                    mat_stateB.B    = mat_stateB.B      + lambdaB * BRate;
                    mat_stateB.phi  = mat_stateB.phi    + lambdaB * evRate * (1.0 - mat_stateB.phi);
                    BeB             = BeS               - lambdaB * BeRate;
                    Ee              = 0.5 * (BeB - MaterialTensor::Identity());
                    Ee0             = Ee - Ee.trace()/3.0 * MaterialTensor::Identity();
                    mat_stateB.ev   = Ee.trace();
                    mat_stateB.es   = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

                    yB = YieldFunctionFromMaterialStateandDeformation(mat_stateB, BeB);

                }

                // Find Plastic Step Length that Solves y = 0
                while (std::abs(yA) > AbsTOL && std::abs(yB) > AbsTOL && k < MaxIter){

                    // Increment k
                    k++;

                    // Choose Bisection or Secant Method
                    if (k > MaxIter / 2.0 && yA * yB < 0){
                        // Recursive Relation for Bisection Method
                        lambdaC = 0.5 * (lambdaA + lambdaB);
                    } else {
                        // Recursive Relation for Secant Method
                        lambdaC = (lambdaA * yB - lambdaB * yA) / (yB - yA);
                    }

                    // If yA and yB Equal, Secant Method and Bisection Will Fail
                    if (std::abs(yA) > AbsTOL && std::abs(yB) > AbsTOL && std::abs(yA - yB) < AbsTOL){
                        // Double Larger Lambda
                        lambdaC = 2.0 * std::max(lambdaA, lambdaB);
                    }

                    // Ensure lambdaC Bounded Below
                    if (lambdaC < 0){
                        lambdaC = 0;
                    }

                    // If lambdaC == lambdaB == 0, then Secant Method Fails
                    if (lambdaC == lambdaB && lambdaC == 0){
                        // Double Larger Lambda
                        lambdaC = 2.0 * lambdaA;
                    }
                    
                    // Ensure lambdaC Bounded Above
                    if (lambdaC > lambdaMAX){
                        lambdaC = lambdaMAX;
                    } 

                    // Compute Material State and Yield Function for C
                    mat_stateC      = mat_state;
                    mat_stateC.B    = mat_stateC.B      + lambdaC * BRate;
                    mat_stateC.phi  = mat_stateC.phi    + lambdaC * evRate * (1.0 - mat_stateC.phi);
                    BeC             = BeS               - lambdaC * BeRate;
                    Ee              = 0.5 * (BeC - MaterialTensor::Identity());
                    Ee0             = Ee - Ee.trace()/3.0 * MaterialTensor::Identity();
                    mat_stateC.ev   = Ee.trace();
                    mat_stateC.es   = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

                    yC = YieldFunctionFromMaterialStateandDeformation(mat_stateC, BeC);

                    // Update Guesses
                    if (yA * yC < 0){
                        lambdaA = lambdaA;
                        yA = yA;
                        lambdaB = lambdaC;
                        yB = yC;
                    } else {
                        lambdaA = lambdaB;
                        yA = yB;
                        lambdaB = lambdaC;
                        yB = yC;
                    }

                }

                // Converged or Reached MaxIter
                if (std::abs(yA) < std::abs(yC)){
                    lambda = lambdaA;
                    y = yA;
                } else {
                    lambda = lambdaC;
                    y = yC;
                }

                // Check for Valid Final State
                if (!std::isfinite(y)){
                    // this implies y = yC and yC is not finite
                    // need to check yA
                    if ((std::abs(y0) < std::abs(yA)) || !std::isfinite(yA)){
                        // either y0 is less than yA or yA is not finite
                        lambda = 0;
                        y = y0;
                    } else {
                        // yA is finite and less than y0
                        lambda = lambdaA;
                        y = yA;
                    }
                }
                
                // Check for Maximum Iterations
                /*
                if (k >= MaxIter && std::abs(y) > 1e3 * AbsTOL){
                    std::cout << "MaxIter Reached: yA = " << yA
                              << ", yB = " << yB 
                              << ", lambdaA = " << lambdaA
                              << ", lambdaB = " << lambdaB
                              << ", lambdaMAX = " << lambdaMAX
                              << ", BRate = " << BRate
                              << ", evRate = " << evRate
                              << ", esRate = " << esRate
                              << ", i = " << i << std::endl;
                }
                */
                

            } else {
                // No Plastic Deformation
            }

            // Final State Update
            mat_state.B     = mat_state.B      + lambda * BRate;
            mat_state.phi   = mat_state.phi    + lambda * evRate * (1.0 - mat_state.phi);
            Be(i)           = BeS              - lambda * BeRate;
            Ee              = 0.5 * (Be(i) - MaterialTensor::Identity());
            Ee0             = Ee - Ee.trace()/3.0 * MaterialTensor::Identity();
            mat_state.ev    = Ee.trace();
            mat_state.es    = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

            y = YieldFunctionFromMaterialStateandDeformation(mat_state, Be(i));

            // Compute Stress
            body->points->T[i] = CauchyStressFromMaterialStateandDeformation(mat_state, Be(i));

            // Porosity and Breakage
            phi(i)  = mat_state.phi;
            B(i)    = mat_state.B;

            // Set Scalar Outputs and History Variables
            evDot(i)    = lambda * evRate / job->dt;
            esDot(i)    = lambda * esRate / job->dt;
            BDot(i)     = lambda * BRate / job->dt;
            yVec(i)     = y;
            kVec(i)     = k;
        }

    }

    return;
}


/*----------------------------------------------------------------------------*/
//define stress assignement for consistency with history dependent materials
void BreakageMechanicsSand::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    std::cout << "WARNING: BreakageMechanicsSand does not implement assignStress(). Assigning pressure instead." << std::endl;
    assignPressure(job, body, (-stressIN.trace()/3.0), idIN, SPEC);
    return;
}


/*----------------------------------------------------------------------------*/
//define pressure assignement for consistency with history dependent materials
void BreakageMechanicsSand::assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    MaterialTensor tmp = body->points->T[idIN];
    body->points->T[idIN] = tmp - (1.0/3.0 * tmp.trace() + pressureIN)*MaterialTensor::Identity();

    // Approximate Pressure Assignment
    double rho = body->points->m(idIN) / body->points->v(idIN);
    double ev = -std::sqrt(pressureIN * rho_0 * 4 / (K * K * pr * rho * (1.0 - theta * B(idIN))));
    Be(idIN) = (2.0 / 3.0 * ev + 1) * MaterialTensor::Identity();

    // Create MaterialState Structure
    MaterialState mat_state = MaterialState();
    mat_state.B = B(idIN);
    mat_state.phi = phi(idIN);
    mat_state.rho = rho;
    mat_state.ev = ev;
    mat_state.es = 0;

    // Compute Stress
    body->points->T[idIN] = CauchyStressFromMaterialStateandDeformation(mat_state, Be(idIN));

    return;
}


/*----------------------------------------------------------------------------*/
//frame and state writing
void BreakageMechanicsSand::writeFrame(Job* job, Body* body, Serializer* serializer){

    // Compute Material Strain
    MaterialTensorArray Ee = MaterialTensorArray(body->points->x.size());
    for (int i=0; i<body->points->x.size(); i++){
        Ee[i] = 0.5 * (Be[i] - MaterialTensor::Identity());
    }
    
    // Compute Stress Invaraints
    Eigen::VectorXd p = Eigen::VectorXd(body->points->x.size());
    Eigen::VectorXd q = Eigen::VectorXd(body->points->x.size());
    MaterialTensor T0;
    for (int i=0; i<body->points->x.size(); i++){
        p(i) = -body->points->T[i].trace()/3.0;
        T0 = body->points->T[i] + p(i) * MaterialTensor::Identity();
        q(i) = std::sqrt(3.0/2.0 * T0.dot(T0));
    }

    // Write Scalar Quantities
    serializer->writeScalarArray(p, "p");
    serializer->writeScalarArray(q, "q");
    serializer->writeScalarArray(phi, "phi");
    serializer->writeScalarArray(B, "B");
    serializer->writeScalarArray(evDot, "evDot");
    serializer->writeScalarArray(esDot, "esDot");
    serializer->writeScalarArray(BDot, "BDot");
    serializer->writeScalarArray(yVec, "yVec");
    serializer->writeScalarArray(kVec, "kVec");

    // Write Tensor Quantities
    serializer->writeTensorArray(Be, "Be");
    serializer->writeTensorArray(Ee, "Ee");

    // Double Projection
    Eigen::VectorXd v_i(body->nodes->x.size());
    Eigen::VectorXd nvec(body->nodes->x.size());
    Eigen::VectorXd qvec(body->points->x.size());
    Eigen::VectorXd pvec(body->points->x.size());

    // Multiply P by Volume
    for (int i=0;i<body->points->x.size();i++) {
        if (body->points->active[i] == 0) {
            pvec(i) = 0;
            continue;
        }

        pvec(i) = p(i) * body->points->v(i);
    }
    // Integrate P onto Nodes
    nvec = body->S * pvec;
    // Integrate Volume onto Nodes
    v_i = body->S * body->points->v;
    // Volume Average P on Nodes
    for (int i=0;i<nvec.size();i++){
        nvec(i) = nvec(i) / v_i(i);
    }
    // Interpolate P onto Points
    pvec = body->S.operate(nvec, MPMSparseMatrixBase::TRANSPOSED);

    // Multiply Q by Volume
    for (int i=0;i<body->points->x.size();i++) {
        if (body->points->active[i] == 0) {
            qvec(i) = 0;
            continue;
        }

        qvec(i) = q(i) * body->points->v(i);
    }
    // Integrate Q onto Nodes
    nvec = body->S * qvec;
    // Integrate Volume onto Nodes
    v_i = body->S * body->points->v;
    // Volume Average P on Nodes
    for (int i=0;i<nvec.size();i++){
        nvec(i) = nvec(i) / v_i(i);
    }
    // Interpolate Q onto Points
    qvec = body->S.operate(nvec, MPMSparseMatrixBase::TRANSPOSED);

    serializer->writeScalarArray(pvec, "p_smooth");
    serializer->writeScalarArray(qvec, "q_smooth");

    return;
}

std::string BreakageMechanicsSand::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int BreakageMechanicsSand::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//MATLAB Helper Functions
double BreakageMechanicsSand::CriticalStatePorosityFromMaterialState(MaterialState stateIN){
//Compute phi_cs from Material State
//  Compute critical state porosity by inverting critical state function from Tengattini et al (2016)

    // Compute EB
    double EB = EBFromMaterialState(stateIN);

    // Compute Limiting Porosities
    double phi_min = phi_l * std::pow(1.0 - stateIN.B, l);
    double phi_max = phi_u * std::pow(1.0 - stateIN.B, u);

    // Compute Critical State
    double phiOUT = phi_max - (phi_max - phi_min) * (1.0 - stateIN.B) * std::sqrt(EB/Ec) / g;

    return phiOUT;
}

double BreakageMechanicsSand::EBFromMaterialState(MaterialState stateIN){
//Compute EB from Scalar Elastic Strains
//  Returns the breakage energy value EB computed from Material State

    //double A = -K * stateIN.ev / 2.0;
    double A = -K / 2.0 * (stateIN.ev + e0 * std::log(1 - stateIN.ev/e0));
    if (A < 0) {
        A = 0;
    }

    double EB = (stateIN.rho / rho_0) * theta * pr * (2.0 * A * A * A / (3.0 * K)
                                                      + 3.0 * G * A * stateIN.es * stateIN.es / 2.0
                                                      + G * std::sqrt(3.0 * K * G) * stateIN.es * stateIN.es * stateIN.es / 4.0);

    return EB;

}


double BreakageMechanicsSand::FFromMaterialState(MaterialState stateIN){
//Compute F from Material State
//  Compute critical state variable F from Tengattini et al. (2016) using material state

    //Compute EB
    double EB = EBFromMaterialState(stateIN);

    // Compute Limiting Porosities
    double phi_min = phi_l * std::pow(1.0 - stateIN.B, l);
    double phi_max = phi_u * std::pow(1.0 - stateIN.B, u);

    // Compute Relative Porosity
    double tau = (phi_max - stateIN.phi) / (phi_max - phi_min);

    return std::sqrt(EB/Ec) * (1.0 - stateIN.B) - g * tau;
}

std::vector<double> BreakageMechanicsSand::PQFromMaterialState(MaterialState stateIN){
//Compute P,Q from Scalar Elastic Strains
//  Compute pressure and shear stress from material state

    // Output Vector
    std::vector<double> pq(2);

    // Compression Coeff.
    //double A = -K * stateIN.ev / 2.0;
    double A = -K / 2.0 * (stateIN.ev + e0 * std::log(1 - stateIN.ev/e0));
    if (A < 0) {
        A = 0;
    }

    // Compute Pressure
    //pq[0] = (stateIN.rho/rho_0) * pr * (1.0 - theta * stateIN.B) * (A*A + 3.0 * G * K * stateIN.es * stateIN.es / 4.0);
    pq[0] = (stateIN.rho/rho_0) * pr * (1.0 - theta * stateIN.B) * (2.0 * A * A / K + 3.0 * G * stateIN.es * stateIN.es / 2.0)
            * (-K / 2.0 * stateIN.ev / (e0 - stateIN.ev));

    // Compute Shear Stress
    pq[1] = (stateIN.rho/rho_0) * pr * (1.0 - theta * stateIN.B) * 3.0 * G * (A * stateIN.es + std::sqrt(3.0 * G * K) * stateIN.es * stateIN.es / 4.0);

    // Zero Out Stresses for Positive Volumetric Strain
    if (stateIN.ev > 0){
        pq[0] = 0;
        pq[1] = 0;
    }

    return pq;
}

std::vector<double> BreakageMechanicsSand::RelativePlasticityRatesFromMaterialStateandDeformation(MaterialState stateIN, MaterialTensor Be){
//Compute relative rates of deformation for ev, es, and B from material state and material deformation

    // Output Vector
    std::vector<double> Rates(3);

    // Compute Material Stress State
    MaterialTensor S = CauchyStressFromMaterialStateandDeformation(stateIN, Be);

    // Compute Pressure and Shear Stress in Deformed Reference Frame
    double p = -S.trace()/3.0;
    MaterialTensor S0 = S + p*MaterialTensor::Identity();
    double q = std::sqrt(3.0/2.0 * S0.dot(S0));

    // Compute EB
    double EB = EBFromMaterialState(stateIN);

    // Compute F
    double F = FFromMaterialState(stateIN);

    // Compute Limiting Porosities
    double phi_min = phi_l * std::pow(1.0 - stateIN.B, l);
    double phi_max = phi_u * std::pow(1.0 - stateIN.B, u);

    // Compute Relative Porosity
    double tau = (phi_max - stateIN.phi) / (phi_max - phi_min);

    // Compute Critical Porosity
    double phi_cs = CriticalStatePorosityFromMaterialState(stateIN);

    // Compute Maximum Friction Angle
    double tp = M_PI / 15.0 + std::asin(3.0 * M_0 / (6.0 + M_0));

    // Ensure Positive Pressure
    if (p <= 0){
        // Zero Pressure Yielding
        Rates[0] = 1; //evRate
        Rates[1] = 1; //esRate
        Rates[2] = 0; //BRate
    } else {
        // Dilation Angle
        double M_d = - g * (stateIN.phi - phi_cs) / (phi_max - phi_min) * (6.0 * std::sin(tp) / (3.0 - std::sin(tp)) - M_0);

        // Coupling Angle
        double w = M_PI / 2.0 * (1 - tau);

        // Plasticity Rates
        Rates[0] = -EB * (1.0 - stateIN.B) * (1.0 - stateIN.B) * std::sin(w) * std::sin(w) / Ec
                   * p / (p*p + q*q)
                   + q * M_d / ((M_0 + M_d)*p * (M_0 + M_d)*p);
        Rates[1] = EB * (1.0 - stateIN.B) * (1.0 - stateIN.B) * std::sin(w) * std::sin(w) / Ec
                   * q / (p*p + q*q)
                   + q / ((M_0 + M_d)*p * (M_0 + M_d)*p);
        Rates[2] = (1.0 - stateIN.B) * (1.0 - stateIN.B) * std::cos(w) * std::cos(w) / Ec;

        // Check for Compaction
        /*
        if (F <= 0) {
            Rates[0] = q * M_d / ((M_0 + M_d)*p * (M_0 + M_d)*p);   //evRate
            Rates[1] = q / ((M_0 + M_d)*p * (M_0 + M_d)*p);         //esRate
            Rates[2] = 0;                                           //BRate
        } else {
            Rates[0] = -EB * (1.0 - stateIN.B) * (1.0 - stateIN.B) * F * std::sin(w) * std::sin(w) / (p * Ec)
                       + q * M_d / ((M_0 + M_d)*p * (M_0 + M_d)*p);
            Rates[1] = q / ((M_0 + M_d)*p * (M_0 + M_d)*p);
            Rates[2] = (1.0 - stateIN.B) * (1.0 - stateIN.B) * F * std::cos(w) * std::cos(w) / Ec;
        }
         */
    }

    return Rates;

}


double BreakageMechanicsSand::YieldFunctionFromMaterialStateandDeformation(MaterialState stateIN, MaterialTensor Be){
//Compute y from material state

    // Compute Material Stress State
    MaterialTensor S = CauchyStressFromMaterialStateandDeformation(stateIN, Be);

    // Compute Pressure and Shear Stress in Deformed Reference Frame
    double p = -S.trace()/3.0;
    MaterialTensor S0 = S + p*MaterialTensor::Identity();
    double q = std::sqrt(3.0/2.0 * S0.dot(S0));

    // Compute EB
    double EB = EBFromMaterialState(stateIN);

    // Compute F
    double F = FFromMaterialState(stateIN);

    // Compute Limiting Porosities
    double phi_min = phi_l * std::pow(1.0 - stateIN.B, l);
    double phi_max = phi_u * std::pow(1.0 - stateIN.B, u);

    // Compute Relative Porosity
    double tau = (phi_max - stateIN.phi) / (phi_max - phi_min);

    // Compute Critical Porosity
    double phi_cs = CriticalStatePorosityFromMaterialState(stateIN);

    // Compute Maximum Friction Angle
    double tp = M_PI / 15.0 + std::asin(3.0 * M_0 / (6.0 + M_0));

    // Dilation Angle
    double M_d = - g * (stateIN.phi - phi_cs) / (phi_max - phi_min) * (6.0 * std::sin(tp) / (3.0 - std::sin(tp)) - M_0);

    // Coupling Angle
    double w = M_PI / 2.0 * (1 - tau);

    // Yield f'n
    return EB * (1.0 - stateIN.B) * (1.0 - stateIN.B) / Ec + q * q / ((M_0 + M_d)*p * (M_0 + M_d)*p) - 1;

}

MaterialTensor BreakageMechanicsSand::CauchyStressFromMaterialStateandDeformation(MaterialState stateIN, MaterialTensor Be){
    // Compute Cauchy Stress from Material State and Deformation

    // Compute Deviator of Be
    MaterialTensor Be0 = Be - Be.trace()/3.0 * MaterialTensor::Identity();

    // Compute Left-Handed Green--St. Venant Strain Tensor
    // Invariants Identical to Invariants of Right-Handed Green--St. Venant Strain Tensor
    MaterialTensor Ee = 0.5 * (Be - MaterialTensor::Identity());
    MaterialTensor Ee0 = Ee - Ee.trace()/3.0 * MaterialTensor::Identity();

    // Compute Strain Invariants (Override Inputs from stateIN)
    stateIN.ev = Ee.trace();
    stateIN.es = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

    // Compute P,Q from Helmholtz free energy
    std::vector<double> pq = PQFromMaterialState(stateIN);
    double p = pq[0];
    double q = pq[1];

    // Compute S
    MaterialTensor S = MaterialTensor();
    if (pq[1] > 0){
        //Non-Zero Shear Stress
        S = (2.0 * q / (3.0 * stateIN.es) * Ee0 - p * MaterialTensor::Identity()) * Be;
    } else {
        //Zero Shear Stress
        S = -p * Be;
    }

    return S;
}

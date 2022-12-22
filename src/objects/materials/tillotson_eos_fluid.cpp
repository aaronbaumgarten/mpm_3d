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
void TillotsonEOSFluid::init(Job* job, Body* body){
    if (fp64_props.size() < 14){
        std::cout << fp64_props.size() << " != 14\n";
        fprintf(stderr, "%s:%s:", __FILE__, __func__);

        std::cerr << "Need at least 14 properties defined:\n";
        std::cerr << "    r0     [kg/m^3] Reference Density\n";
        std::cerr << "    rIV    [kg/m^3] Density of Incipient Vaporization\n";
        std::cerr << "    a      [--]     (?)\n";
        std::cerr << "    b      [--]     (?)\n";
        std::cerr << "    A      [Pa]     First-Order Bulk Modulus\n";
        std::cerr << "    B      [Pa]     Second-Order Bulk Modulus\n";
        std::cerr << "    E0     [J/kg]   Reference Energy\n";
        std::cerr << "    alfa   [--]     (?)\n";
        std::cerr << "    beta   [--]     (?)\n";
        std::cerr << "    EIV    [J/kg]   Energy of Incipient Vaporization\n";
        std::cerr << "    ECV    [J/kg]   Energy of Complete Vaporization\n";
        std::cerr << "    cv     [J/kg*K] Heat Capacity at Constant Volume\n";
        std::cerr << "    T0     [K]      Reference Temperature\n";
        std::cerr << "    eta    [Pa*s]   Viscosity\n";

        exit(0);
    } else {
        r0      = fp64_props[0];
        rIV     = fp64_props[1];
        a       = fp64_props[2];
        b       = fp64_props[3];
        A       = fp64_props[4];
        B       = fp64_props[5];
        E0      = fp64_props[6];
        alfa    = fp64_props[7];
        beta    = fp64_props[8];
        EIV     = fp64_props[9];
        ECV     = fp64_props[10];

        cv      = fp64_props[11];
        T0      = fp64_props[12];
        eta     = fp64_props[13];

        std::cout << "Material Propoerties:\n";


        std::cout << "Need at least 14 properties defined:\n";
        std::cout << "    r0     = " << r0 << " [kg/m^3]\n";
        std::cout << "    rIV    = " << rIV << " [kg/m^3]\n";
        std::cout << "    a      = " << a << " [--]\n";
        std::cout << "    b      = " << b << " [--]\n";
        std::cout << "    A      = " << A << " [Pa]\n";
        std::cout << "    B      = " << B << " [Pa]\n";
        std::cout << "    E0     = " << E0 << " [J/kg]\n";
        std::cout << "    alfa   = " << alfa << " [--]\n";
        std::cout << "    beta   = " << beta << " [--]\n";
        std::cout << "    EIV    = " << EIV << " [J/kg]\n";
        std::cout << "    ECV    = " << ECV << " [J/kg]\n";
        std::cout << "    cv     = " << cv << " [J/kg*K]\n";
        std::cout << "    T0     = " << T0 << " [K]\n";
        std::cout << "    eta    = " << eta << " [Pa*s]\n";

        //check for string flags
        for (int s=0; s<str_props.size(); s++){
            if (str_props[s].compare("ISOTHERMAL") == 0){
                //isothermal simulation
                is_adiabatic = false;
            } else if (str_props[s].compare("USE_ARTIFICIAL_VISCOSITY") == 0){
                //artificial velocity
                use_artificial_viscosity = true;
            }
        }

        if (use_artificial_viscosity && fp64_props.size() > 14){
            h_i = fp64_props[15];
            std::cout << "  Using Artificial Viscosity:\n";
            std::cout << "    h_i    = " << h_i << " [m]\n";
        } else if (use_artificial_viscosity){
            std::cout << "  *NOT* Using Artificial Viscosity! No h_i parameter defined.\n";
        }

        if (is_adiabatic){
            std::cout << "  Using Adiabatic Model.\n";
        } else {
            std::cout << "  Using Isothermal Model.\n";
        }

        //initialize internal state
        rho = r0 * Eigen::VectorXd::Ones(body->points->x.size());
        J = Eigen::VectorXd::Ones(body->points->x.size());

        Ef = cv * T0 * Eigen::VectorXd::Ones(body->points->x.size());
        Tf = T0 * Eigen::VectorXd::Ones(body->points->x.size());

        //compute cold energy reference states
        ComputeColdEnergyReferenceStates();
    }

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//calculate stress state given prior state and strain-rate
void TillotsonEOSFluid::calculateStress(Job* job, Body* body, int SPEC){
    //convergence criteria
    double AbsTOL = 1e-7;
    int MaxIter = 50;

    //tensors for matrix calculation
    MaterialTensor S, s0, Sn, L, D, W;
    MaterialTensor Ben;
    MaterialTensor Ee, Ee0, Be0, BeRate, BeRate0;

    //state for matrix calculations
    MaterialState mat_state, mat_stateA, mat_stateB, mat_stateC;

    //scalars for matrix calculation
    std::vector<double> pq, Rates;
    double p, q, EB, evRate, esRate, BRate, phi_max;
    double lambda, lambdaA, lambdaB, lambdaC, lambdaMAX, lambdaMAXEST;
    double y, y0, yA, yB, yC;
    double trDe, a, dadp, phiSTemp, JeTemp;

    //iteration variables
    int k;
    bool tensileYield, compactionYield, validB;

    for (int i=0;i<body->points->x.size();i++){
        if (body->points->active[i] == 0){
            continue;
        }
        
        // Set Plasticity Flags
        tensileYield = false;
        compactionYield = false;

        // Initialize Material State Structure
        mat_state = MaterialState();
        mat_state.B    = B(i);
        mat_state.phi  = 1.0 - phiS(i);
        mat_state.phiP = phiP(i);
        mat_state.rho  = body->points->m(i) / body->points->v(i);
        mat_state.F    = F[i];
        mat_state.Be   = Be[i];
        mat_state.Ts   = Ts(i);
        mat_state.Es   = Es(i);

        // Strain-Rate
        L = body->points->L[i];

        // Stress Tensor
        S = body->points->T[i];
        Sn = S;

        // Deformation and Spin Rate
        D = 0.5*(L+L.transpose());
        W = 0.5*(L-L.transpose());

        // Deformation at n
        Ben = Be[i];

        // Estimated Deformation w/o Plasticity
        mat_state.Be = Be[i] + job->dt * (L * Ben + Ben * L.transpose());
        mat_state.F  = F[i] + job->dt * L * F[i];

        // Estimate Porosity
        mat_state.phi = PorosityFromMaterialState(mat_state);

        // Check for Yielding in Tension
        if (mat_state.Be.trace() > 3 && mat_state.Be.det() > 1){
            // No Stresses in Tension
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
            phiP(i) += job->dt * D.trace() * (1 - phiP(i));
            phiS(i) = 1.0 - mat_state.phi;
            // - energy is only thermal
            Es(i) = cv * Ts(i);
            // - save deformation
            F[i]       = mat_state.F;

            // Set Scalar Outputs and History Variables
            evDot(i) = D.trace();
            esDot(i) = std::sqrt(2.0 / 3.0 * D.deviator().dot(D.deviator()));
            BDot(i) = 0;
            yVec(i) = 0;
            kVec(i) = -1;
            lVec(i) = 0;

        } else if (compactionYield) {
            // Material has Yielded in Under-Compacted State
            // - no stress
            body->points->T[i].setZero();
            // - zero elastic strain
            Be[i].setIdentity();
            // - porosity evolution continues
            phiP(i) += job->dt * D.trace() * (1 - phiP(i));
            phiS(i) = 1.0 - mat_state.phi;
            // - energy is only thermal
            Es(i) = cv * Ts(i);
            // - save deformation
            F[i]       = mat_state.F;

            // Set Scalar Outputs and History Variables
            evDot(i) = D.trace();
            esDot(i) = std::sqrt(2.0 / 3.0 * D.deviator().dot(D.deviator()));
            BDot(i) = 0;
            yVec(i) = 0;
            kVec(i) = -2;
            lVec(i) = 0;

        } else {
            // Material is NOT in Tension or Under-Compacted

            // Estimate Stress State
            S = CauchyStressFromMaterialState(mat_state);
            s0 = S - S.trace() / 3.0 * MaterialTensor::Identity();

            // Estimate True P, Q
            p = -S.trace() / 3.0;
            q = std::sqrt(3.0 / 2.0 * s0.dot(s0));

            // Get Relative Deformation Rates
            Rates = PlasticFlowRulesFromMaterialState(mat_state);
            evRate = Rates[0];
            esRate = Rates[1];
            BRate = Rates[2];

            // Estimate Maximum Plastic Step Size
            BeRate = MaterialTensor();
            if (q > 0) {
                BeRate = 3.0 / 2.0 * esRate / q * (s0 * Ben + Ben * s0) + 2.0 / 3.0 * evRate * Ben;
            } else {
                BeRate = 2.0 / 3.0 * evRate * Ben;
            }

            /*
            if (!std::isfinite(BeRate.det())){
                std::cout << i << " : Uh oh!" << std::endl;
                std::cout << evRate << ", " << esRate << std::endl;
                std::cout << p << ", " << q << std::endl;
                std::cout << Ben(0,0) << ", " << Ben(0,1) << ", " << Ben(0,2) << ", ";
                std::cout << Ben(1,0) << ", " << Ben(1,1) << ", " << Ben(1,2) << ", ";
                std::cout << Ben(2,0) << ", " << Ben(2,1) << ", " << Ben(2,2) << std::endl;
            }
             */

            lambdaMAX = 1e10; //a ludicrous number
            lambdaMAXEST = lambdaMAX;

            // Check for Plastic Compaction/Dilation
            if (BeRate.trace() < 0) {
                lambdaMAX = (mat_state.Be.trace() - 3) / BeRate.trace();
            } else {
                lambdaMAX = mat_state.Be.trace() / BeRate.trace();
            }

            // Check for Shearing Reversal
            BeRate0 = BeRate - BeRate.trace() / 3.0 * MaterialTensor::Identity();
            Be0 = mat_state.Be - mat_state.Be.trace() / 3.0 * MaterialTensor::Identity();
            if (Be0.dot(BeRate0) > 0){
                lambdaMAXEST = (Be0.dot(Be0)) / (Be0.dot(BeRate0));
            }
            if (lambdaMAXEST < lambdaMAX){
                lambdaMAX = lambdaMAXEST;
            }

            // Check for Over-Breakage
            if (BRate > 0){
                lambdaMAXEST = (1.0 - mat_state.B) / BRate;
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
            y = YieldFunctionFromMaterialState(mat_state);
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
                mat_stateA.phiP = mat_stateA.phiP   + lambdaA * evRate * (1.0 - mat_stateA.phiP);
                mat_stateA.Be   = mat_stateA.Be     - lambdaA * BeRate;
                mat_stateA.phi  = PorosityFromMaterialState(mat_stateA);

                yA = YieldFunctionFromMaterialState(mat_stateA);

                // Compute Material State and Yield Function for B
                mat_stateB      = mat_state;
                mat_stateB.B    = mat_stateB.B      + lambdaB * BRate;
                mat_stateB.phiP = mat_stateB.phiP   + lambdaB * evRate * (1.0 - mat_stateB.phiP);
                mat_stateB.Be   = mat_stateB.Be     - lambdaB * BeRate;
                mat_stateB.phi  = PorosityFromMaterialState(mat_stateB);

                yB = YieldFunctionFromMaterialState(mat_stateB);

                // Ensure yB is Valid
                validB = false;
                while (!validB && k<MaxIter){

                    // Set Validity Flag
                    validB = true;

                    // Ensure Finite Yield Function Value
                    if (!std::isfinite(yB)){
                        // Increment Counter
                        k++;

                        // Bisect Back
                        lambdaB = 0.5 * (lambdaA + lambdaB);

                        // Set Validity Flag
                        validB = false;

                        // Recompute Material State and Yield Function for B
                        mat_stateB      = mat_state;
                        mat_stateB.B    = mat_stateB.B      + lambdaB * BRate;
                        mat_stateB.phiP = mat_stateB.phiP   + lambdaB * evRate * (1.0 - mat_stateB.phiP);
                        mat_stateB.Be   = mat_stateB.Be     - lambdaB * BeRate;
                        mat_stateB.phi  = PorosityFromMaterialState(mat_stateB);

                        yB = YieldFunctionFromMaterialState(mat_stateB);

                    }

                    // Ensure Yield Function Values Differ
                    if (std::abs(yA) > AbsTOL && std::abs(yB) > AbsTOL && std::abs(yA - yB) < AbsTOL && k < MaxIter){
                        // Increment Counter
                        k++;

                        // Double lambdaB
                        lambdaB *= 2.0;

                        // Set Validity Flag
                        validB = false;

                        // Recompute Material State and Yield Function for B
                        mat_stateB      = mat_state;
                        mat_stateB.B    = mat_stateB.B      + lambdaB * BRate;
                        mat_stateB.phiP = mat_stateB.phiP   + lambdaB * evRate * (1.0 - mat_stateB.phiP);
                        mat_stateB.Be   = mat_stateB.Be     - lambdaB * BeRate;
                        mat_stateB.phi  = PorosityFromMaterialState(mat_stateB);

                        yB = YieldFunctionFromMaterialState(mat_stateB);

                    }
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
                    mat_stateC.phiP = mat_stateC.phiP   + lambdaC * evRate * (1.0 - mat_stateC.phiP);
                    mat_stateC.Be   = mat_stateC.Be     - lambdaC * BeRate;
                    mat_stateC.phi  = PorosityFromMaterialState(mat_stateC);

                    yC = YieldFunctionFromMaterialState(mat_stateC);


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
            mat_state.B    = mat_state.B      + lambda * BRate;
            mat_state.phiP = mat_state.phiP   + lambda * evRate * (1.0 - mat_state.phiP);
            mat_state.Be   = mat_state.Be     - lambda * BeRate;
            mat_state.phi  = PorosityFromMaterialState(mat_state);
            mat_state.S    = CauchyStressFromMaterialState(mat_state);

            // Compute Elastic Volumetric Strain-Rate
            trDe = D.trace() - lambda * evRate;

            // Add Artificial Viscosity in Compression
            if (use_artificial_viscosity && trDe < 0){
                phiSTemp = 1.0 - mat_state.phi;
                JeTemp = std::sqrt(mat_state.Be.det());
                a = std::pow(phiSTemp, b);
                dadp = b * a / phiSTemp;
                mat_state.S += rho_0 * C0 * h_i * trDe *
                        (a + phiSTemp * dadp * (1.0 - JeTemp)) /
                        (JeTemp + (a + phiSTemp * dadp) * (1 - JeTemp)) *
                        MaterialTensor::Identity();
            }

            // Compute Temperature/Energy Increment
            mat_state.Es += 0.5 * (Sn.dot(D) + mat_state.S.dot(D)) * job->dt / mat_state.rho;
            mat_state.Ts  = TemperatureFromMaterialState(mat_state);

            y = YieldFunctionFromMaterialState(mat_state);

            // Final Stress Update
            body->points->T[i] = CauchyStressFromMaterialState(mat_state);

            // Add Artificial Viscosity in Compression
            if (use_artificial_viscosity && trDe < 0){

                phiSTemp = 1.0 - mat_state.phi;
                JeTemp = std::sqrt(mat_state.Be.det());
                a = std::pow(phiSTemp, b);
                dadp = b * a / phiSTemp;
                body->points->T[i] += rho_0 * C0 * h_i * trDe *
                               (a + phiSTemp * dadp * (1.0 - JeTemp)) /
                               (JeTemp + (a + phiSTemp * dadp) * (1 - JeTemp)) *
                               MaterialTensor::Identity();

                //body->points->T[i] += rho_0 * C0 * h_i * D.trace() * MaterialTensor::Identity();
            }

            // Porosity and Breakage
            B(i)    = mat_state.B;
            phiS(i)    = 1.0 - mat_state.phi;
            phiP(i)    = mat_state.phiP;
            Es(i)      = mat_state.Es;
            Ts(i)      = mat_state.Ts;
            F[i]       = mat_state.F;
            Be[i]      = mat_state.Be;



            // Set Scalar Outputs and History Variables
            evDot(i)    = lambda * evRate / job->dt;
            esDot(i)    = lambda * esRate / job->dt;
            BDot(i)     = lambda * BRate / job->dt;
            yVec(i)     = y;
            kVec(i)     = k;
            lVec(i)     = lambda;
        }

    }

    return;
}


/*----------------------------------------------------------------------------*/
//define stress assignement for consistency with history dependent materials
void TillotsonEOSFluid::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    std::cout << "WARNING: TillotsonEOSFluid does not implement assignStress(). Assigning pressure instead." << std::endl;
    assignPressure(job, body, (-stressIN.trace()/3.0), idIN, SPEC);
    return;
}


/*----------------------------------------------------------------------------*/
//define pressure assignement for consistency with history dependent materials
void TillotsonEOSFluid::assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    MaterialTensor tmp = body->points->T[idIN];
    body->points->T[idIN] = tmp - (1.0/3.0 * tmp.trace() + pressureIN)*MaterialTensor::Identity();

    // Approximate Pressure Assignment
    double rho = body->points->m(idIN) / body->points->v(idIN);
    double ev = -std::sqrt(pressureIN * rho_0 * 4 / (K * K * pr * rho * (1.0 - theta * B(idIN))));
    Be(idIN) = (2.0 / 3.0 * ev + 1) * MaterialTensor::Identity();

    // Create MaterialState Structure
    MaterialState mat_state = MaterialState();
    mat_state.B = B(idIN);
    mat_state.phi = 1.0 - phiS(idIN);
    mat_state.phiP = phiP(idIN);
    mat_state.rho = rho;
    mat_state.Be = Be[idIN];
    mat_state.F  = F[idIN];
    mat_state.Ts = Ts[idIN];
    mat_state.Es = Es[idIN];

    // Compute Stress
    body->points->T[idIN] = CauchyStressFromMaterialState(mat_state);

    return;
}


/*----------------------------------------------------------------------------*/
//frame and state writing
void TillotsonEOSFluid::writeFrame(Job* job, Body* body, Serializer* serializer){

    // Compute Material Strain
    MaterialTensor Ee = MaterialTensor();
    MaterialTensor Ee0;
    Eigen::VectorXd eve = Eigen::VectorXd(body->points->x.size());
    Eigen::VectorXd ese = Eigen::VectorXd(body->points->x.size());
    Eigen::VectorXd Je = Eigen::VectorXd(body->points->x.size());
    for (int i=0; i<body->points->x.size(); i++){
        Ee = 0.5 * (Be[i] - MaterialTensor::Identity());
        Ee0 = Ee.deviator();
        eve(i) = Ee.trace();
        ese(i) = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));
        Je(i)  = std::sqrt(Be[i].det());
    }
    
    // Compute Stress Invaraints
    Eigen::VectorXd p = Eigen::VectorXd(body->points->x.size());
    Eigen::VectorXd q = Eigen::VectorXd(body->points->x.size());
    MaterialTensor s0;
    for (int i=0; i<body->points->x.size(); i++){
        p(i) = -body->points->T[i].trace()/3.0;
        s0 = body->points->T[i] + p(i) * MaterialTensor::Identity();
        q(i) = std::sqrt(3.0 / 2.0 * s0.dot(s0));
    }

    // Write Scalar Quantities
    serializer->writeScalarArray(p, "p");
    serializer->writeScalarArray(q, "q");
    serializer->writeScalarArray(phiS, "phiS");
    serializer->writeScalarArray(phiP, "phiP");
    serializer->writeScalarArray(B, "B");
    serializer->writeScalarArray(Es, "Es");
    serializer->writeScalarArray(Ts, "Ts");
    serializer->writeScalarArray(evDot, "evDot");
    serializer->writeScalarArray(esDot, "esDot");
    serializer->writeScalarArray(BDot, "BDot");
    serializer->writeScalarArray(yVec, "yVec");
    serializer->writeScalarArray(kVec, "kVec");
    serializer->writeScalarArray(lVec, "lVec");
    serializer->writeScalarArray(eve, "eve");
    serializer->writeScalarArray(ese, "ese");
    serializer->writeScalarArray(Je, "Je");

    // Write Tensor Quantities
    serializer->writeTensorArray(Be, "Be");
    serializer->writeTensorArray(F, "F");

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

std::string TillotsonEOSFluid::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int TillotsonEOSFluid::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//Helper Functions

void TillotsonEOSFluid::ComputeColdEnergyReferenceStates(){
//Compute Reference Values for Cold Energy
//  Compute discrete reference cold energy values between Jmin and Jmax to save computation time during simulation

    // State What You Are Doing
    std::cout << "Computing cold energy curve for TillotsonEOSFluid..." << std::endl;

    // Create Cold Energy References and Volumes
    int nBins = 100;
    JVec = std::vector<double>(nBins + 1);
    eCVec = std::vector<double>(nBins + 1);

    // Initialize Jmin, Jmax, and dJ
    Jmin = 0.2;                          // 5x Density
    Jmax = 2.0;                          // 0.5x Density
    dJ   = (Jmax - Jmin) / nBins;        // Volume Ratio Increment

    // Correct Jmin, Jmax s.t. J = 1 gives eCVec = 0
    int zBin = (int)((1.0 - Jmin)/dJ);
    Jmin = 1.0 - dJ * zBin;             // Corrected Minimum
    Jmax = 1.0 + dJ * (nBins - zBin);   // Corrected Maximum
    JVec[zBin] = 1.0;
    eCVec[zBin] = 0.0;

    // Intermediate Variables
    double Jp, pH, cvTH, cvT0;
    cvT0 = cv * T0;
    cvTH = 0.0;                         // Hugoniot Temperature at J = 1.0

    // Compute Cold Curve for J < 1.0 (Compression)
    for (int n = zBin-1; n > 0; n--){

        JVec[n] = JVec[n + 1] - dJ;           //Reference Volume Ratio

        Jp = 0.5 * (JVec[n + 1] + JVec[n]);   //Volume Ratio for Integration

        //Hugoniot Pressure
        pH = rho_0 * C0 * C0 * (1.0 - Jp) /
                ((1.0 - S0 *(1.0 - Jp)) * (1.0 - S0 *(1.0 - Jp)));

        //Hugoniot Temperature
        cvTH += dJ * pH / rho_0 * (1.0 - G0/2.0 * (1.0 - Jp)) *
                std::exp(-G0 * (1.0 - Jp));

        //Cold Energy
        eCVec[n] = std::exp(G0 * (1.0 - JVec[n])) * (cvTH - cvT0) + cvT0;
    }

    // Compute Cold Curve for J > 1.0 (Extension)
    cvTH = 0.0;
    for (int n = zBin+1; n < nBins+1; n++){

        JVec[n] = JVec[n - 1] + dJ;           //Reference Volume Ratio

        Jp = 0.5 * (JVec[n - 1] + JVec[n]);   //Volume Ratio for Integration

        //Hugoniot Pressure
        pH = rho_0 * C0 * C0 * (1.0 - Jp);

        //Hugoniot Temperature
        cvTH -= dJ * pH / rho_0 * (1.0 - G0/2.0 * (1.0 - Jp)) *
                std::exp(-G0 * (1.0 - Jp));

        //Cold Energy
        eCVec[n] = std::exp(G0 * (1.0 - JVec[n])) * (cvTH - cvT0) + cvT0;
    }


    // List Should Be Complete
    std::cout << "Cold Energy Curve Computed!" << std::endl;

    return;

}


double TillotsonEOSFluid::EBFromMaterialState(MaterialState& stateIN){
//Compute EB from Scalar Elastic Strains
//  Returns the breakage energy value EB computed from Material State

    MaterialTensor Ee = 0.5 * (stateIN.Be - MaterialTensor::Identity());
    MaterialTensor Ee0 = Ee.deviator();
    double eve = Ee.trace();
    double ese = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

    double A = -K * eve / 2.0;
    if (A < 0) {
        A = 0;
    }

    double EB = (stateIN.rho / rho_0) * theta * pr * (2.0 * A * A * A / (3.0 * K)
                                                      + 3.0 * G * A * ese * ese / 2.0
                                                      + G * std::sqrt(3.0 * K * G) * ese * ese * ese / 4.0);

    return EB;

}


std::vector<double> TillotsonEOSFluid::PQFromMaterialState(MaterialState& stateIN){
//Compute P,Q from Scalar Elastic Strains
//  Compute pressure and shear stress from material state

    // Output Vector
    // pBAR, qBAR, pSTAR
    std::vector<double> pq(3);
    double pBAR, qBAR, pSTAR;

    // Strain Invaraints
    MaterialTensor Ee = 0.5 * (stateIN.Be - MaterialTensor::Identity());
    MaterialTensor Ee0 = Ee.deviator();
    double eve = Ee.trace();
    double ese = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));
    double Je  = std::sqrt(stateIN.Be.det());

    // Compression Coefficient
    double A = -K * eve / 2.0;
    double dAde = -K / 2.0;
    if (A < 0) {
        A = 0; dAde = 0;
    }

    // Compute Pressure
    pBAR = (-stateIN.rho/rho_0) * (1.0 - theta * stateIN.B) * pr * dAde *
            (2.0 * A * A / K + 3.0 * G * ese * ese / 2.0);

    // Compute Shear Stress
    qBAR = (stateIN.rho/rho_0) * (1.0 - theta * stateIN.B) * pr *
            3.0 * G * (A * ese + std::sqrt(3.0 * G * K) * ese * ese / 4.0);

    // Zero Out Stresses for Positive Volumetric Strain
    if (eve > 0){
        pBAR = 0;
        qBAR = 0;
    }

    // Compute Constituent Pressure
    double J = 1.0;
    double a, eC, pH;
    if (is_compressible){
        // Use Compressible Model

        // Compute Porosity
        stateIN.phi = PorosityFromMaterialState(stateIN);

        // Compute Density Ratio
        a = std::pow(1.0 - stateIN.phi, b);
        J = Je / (a + Je * (1.0 - a));

        // Compute Cold Energy
        eC = ColdEnergyFromMaterialState(stateIN);

        // Compute Hugoniot Pressure
        pH = rho_0 * C0 * C0 * (1.0 - J);
        if (J < 1.0){
            // Compressed
            pH /= (1.0 - S0 * (1.0 - J)) * (1.0 - S0 * (1.0 - J));
        }

        // Compute Pressure
        pSTAR = pH * (1.0 - G0/2.0 * (1.0 - J)) +
                rho_0 * G0 * (eC + cv * (stateIN.Ts - T0));

        // Only Consider Positive Pressures
        if (pSTAR < 0){
            pSTAR = 0.0;
        }

    } else {
        // Use Incompressible Model
        pSTAR = 0.0;
    }

    // Construct Output
    pq[0] = pBAR;
    pq[1] = qBAR;
    pq[2] = pSTAR;

    return pq;
}

std::vector<double> TillotsonEOSFluid::PlasticFlowRulesFromMaterialState(MaterialState& stateIN){
//Compute relative rates of deformation for ev, es, and B from material state and material deformation

    // Output Vector
    std::vector<double> Rates(3);

    // Compute Material Stress State
    MaterialTensor Sy = YieldStressFromMaterialState(stateIN);

    // Compute Pressure and Shear Stress in Deformed Reference Frame
    double p = -Sy.trace()/3.0;
    MaterialTensor s0 = Sy + p * MaterialTensor::Identity();
    double q = std::sqrt(3.0 / 2.0 * s0.dot(s0));

    // Compute EB
    double EB = EBFromMaterialState(stateIN);

    // Compute Limiting Porosities
    double phi_min = phi_l * std::pow(1.0 - stateIN.B, l);
    double phi_max = phi_u * std::pow(1.0 - stateIN.B, u);

    // Compute Relative Porosity
    double tau = (phi_max - stateIN.phiP) / (phi_max - phi_min);

    // Compute Critical Porosity
    double tau_cs = std::sqrt(EB/Ec) * (1.0 - stateIN.B) / g;

    // Compute Maximum Friction Angle
    double tp = M_PI / 15.0 + std::asin(3.0 * M_0 / (6.0 + M_0));

    // Ensure Positive Pressure
    if (p <= AbsTOL){
        // Zero Pressure Yielding
        Rates[0] = 1; //evRate
        Rates[1] = 1; //esRate
        Rates[2] = 0; //BRate
    } else {

        // Dilation Angle
        double M_d = g * (tau - tau_cs) * (6.0 * std::sin(tp) / (3.0 - std::sin(tp)) - M_0);

        // Coupling Angle
        double w = M_PI / 2.0 * (1.0 - tau);

        // Plasticity Rates
        // evRate
        Rates[0] = -EB * (1.0 - stateIN.B) * (1.0 - stateIN.B) * std::sin(w) * std::sin(w) / Ec
                   * p / (p*p + q*q)
                   + q * M_d / ((M_0 + M_d)*p * (M_0 + M_d)*p);
        // esRate
        Rates[1] = EB * (1.0 - stateIN.B) * (1.0 - stateIN.B) * std::sin(w) * std::sin(w) / Ec
                   * q / (p*p + q*q)
                   + q / ((M_0 + M_d)*p * (M_0 + M_d)*p);
        // BRate
        Rates[2] = (1.0 - stateIN.B) * (1.0 - stateIN.B) * std::cos(w) * std::cos(w) / Ec;

    }

    return Rates;

}


double TillotsonEOSFluid::YieldFunctionFromMaterialState(MaterialState& stateIN){
//Compute y from material state

    // Compute Material Stress State
    MaterialTensor Sy = YieldStressFromMaterialState(stateIN);

    // Compute Pressure and Shear Stress in Deformed Reference Frame
    double p = -Sy.trace()/3.0;
    MaterialTensor s0 = Sy + p * MaterialTensor::Identity();
    double q = std::sqrt(3.0 / 2.0 * s0.dot(s0));

    // Compute EB
    double EB = EBFromMaterialState(stateIN);

    // Compute Limiting Porosities
    double phi_min = phi_l * std::pow(1.0 - stateIN.B, l);
    double phi_max = phi_u * std::pow(1.0 - stateIN.B, u);

    // Compute Relative Porosity
    double tau = (phi_max - stateIN.phiP) / (phi_max - phi_min);

    // Compute Critical Porosity
    double tau_cs = std::sqrt(EB/Ec) * (1.0 - stateIN.B) / g;

    // Compute Maximum Friction Angle
    double tp = M_PI / 15.0 + std::asin(3.0 * M_0 / (6.0 + M_0));

    // Dilation Angle
    double M_d = g * (tau - tau_cs) * (6.0 * std::sin(tp) / (3.0 - std::sin(tp)) - M_0);

    // Yield f'n
    return EB * (1.0 - stateIN.B) * (1.0 - stateIN.B) / Ec + q * q / ((M_0 + M_d)*p * (M_0 + M_d)*p) - 1;

}

MaterialTensor TillotsonEOSFluid::CauchyStressFromMaterialState(MaterialState& stateIN){
    // Compute Cauchy Stress from Material State and Deformation

    // Strain Invaraints
    MaterialTensor Ee = 0.5 * (stateIN.Be - MaterialTensor::Identity());
    MaterialTensor Be0 = stateIN.Be.deviator();
    MaterialTensor Ee0 = Ee.deviator();
    double eve = Ee.trace();
    double ese = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));
    double Je  = std::sqrt(stateIN.Be.det());

    // Compute P,Q from Helmholtz free energy
    std::vector<double> pq = PQFromMaterialState(stateIN);
    double pBAR = pq[0];
    double qBAR = pq[1];
    double pSTAR = pq[2];

    // Compute S
    MaterialTensor S = MaterialTensor();
    if (qBAR > 0){
        //Non-Zero Shear Stress
        S = qBAR / (3.0 * ese) * Be0 * stateIN.Be - pBAR * stateIN.Be;
    } else {
        //Zero Shear Stress
        S = -pBAR * stateIN.Be;
    }

    // Compute Constituent Pressure Term
    double J = 1.0;
    double a, dadp, phiSTemp;
    if (is_compressible){
        // Use Compressible Model

        // Compute Porosity
        stateIN.phi = PorosityFromMaterialState(stateIN);

        // Compute Solid Volume Fraction
        phiSTemp = 1.0 - stateIN.phi;

        // Compute Density Ratio
        a = std::pow(phiSTemp, b);
        dadp = b * a / phiSTemp;

        // Add pSTAR Constribution to S
        S -= phiSTemp * pSTAR * (a + phiSTemp * dadp * (1.0 - Je)) /
                (Je + a * (1.0 - Je) + phiSTemp * dadp * (1.0 - Je)) *
                MaterialTensor::Identity();

    } else {
        // Do Nothing
    }

    return S;
}

MaterialTensor TillotsonEOSFluid::YieldStressFromMaterialState(MaterialState& stateIN){
    // Stress Conjugate to Plastic Deformation
    // This is NOT the Cauchy Stress, but is almost the Cauchy Stress

    // Strain Invaraints
    MaterialTensor Ee = 0.5 * (stateIN.Be - MaterialTensor::Identity());
    MaterialTensor Be0 = stateIN.Be.deviator();
    MaterialTensor Ee0 = Ee.deviator();
    double eve = Ee.trace();
    double ese = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));
    double Je  = std::sqrt(stateIN.Be.det());

    // Compute P,Q from Helmholtz free energy
    std::vector<double> pq = PQFromMaterialState(stateIN);
    double pBAR = pq[0];
    double qBAR = pq[1];
    double pSTAR = pq[2];

    // Compute S
    MaterialTensor S = MaterialTensor();
    if (qBAR > 0){
        //Non-Zero Shear Stress
        S = qBAR / (3.0 * ese) * Be0 * stateIN.Be - pBAR * stateIN.Be;
    } else {
        //Zero Shear Stress
        S = -pBAR * stateIN.Be;
    }

    // Compute Constituent Pressure Term
    double a, dadp, phiSTemp;
    if (is_compressible){
        // Use Compressible Model

        // Compute Porosity
        stateIN.phi = PorosityFromMaterialState(stateIN);

        // Compute Solid Volume Fraction
        phiSTemp = 1.0 - stateIN.phi;

        // Compute Density Ratio
        a = std::pow(phiSTemp, b);
        dadp = b * a / phiSTemp;

        // Add pSTAR Constribution to S
        S -= phiSTemp * pSTAR * a /
             (Je + a * (1.0 - Je) + phiSTemp * dadp * (1.0 - Je)) *
             MaterialTensor::Identity();

    } else {
        // Do Nothing
    }

    return S;
}

double TillotsonEOSFluid::PorosityFromMaterialState(MaterialState &stateIN) {
    // Estimate Porosity from Material Deformation

    // Elastic Volume Ratio
    double Je  = std::sqrt(stateIN.Be.det());

    // Compute Porosity Using Bisection
    double J = 1.0;
    double a, phiSTemp;
    double phiA, phiB, phiC;
    double yA, yB, yC;

    // Bisection
    int n = 0;
    double A = (1.0/Je - 1.0);
    phiA = 0; yA = -1;
    phiB = 1; yB =  1;
    phiC = 1; yC =  1;
    if (is_compressible){
        // Use Compressible Model

        while (std::abs(yC) > AbsTOL && n < MaxIter){
            //increment n
            n++;

            //set phiC to average of phiA and phiB
            phiC = 0.5 * (phiA + phiB);

            //evaluate yC
            yC = phiC + A * std::pow(phiC, b+1) - stateIN.rho / rho_0;

            // assign phiA, phiB
            if (yC * yA > 0){
                phiA = phiC; yA = yC;
            } else {
                phiB = phiC; yB = yC;
            }
        }

        // Compute Solid Volume Fraction
        phiSTemp = phiC;
    } else {
        // Solid is Incompressible
        phiSTemp = stateIN.rho / rho_0;
    }

    // Porosity = 1 - Solid Volume Fraction
    stateIN.phi = 1.0 - phiSTemp;
    return 1.0 - phiSTemp;
}


double TillotsonEOSFluid::ContactEnergyFromMaterialState(MaterialState &stateIN) {
    // Compute Contact Strain Energy (J/kg) from Material Deformation

    // Strain Invaraints
    MaterialTensor Ee = 0.5 * (stateIN.Be - MaterialTensor::Identity());
    MaterialTensor Be0 = stateIN.Be.deviator();
    MaterialTensor Ee0 = Ee.deviator();
    double eve = Ee.trace();
    double ese = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

    double A = -K * eve / 2.0;
    if (A < 0){
        A = 0;
    }

    double psiC = (1.0 - theta * stateIN.B) / rho_0 * pr *
            (2.0 * A * A * A / (3.0 * K) +
            3.0 * G * A * ese * ese / 2.0 +
            G * std::sqrt(3.0 * G * K) * ese * ese * ese / 4.0);

    return psiC;
}


double TillotsonEOSFluid::ColdEnergyFromMaterialState(MaterialState &stateIN) {
    // Compute Cold Energy (J/kg) from Material State
    double eC = 0;

    // Elastic Volume Ratio
    double Je  = std::sqrt(stateIN.Be.det());

    // Porosity
    stateIN.phi = PorosityFromMaterialState(stateIN);

    // Solid Volume Fraction
    double phiSTemp = 1.0 - stateIN.phi;

    // Cold Energy Computation
    double a, J;
    double cvTH, cvT0, dJp, Jp, pHTemp;
    double Jminus, Jplus, eCminus, eCplus;
    int N = 25;
    int jBin;
    if (is_compressible){
        // Use Compressible Model

        // Constituent Solid Volume Ratio
        a = std::pow(phiSTemp, b);
        J = Je / (a + Je * (1.0 - a));

        // Look Up or Compute eC
        if (J > Jmin && J < Jmax){
            // Use Look Up Table

            // Identify Bin for J
            jBin = (int) ((J - Jmin) / dJ);

            // Identify Bin Variables
            Jminus  = JVec[jBin];
            Jplus   = JVec[jBin+1];
            eCminus = eCVec[jBin];
            eCplus  = eCVec[jBin+1];

            // Linear Interpolation w/in Bin
            eC = eCminus + (J - Jminus) * (eCplus - eCminus) / (Jplus - Jminus);

        } else {
            // Compute Numerically
            cvTH = 0; cvT0 = cv * T0;
            for (int i=0; i<N; i++){
                dJp = (1 - J)/N;
                Jp = 1.0 - (i + 0.5) * dJp;

                pHTemp = rho_0 * C0 * C0 * (1.0 - Jp);
                if (Jp < 1.0){
                    pHTemp /= ((1.0 - S0 *(1.0 - Jp)) * (1.0 - S0 *(1.0 - Jp)));
                }

                cvTH += dJp * pHTemp / rho_0 * (1.0 - G0/2.0 * (1.0 - Jp)) *
                        std::exp(-G0 * (1.0 - Jp));
            }
            eC = std::exp(G0 * (1.0 - J)) * (cvTH - cvT0) + cvT0;
        }

    } else {
        // Use Incompressible Model
        // Do Nothing
    }

    return eC;
}


double TillotsonEOSFluid::TemperatureFromMaterialState(MaterialState &stateIN) {
    // Compute Temperature from Material State

    // Compute Contact Strain Energy
    double psiC = ContactEnergyFromMaterialState(stateIN);

    // Compute Cold Strain Energy
    double eC = ColdEnergyFromMaterialState(stateIN);

    // Compute Temperature or Energy
    if (is_adiabatic){
        // Use Thermal Version of Model
        stateIN.Ts = (stateIN.Es - psiC - eC)/cv;
        return stateIN.Ts;

    } else {
        // Use Isothermal Version of Model
        stateIN.Es = psiC + eC + cv * T0;
        stateIN.Ts = T0;
        return T0;
    }

}
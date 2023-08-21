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
void CompressibleDamageMechanicsSandstone::init(Job* job, Body* body){
    if (fp64_props.size() < 25){
        std::cout << fp64_props.size() << " != 25\n";
        fprintf(stderr, "%s:%s:", __FILE__, __func__);

        std::cerr << "Need at least 24 properties defined:\n";
        std::cerr << "    pr     [Pa]     Nonlinear Reference Pressure\n";
        std::cerr << "    Kg     [--]     Nonlinear Granular Bulk Modulus Coefficient\n";
        std::cerr << "    Gg     [--]     Nonlinear Granular Shear Modulus Coefficient\n";
        std::cerr << "    Kd     [Pa]     Linear Porous Rock Bulk Modulus\n";
        std::cerr << "    Gd     [Pa]     Linear Porous Rock Shear Modulus\n";
        std::cerr << "    EBc    [J/m^3]  Critical Breakage Energy\n";
        std::cerr << "    EDc    [J/m^3]  Critical Damage Energy\n";
        std::cerr << "    c      [Pa]     Cohesion\n";
        std::cerr << "    M_0    [--]     Friction Angle\n";
        std::cerr << "    g      [--]     Dilation Factor (0 to 1)\n";
        std::cerr << "    phi_l  [--]     Lower Porosity Limit at (B=0)\n";
        std::cerr << "    phi_u  [--]     Upper Porosity Limit at (B=0)\n";
        std::cerr << "    l      [--]     Lower Porosity Power Law Rate\n";
        std::cerr << "    u      [--]     Upper Porosity Power Law Rate\n";
        std::cerr << "    theta  [--]     Grading Index\n";
        std::cerr << "    rho_0  [kg/m^3] Reference Density of Constituent Solid\n";
        std::cerr << "    C0     [m/s]    Bulk Sound Speed at Reference Density\n";
        std::cerr << "    S0     [--]     Slope of Us--Up Curve\n";
        std::cerr << "    G0     [--]     Gruneisen Parameter\n";
        std::cerr << "    cv     [J/kg*K] Heat Capacity at Constant Volume\n";
        std::cerr << "    b      [--]     Porosity Function Power\n";
        std::cerr << "    T0     [K]      Reference Temperature\n";
        std::cerr << "    B_0    [--]     Initial Value of Breakage\n";
        std::cerr << "    D_0    [--]     Initial Value of Damage\n";
        std::cerr << "    phi_0  [--]     Initial Porosity\n";

        exit(0);
    } else {
        pr      = fp64_props[0];
        Kg      = fp64_props[1];
        Gg      = fp64_props[2];
        Kd      = fp64_props[3];
        Gd      = fp64_props[4];
        EBc     = fp64_props[5];
        EDc     = fp64_props[6];
        c       = fp64_props[7];
        M_0     = fp64_props[8];
        g       = fp64_props[9];
        phi_l   = fp64_props[10];
        phi_u   = fp64_props[11];
        l       = fp64_props[12];
        u       = fp64_props[13];
        theta   = fp64_props[14];

        rho_0   = fp64_props[15];
        C0      = fp64_props[16];
        S0      = fp64_props[17];
        G0      = fp64_props[18];
        cv      = fp64_props[19];
        b       = fp64_props[20];
        T0      = fp64_props[21];

        double B_0 = fp64_props[22];
        double D_0 = fp64_props[23];
        double phi_0 = fp64_props[24];

        std::cout << "Material Properties:\n";
        std::cout << "    pr     = " << pr << " [Pa]\n";
        std::cout << "    Kg     = " << Kg << " [--]\n";
        std::cout << "    Gg     = " << Gg << " [--]\n";
        std::cout << "    Kd     = " << Kd << " [Pa]\n";
        std::cout << "    Gd     = " << Gd << " [Pa]\n";
        std::cout << "    EBc    = " << EBc << " [J/m^3]\n";
        std::cout << "    EDc    = " << EDc << " [J/m^3]\n";
        std::cout << "    c      = " << c << " [Pa]\n";
        std::cout << "    M_0    = " << M_0 << " [--]\n";
        std::cout << "    g      = " << g << " [--]\n";
        std::cout << "    phi_l  = " << phi_l << " [--]\n";
        std::cout << "    phi_u  = " << phi_u << " [--]\n";
        std::cout << "    l      = " << l << " [--]\n";
        std::cout << "    u      = " << u << " [--]\n";
        std::cout << "    theta  = " << theta << " [--]\n";
        std::cout << "    rho_0  = " << rho_0 << " [kg/m^3]\n";
        std::cout << "    C0     = " << C0 << " [m/s]\n";
        std::cout << "    S0     = " << S0 << " [--]\n";
        std::cout << "    G0     = " << G0 << " [--]\n";
        std::cout << "    cv     = " << cv << " [J/kg*K]\n";
        std::cout << "    b      = " << b << " [--]\n";
        std::cout << "    T0     = " << T0 << " [K]\n";
        std::cout << "    B_0    = " << B_0 << " [--]\n";
        std::cout << "    D_0    = " << D_0 << " [--]\n";
        std::cout << "    phi_0  = " << phi_0 << " [--]\n";

        //check for string flags
        for (int s=0; s<str_props.size(); s++){
            if (str_props[s].compare("ISOTHERMAL") == 0){
                //isothermal simulation
                is_adiabatic = false;
            } else if (str_props[s].compare("INCOMPRESSIBLE") == 0){
                //incompressible model
                is_compressible = false;
            } else if (str_props[s].compare("USE_ARTIFICIAL_VISCOSITY") == 0){
                //artificial velocity
                use_artificial_viscosity = true;
            } else if (str_props[s].compare("USE_NEWTONS_METHOD") == 0){
                //newtons method for porosity calculation
                use_newtons_method = true;
            }
        }

        if (use_artificial_viscosity && fp64_props.size() > 20){
            h_i = fp64_props[20];
            std::cout << "  Using Artificial Viscosity:\n";
            std::cout << "    h_i    = " << h_i << " [m]\n";
        } else if (use_artificial_viscosity){
            std::cout << "  *NOT* Using Artificial Viscosity! No h_i parameter defined.\n";
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

        if (use_newtons_method){
            std::cout << "  Using Newton's Method for Porosity Computation.\n";
        } else {
            std::cout << "  Using Bisection Method for Porosity Computation.\n";
        }

        //initialize B, phiS, phiP, Be, F
        B = B_0 * Eigen::VectorXd::Ones(body->points->x.size());
        D = D_0 * Eigen::VectorXd::Ones(body->points->x.size());
        phiS = (1 - phi_0) * Eigen::VectorXd::Ones(body->points->x.size());
        phiP = phi_0 * Eigen::VectorXd::Ones(body->points->x.size());
        Es = cv * T0 * Eigen::VectorXd::Ones(body->points->x.size());
        Ts = T0 * Eigen::VectorXd::Ones(body->points->x.size());
        Be = MaterialTensorArray(body->points->x.size());
        F = MaterialTensorArray(body->points->x.size());

        //initialize scalar outputs
        evDot = Eigen::VectorXd(body->points->x.size());
        esDot = Eigen::VectorXd(body->points->x.size());
        BDot = Eigen::VectorXd(body->points->x.size());
        DDot = Eigen::VectorXd(body->points->x.size());

        //debug variables
        kVec = Eigen::VectorXd(body->points->x.size());
        y1Vec = Eigen::VectorXd(body->points->x.size());
        y2Vec = Eigen::VectorXd(body->points->x.size());
        y3Vec = Eigen::VectorXd(body->points->x.size());
        lVec = Eigen::VectorXd(body->points->x.size());

        //fill Be, F with initial deformation
        for (int i=0; i<body->points->x.size(); i++){
            Be(i) = MaterialTensor::Identity();
            F(i) = MaterialTensor::Identity();
        }

        //compute cold energy reference states
        ComputeColdEnergyReferenceStates();
    }

    std::cout << "Material Initialized: [" << body->name << "]." << std::endl;

    return;
}

/*----------------------------------------------------------------------------*/
//calculate stress state given prior state and strain-rate
void CompressibleDamageMechanicsSandstone::calculateStress(Job* job, Body* body, int SPEC){
    //convergence criteria
    double AbsTOL = 1e-7;
    int MaxIter = 50;

    //tensors for matrix calculation
    MaterialTensor S, s0, Sn, L, DTensor, W;
    MaterialTensor Ben;
    MaterialTensor Ee, Ee0, Be0, BeRate, BeRate0;

    //state for matrix calculations
    MaterialState mat_state, mat_stateA, mat_stateB, mat_stateC;

    //scalars for matrix calculation
    std::vector<double> Rates, ys;
    double p, q, EB, ED, evRate, esRate, BRate, DRate, phi_max;
    double lambda, lambdaA, lambdaB, lambdaC, lambdaMAX, lambdaMAXEST;
    double y, y0, yA, yB, yC;
    double y1, y2, y3;
    double trDe, a, dadp, phiSTemp, JeTemp, detBe;

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
        mat_state.D    = D(i);
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
        DTensor = 0.5*(L+L.transpose());
        W = 0.5*(L-L.transpose());

        // Deformation at n
        Ben = Be[i];

        // Estimated Deformation w/o Plasticity
        mat_state.Be = Be[i] + job->dt * (L * Ben + Ben * L.transpose());
        mat_state.F  = F[i] + job->dt * L * F[i];

        // Estimate Porosity
        mat_state.phi = PorosityFromMaterialState(mat_state);

        // Compute Yield Function Values
        ys = YieldFunctionsFromMaterialState(mat_state);
        y1 = ys[0]; y2 = ys[1]; y3 = ys[2];

        // Check for Yielding in Tension
        if (y2 >= 0){
            // Set Flag
            tensileYield = true;

            // No Deviatoric Stresses
            mat_state.Be = mat_state.Be.trace() / 3.0 * MaterialTensor::Identity();
        }

        // Check for Yielding in Compaction
        if (y3 >= 0){
            // No Stresses in Under-Compacted State
            compactionYield = true;
        }


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
        DRate = Rates[2];
        BRate = Rates[3];

        // Assign Plasticity Multiplier
        lambda = 0;

        // Estimate Maximum Plastic Step Size
        BeRate = MaterialTensor();
        if (q > 0) {
            BeRate = 3.0 / 2.0 * esRate / q * (s0 * Ben + Ben * s0) + 2.0 / 3.0 * evRate * Ben;
        } else {
            BeRate = 2.0 / 3.0 * evRate * Ben;
        }

        lambdaMAX = 1e10; //a ludicrous number
        lambdaMAXEST = lambdaMAX;

        // Check for Plastic Compaction/Dilation
        if (BeRate.trace() < 0 && mat_state.Be.trace() < 3) {
            lambdaMAX = (mat_state.Be.trace() - 3) / BeRate.trace();
        } else if (BeRate.trace() > 0) {
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

        // Check for Over-Damage
        if (DRate > 0 && mat_state.D < 1){
            lambdaMAXEST = (1.0 - mat_state.D) / DRate;
        }
        if (lambdaMAXEST < lambdaMAX){
            lambdaMAX = lambdaMAXEST;
        }

        // Check for Porosity Estimate
        if (evRate > 0){
            lambdaMAXEST = 1.0 / evRate;
        } else if (evRate < 0) {
            lambdaMAXEST = -mat_state.phiP / (evRate * (1.0 - mat_state.phiP));
        }
        if (lambdaMAXEST < lambdaMAX){
            lambdaMAX = lambdaMAXEST;
        }

        // Check for Valid Limits
        /*
        if (lambdaMAX < 0){
            std::cout << "Hmmm...\n"
                      << "    " << "lambdaMAX = " << lambdaMAX << "\n"
                      << "    " << "D = " << mat_state.D << "\n"
                      << "    " << "B = " << mat_state.B << "\n"
                      << "    " << "phiP = " << mat_state.phiP << "\n"
                      << "    " << "evRate = " << evRate << "\n"
                      << "    " << "esRate = " << esRate << "\n"
                      << "    " << "DRate = " << DRate << "\n"
                      << "    " << "BRate = " << BRate << "\n"
                      << "    " << "Be.trace() = " << mat_state.Be.trace() << "\n"
                      << "    " << Ben(0,0) << ", " << Ben(0,1) << ", " << Ben(0,2) << ", "
                      << "    " << Ben(1,0) << ", " << Ben(1,1) << ", " << Ben(1,2) << ", "
                      << "    " << Ben(2,0) << ", " << Ben(2,1) << ", " << Ben(2,2) << "\n"
                      << "    " << "i = " << i << std::endl;
        }
         */

        // Update Material State
        if (tensileYield){
            // Material has Yielded in Tension

            // No Deviatoric Stresses
            mat_state.Be = mat_state.Be.trace()/3.0 * MaterialTensor::Identity();
            BeRate = BeRate.trace()/3.0 * MaterialTensor::Identity();
            BRate = 0;

            // Use Bisection Method
            k = 0; lambda = 0;
            lambdaA = lambda;
            lambdaB = lambdaMAX;
            lambdaC = lambdaB;
            yA = 1;
            yB = 10;
            yC = 100;

            // Compute Material State and Yield Function for A
            mat_stateA      = mat_state;
            mat_stateA.D    = mat_stateA.D      + lambdaA * DRate;
            mat_stateA.B    = mat_stateA.B      + lambdaA * BRate;
            mat_stateA.phiP = mat_stateA.phiP   + lambdaA * evRate * (1.0 - mat_stateA.phiP);
            mat_stateA.Be   = mat_stateA.Be     - lambdaA * BeRate;
            mat_stateA.phi  = PorosityFromMaterialState(mat_stateA);


            ys = YieldFunctionsFromMaterialState(mat_stateA);
            yA = ys[1]; //y2

            // Compute Material State and Yield Function for B
            mat_stateB      = mat_state;
            mat_stateB.D    = mat_stateB.D      + lambdaB * DRate;
            mat_stateB.B    = mat_stateB.B      + lambdaB * BRate;
            mat_stateB.phiP = mat_stateB.phiP   + lambdaB * evRate * (1.0 - mat_stateB.phiP);
            mat_stateB.Be   = mat_stateB.Be     - lambdaB * BeRate;
            mat_stateB.phi  = PorosityFromMaterialState(mat_stateB);

            ys = YieldFunctionsFromMaterialState(mat_stateB);
            yB = ys[1]; //y2

            // Find Plastic Step Length that Solves y = 0
            while (std::abs(yA) > AbsTOL && std::abs(yB) > AbsTOL && k < MaxIter) {

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
                mat_stateC.D    = mat_stateC.D      + lambdaC * DRate;
                mat_stateC.B    = mat_stateC.B      + lambdaC * BRate;
                mat_stateC.phiP = mat_stateC.phiP   + lambdaC * evRate * (1.0 - mat_stateC.phiP);
                mat_stateC.Be   = mat_stateC.Be     - lambdaC * BeRate;
                mat_stateC.phi  = PorosityFromMaterialState(mat_stateC);

                ys = YieldFunctionsFromMaterialState(mat_stateC);
                yC = ys[1]; //y2

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

            // Final State Update
            mat_state.D    = mat_state.D      + lambda * DRate;
            mat_state.B    = mat_state.B      + lambda * BRate;
            mat_state.phiP = mat_state.phiP   + lambda * evRate * (1.0 - mat_state.phiP);
            mat_state.Be   = mat_state.Be     - lambda * BeRate;
            mat_state.phi  = PorosityFromMaterialState(mat_state);
            mat_state.S    = CauchyStressFromMaterialState(mat_state);

            // Compute Elastic Volumetric Strain-Rate
            trDe = DTensor.trace() - lambda * evRate;

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
            mat_state.Es += 0.5 * (Sn.dot(DTensor) + mat_state.S.dot(DTensor)) * job->dt / mat_state.rho;
            mat_state.Ts  = TemperatureFromMaterialState(mat_state);

            ys = YieldFunctionsFromMaterialState(mat_state);
            y1 = ys[0]; y2 = ys[1]; y3 = ys[2];

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

                //body->points->T[i] += rho_0 * C0 * h_i * DTensor.trace() * MaterialTensor::Identity();
            }

            // Porosity and Breakage
            D(i)    = mat_state.D;
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
            DDot(i)     = lambda * DRate / job->dt;
            y1Vec(i)     = y1;
            y2Vec(i)     = y2;
            y3Vec(i)     = y3;
            kVec(i)     = k;
            lVec(i)     = lambda;

        } else if (compactionYield) {
            // Material has Yielded in Under-Compacted State
            // - no stress
            body->points->T[i].setZero();
            // - zero stress elastic strain
            Be[i] = MaterialTensor::Identity();
            mat_state.Be = Be[i];
            // - porosity evolution continues
            phiP(i) += job->dt * DTensor.trace() * (1 - phiP(i));
            phiS(i) = 1.0 - PorosityFromMaterialState(mat_state);
            // - energy is only thermal
            Es(i) = cv * Ts(i);
            // - save deformation
            F[i]       = mat_state.F;

            // Set Scalar Outputs and History Variables
            evDot(i) = DTensor.trace();
            esDot(i) = std::sqrt(2.0 / 3.0 * DTensor.deviator().dot(DTensor.deviator()));
            BDot(i) = 0;
            DDot(i) = 0;
            y1Vec(i) = y1;
            y2Vec(i) = y2;
            y3Vec(i) = 0;
            kVec(i) = -2;
            lVec(i) = 0;

        } else {
            // Material is NOT in Tension or Under-Compacted

            // Compute Yield Function Value
            ys = YieldFunctionsFromMaterialState(mat_state);
            y0 = ys[0]; //y1

            // Initial Guess for Plastic Deformation
            lambda = 0;

            // Initialize Counter
            k = 0;

            // Check if Material has Yielded
            if (y0 > 0){

                // Use Secant Method
                lambdaA = lambda;               //Plastic Step Lengths
                lambdaB = 1e-2 * lambdaMAX;
                lambdaC = lambdaB;
                yA = 1;                         //Yield Function Values
                yB = 10;
                yC = 100;

                // Compute Material State and Yield Function for A
                mat_stateA      = mat_state;
                mat_stateA.D    = mat_stateA.D      + lambdaA * DRate;
                mat_stateA.B    = mat_stateA.B      + lambdaA * BRate;
                mat_stateA.phiP = mat_stateA.phiP   + lambdaA * evRate * (1.0 - mat_stateA.phiP);
                mat_stateA.Be   = mat_stateA.Be     - lambdaA * BeRate;
                mat_stateA.phi  = PorosityFromMaterialState(mat_stateA);

                ys = YieldFunctionsFromMaterialState(mat_stateA);
                yA = ys[0]; //y1

                // Compute Material State and Yield Function for B
                mat_stateB      = mat_state;
                mat_stateB.D    = mat_stateB.D      + lambdaB * DRate;
                mat_stateB.B    = mat_stateB.B      + lambdaB * BRate;
                mat_stateB.phiP = mat_stateB.phiP   + lambdaB * evRate * (1.0 - mat_stateB.phiP);
                mat_stateB.Be   = mat_stateB.Be     - lambdaB * BeRate;
                mat_stateB.phi  = PorosityFromMaterialState(mat_stateB);

                ys = YieldFunctionsFromMaterialState(mat_stateB);
                yB = ys[0]; //y1

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
                        mat_stateB.D    = mat_stateB.D      + lambdaB * DRate;
                        mat_stateB.B    = mat_stateB.B      + lambdaB * BRate;
                        mat_stateB.phiP = mat_stateB.phiP   + lambdaB * evRate * (1.0 - mat_stateB.phiP);
                        mat_stateB.Be   = mat_stateB.Be     - lambdaB * BeRate;
                        mat_stateB.phi  = PorosityFromMaterialState(mat_stateB);

                        ys = YieldFunctionsFromMaterialState(mat_stateB);
                        yB = ys[0]; //y1

                    }

                    // Ensure Yield Function Values Differ
                    if (std::abs(yA) > AbsTOL && std::abs(yB) > AbsTOL && std::abs(yA - yB) < AbsTOL){
                        // Increment Counter
                        k++;

                        // Double lambdaB
                        lambdaB *= 2.0;

                        // Set Validity Flag
                        validB = false;

                        // Recompute Material State and Yield Function for B
                        mat_stateB      = mat_state;
                        mat_stateB.D    = mat_stateB.D      + lambdaB * DRate;
                        mat_stateB.B    = mat_stateB.B      + lambdaB * BRate;
                        mat_stateB.phiP = mat_stateB.phiP   + lambdaB * evRate * (1.0 - mat_stateB.phiP);
                        mat_stateB.Be   = mat_stateB.Be     - lambdaB * BeRate;
                        mat_stateB.phi  = PorosityFromMaterialState(mat_stateB);

                        ys = YieldFunctionsFromMaterialState(mat_stateB);
                        yB = ys[0]; //y1

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
                    mat_stateC.D    = mat_stateC.D      + lambdaC * DRate;
                    mat_stateC.B    = mat_stateC.B      + lambdaC * BRate;
                    mat_stateC.phiP = mat_stateC.phiP   + lambdaC * evRate * (1.0 - mat_stateC.phiP);
                    mat_stateC.Be   = mat_stateC.Be     - lambdaC * BeRate;
                    mat_stateC.phi  = PorosityFromMaterialState(mat_stateC);

                    ys = YieldFunctionsFromMaterialState(mat_stateC);
                    yC = ys[0]; //y1


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

                // Stabilize Search:
                if ((std::abs(y0) < std::abs(y))){
                    // Don't Take Inelastic Step if Yield F'n Value WORSE!
                    lambda = 0;
                    y = y0;
                }

            } else {
                // No Plastic Deformation
            }

            // Final State Update
            mat_state.B    = mat_state.B      + lambda * BRate;
            mat_state.D    = mat_state.D      + lambda * DRate;
            mat_state.phiP = mat_state.phiP   + lambda * evRate * (1.0 - mat_state.phiP);
            mat_state.Be   = mat_state.Be     - lambda * BeRate;
            mat_state.phi  = PorosityFromMaterialState(mat_state);
            mat_state.S    = CauchyStressFromMaterialState(mat_state);

            // Compute Elastic Volumetric Strain-Rate
            trDe = DTensor.trace() - lambda * evRate;

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
            mat_state.Es += 0.5 * (Sn.dot(DTensor) + mat_state.S.dot(DTensor)) * job->dt / mat_state.rho;
            mat_state.Ts  = TemperatureFromMaterialState(mat_state);

            ys = YieldFunctionsFromMaterialState(mat_state);

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

                //body->points->T[i] += rho_0 * C0 * h_i * DTensor.trace() * MaterialTensor::Identity();
            }

            // Porosity and Breakage
            D(i)    = mat_state.D;
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
            DDot(i)     = lambda * DRate / job->dt;
            BDot(i)     = lambda * BRate / job->dt;
            y1Vec(i)     = ys[0];
            y2Vec(i)     = ys[1];
            y3Vec(i)     = ys[2];
            kVec(i)     = k;
            lVec(i)     = lambda;
        }

        // Check for Maximum Iterations
        /*
        if (k >= MaxIter && std::abs(yA) > 1e3 * AbsTOL){
            std::cout << "MaxIter Reached: yA = " << yA
                      << ", yB = " << yB
                      << ", lambdaA = " << lambdaA
                      << ", lambdaB = " << lambdaB
                      << ", lambdaMAX = " << lambdaMAX
                      << ", DRate = " << DRate
                      << ", BRate = " << BRate
                      << ", evRate = " << evRate
                      << ", esRate = " << esRate
                      << ", i = " << i << std::endl;
        }
         */


    }

    return;
}


/*----------------------------------------------------------------------------*/
//define stress assignement for consistency with history dependent materials
void CompressibleDamageMechanicsSandstone::assignStress(Job* job, Body* body, MaterialTensor& stressIN, int idIN, int SPEC){
    std::cout << "WARNING: CompressibleDamageMechanicsSandstone does not implement assignStress()." << std::endl;
    return;
}


/*----------------------------------------------------------------------------*/
//define pressure assignement for consistency with history dependent materials
void CompressibleDamageMechanicsSandstone::assignPressure(Job* job, Body* body, double pressureIN, int idIN, int SPEC){
    std::cout << "WARNING: CompressibleDamageMechanicsSandstone does not implement assignPressure()." << std::endl;
    return;
}


/*----------------------------------------------------------------------------*/
//frame and state writing
void CompressibleDamageMechanicsSandstone::writeFrame(Job* job, Body* body, Serializer* serializer){

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
    serializer->writeScalarArray(D, "D");
    serializer->writeScalarArray(B, "B");
    serializer->writeScalarArray(Es, "Es");
    serializer->writeScalarArray(Ts, "Ts");
    serializer->writeScalarArray(evDot, "evDot");
    serializer->writeScalarArray(esDot, "esDot");
    serializer->writeScalarArray(DDot, "DDot");
    serializer->writeScalarArray(BDot, "BDot");
    serializer->writeScalarArray(y1Vec, "y1Vec");
    serializer->writeScalarArray(y2Vec, "y2Vec");
    serializer->writeScalarArray(y3Vec, "y3Vec");
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

std::string CompressibleDamageMechanicsSandstone::saveState(Job* job, Body* body, Serializer* serializer, std::string filepath){
    return "err";
}

int CompressibleDamageMechanicsSandstone::loadState(Job* job, Body* body, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//Helper Functions

void CompressibleDamageMechanicsSandstone::ComputeColdEnergyReferenceStates(){
//Compute Reference Values for Cold Energy
//  Compute discrete reference cold energy values between Jmin and Jmax to save computation time during simulation

    // State What You Are Doing
    std::cout << "Computing cold energy curve for CompressibleDamageMechanicsSandstone..." << std::endl;

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


std::vector<double> CompressibleDamageMechanicsSandstone::EDEBFromMaterialState(MaterialState& stateIN){
//Compute EB from Scalar Elastic Strains
//  Returns the breakage energy value EB computed from Material State

    MaterialTensor Ee = 0.5 * (stateIN.Be - MaterialTensor::Identity());
    MaterialTensor Ee0 = Ee.deviator();
    double eve = Ee.trace();
    double ese = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

    // Compute Specific Porous Rock Strain Energy
    double psiD = (1.0 / rho_0) * (Kd * eve * eve / 2.0 + 3.0 * Gd * ese * ese / 2.0);

    // Compute Specific Granular Strain Energy
    if (eve > 0){
        eve = 0;
    }
    double psiG = (1.0/ rho_0) * (1.0 - theta * stateIN.B) * pr * (-Kg * Kg * eve * eve * eve / 12.0
                                                                   - 3.0 * Gg * Kg * eve * ese * ese / 4.0
                                                                   + Gg * std::sqrt(3.0 * Kg * Gg) * ese * ese * ese / 4.0);

    // Compute Damage and Breakage Energy Conjugates
    std::vector<double>EDEB = std::vector<double>(2);
    EDEB[0] = stateIN.rho * (psiD - psiG);

    EDEB[1] = stateIN.D * (stateIN.rho / rho_0) * theta * pr * (-Kg * Kg * eve * eve * eve / 12.0
                                                                   - 3.0 * Gg * Kg * eve * ese * ese / 4.0
                                                                   + Gg * std::sqrt(3.0 * Kg * Gg) * ese * ese * ese / 4.0);

    return EDEB;

}


std::vector<double> CompressibleDamageMechanicsSandstone::PQFromMaterialState(MaterialState& stateIN){
//Compute P,Q from Scalar Elastic Strains
//  Compute pressure and shear stress from material state

    // Output Vector
    // pBARD, qBARD, pBARG, qBARG, pSTAR
    std::vector<double> pq(5);
    double pBARD, qBARD, pBARG, qBARG, pSTAR;

    // Strain Invaraints
    MaterialTensor Ee = 0.5 * (stateIN.Be - MaterialTensor::Identity());
    MaterialTensor Ee0 = Ee.deviator();
    double eve = Ee.trace();
    double ese = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));
    double Je  = std::sqrt(stateIN.Be.det());

    // Compute Porous Rock Pressure
    pBARD = (-stateIN.rho/rho_0) * (Kd * eve);

    // Compute Porous Rock Shear Strength
    qBARD = (stateIN.rho/rho_0) * (3.0 * Gd * ese);

    // Compression Coefficient
    double A = -Kg * eve / 2.0;
    double dAde = -Kg / 2.0;
    if (A < 0) {
        A = 0; dAde = 0;
    }

    // Compute Granular Pressure
    pBARG = (-stateIN.rho/rho_0) * (1.0 - theta * stateIN.B) * pr * dAde *
            (2.0 * A * A / Kg + 3.0 * Gg * ese * ese / 2.0);

    // Compute Granular Shear Stress
    qBARG = (stateIN.rho/rho_0) * (1.0 - theta * stateIN.B) * pr *
            3.0 * Gg * (A * ese + std::sqrt(3.0 * Gg * Kg) * ese * ese / 4.0);

    // Compute Constituent Pressure
    double J = 1.0;
    double a, eC, pH;
    if (is_compressible && std::isfinite(b)){
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

    } else {
        // Use Incompressible Model
        pSTAR = 0.0;
    }

    // Construct Output
    pq[0] = pBARD;
    pq[1] = qBARD;
    pq[2] = pBARG;
    pq[3] = qBARG;
    pq[4] = pSTAR;

    return pq;
}

std::vector<double> CompressibleDamageMechanicsSandstone::PlasticFlowRulesFromMaterialState(MaterialState& stateIN){
//Compute relative rates of deformation for ev, es, and B from material state and material deformation

    // Output Vector
    std::vector<double> Rates(4);

    // Compute Material Stress State
    MaterialTensor Sy = YieldStressFromMaterialState(stateIN);

    // Compute Pressure and Shear Stress in Deformed Reference Frame
    double p = -Sy.trace()/3.0;
    MaterialTensor s0 = Sy + p * MaterialTensor::Identity();
    double q = std::sqrt(3.0 / 2.0 * s0.dot(s0));

    // Compute ED,EB
    std::vector<double> EDEB = EDEBFromMaterialState(stateIN);
    double ED = EDEB[0];
    double EB = EDEB[1];
    if (ED < 0){
        ED = 0;
    }

    // Compute Limiting Porosities
    double phi_min = phi_l * std::pow(1.0 - stateIN.B, l);
    double phi_max = phi_u * std::pow(1.0 - stateIN.B, u);

    // Compute Relative Porosity
    double tau = (phi_max - stateIN.phiP) / (phi_max - phi_min);

    // Compute Critical Porosity
    double tau_cs = std::sqrt(EB/EBc) * (1.0 - stateIN.B) / g;

    // Compute Maximum Friction Angle
    double tp = M_PI / 15.0 + std::asin(3.0 * M_0 / (6.0 + M_0));

    // Compute Yield Function Values
    std::vector<double> y = YieldFunctionsFromMaterialState(stateIN);
    double y1 = y[0];
    double y2 = y[1];
    double y3 = y[2];

    // Check for Free Compression
    if (y3 >= 0){

        // No Damage/Breakage/Shear Evolution (No Shear Stress Either)
        // Zero Pressure Yielding
        Rates[0] = -1; //evRate
        Rates[1] = 0; //esRate
        Rates[2] = 0; //DRate
        Rates[3] = 0; //BRate

    // Check for Tensile Failure
    } else if (y2 >= 0) {

        // Damage OR Volumetric Evolution (No Shear Stress Either)
        // Check if Material is Still Intact
        if (stateIN.D < 1) {
            Rates[0] = 0; //evRate
            Rates[1] = 0; //esRate
            Rates[2] = 1; //DRate
            Rates[3] = 0; //BRate
        } else {
            Rates[0] = 1; //evRate
            Rates[1] = 0; //esRate
            Rates[2] = 0; //DRate
            Rates[3] = 0; //BRate
        }

    // Check for Compressive Yielding
    } else if (y1 >= 0){

        //Damage/Breakage/Shear/Volumetric Evolution

        // Dilation Angle
        double M_d = g * (tau - tau_cs) * (6.0 * std::sin(tp) / (3.0 - std::sin(tp)) - M_0);

        // Coupling Angle
        double w = M_PI / 2.0 * (1.0 - tau);
        if (tau > 1.0){
            w = 0.0;
        } else if (tau < 0.0){
            w = M_PI / 2.0;
        }

        // Plasticity Rates
        // evRate
        Rates[0] = -EB * (1.0 - stateIN.B) * (1.0 - stateIN.B) * std::sin(w) * std::sin(w) / EBc
                   * p / (p*p + q*q)
                   + q * M_d / ((M_0 + M_d)*p * (M_0 + M_d)*p);
        // esRate
        Rates[1] = EB * (1.0 - stateIN.B) * (1.0 - stateIN.B) * std::sin(w) * std::sin(w) / EBc
                   * q / (p*p + q*q)
                   + q / ((M_0 + M_d) * (p + c*(1 - stateIN.D)) * (M_0 + M_d) * (p + c*(1 - stateIN.D)));

        // DRate
        Rates[2] = (1.0 - stateIN.D) * (1.0 - stateIN.D) / EDc;

        // BRate
        Rates[3] = (1.0 - stateIN.B) * (1.0 - stateIN.B) * std::cos(w) * std::cos(w) / EBc;

    // No Yielding
    } else {

        // Zero  Yielding
        Rates[0] = 0; //evRate
        Rates[1] = 0; //esRate
        Rates[2] = 0; //DRate
        Rates[3] = 0; //BRate
    }

    if (ED < 0){
        Rates[2] = 0; //DRate
    }

    // Correct for Degenerate Case
    if (EB < 0){
        Rates[3] = 0; //BRate
    }

    return Rates;

}


std::vector<double> CompressibleDamageMechanicsSandstone::YieldFunctionsFromMaterialState(MaterialState& stateIN){
//Compute y from material state

    // Compute Material Stress State
    MaterialTensor Sy = YieldStressFromMaterialState(stateIN);

    // Compute Pressure and Shear Stress in Deformed Reference Frame
    double p = -Sy.trace()/3.0;
    MaterialTensor s0 = Sy + p * MaterialTensor::Identity();
    double q = std::sqrt(3.0 / 2.0 * s0.dot(s0));

    // Compute ED,EB
    std::vector<double> EDEB = EDEBFromMaterialState(stateIN);
    double ED = EDEB[0];
    double EB = EDEB[1];
    if (ED < 0){
        ED = 0;
    }

    // Compute Limiting Porosities
    double phi_min = phi_l * std::pow(1.0 - stateIN.B, l);
    double phi_max = phi_u * std::pow(1.0 - stateIN.B, u);

    // Compute Relative Porosity
    double tau = (phi_max - stateIN.phiP) / (phi_max - phi_min);

    // Compute Critical Porosity
    double tau_cs = std::sqrt(EB/EBc) * (1.0 - stateIN.B) / g;

    // Compute Maximum Friction Angle
    double tp = M_PI / 15.0 + std::asin(3.0 * M_0 / (6.0 + M_0));

    // Dilation Angle
    double M_d = g * (tau - tau_cs) * (6.0 * std::sin(tp) / (3.0 - std::sin(tp)) - M_0);

    // Handle B = 1
    if (!std::isfinite(tau)){
        M_d = 0;
    }

    // Identify Relevant Yield Functions
    double y1, y2, y3;

    // Check for Disperse Yielding
    if (stateIN.D >= 1 && (1 - stateIN.rho/rho_0) > phi_max){
        y1 = -1; // NOT Breakage/Damage
        y2 = -1; // NOT Tension
        y3 = 0; // YES Free Compression

    // Check for Tensile Yielding
    } else if (p < -c * (1 - stateIN.D)){
        y1 = -1; // NOT Breakage/Damage
        y2 = -c * (1.0 - stateIN.D) - p;  // YES Tension
        y3 = -1; // NOT Free Compression

    // Check for Compressive Yielding
    } else {
        y1 = ED * (1.0 - stateIN.D) * (1.0 - stateIN.B) / EDc
                + EB * (1.0 - stateIN.B) * (1.0 - stateIN.B) / EBc
                + q * q / ((M_0 + M_d)*(p + c * (1.0 - stateIN.D)) * (M_0 + M_d)*(p + c * (1.0 - stateIN.D))) - 1.0;

        y2 = -c * (1.0 - stateIN.D) - p;

        if ((1.0 - stateIN.rho/rho_0) < phi_max){
            y3 = 1.0 - stateIN.rho/rho_0 - phi_max;
        } else {
            y3 = stateIN.D - 1.0;
        }
    }

    // Yield f'n Values
    std::vector<double> y = std::vector<double>(3);
    y[0] = y1;
    y[1] = y2;
    y[2] = y3;

    // Return Something
    return y;

}

MaterialTensor CompressibleDamageMechanicsSandstone::CauchyStressFromMaterialState(MaterialState& stateIN){
    // Compute Cauchy Stress from Material State and Deformation

    // Strain Invariants
    MaterialTensor Ee = 0.5 * (stateIN.Be - MaterialTensor::Identity());
    MaterialTensor Be0 = stateIN.Be.deviator();
    MaterialTensor Ee0 = Ee.deviator();
    double eve = Ee.trace();
    double ese = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));
    double Je  = std::sqrt(stateIN.Be.det());

    // Compute P,Q from Helmholtz free energy
    std::vector<double> pq = PQFromMaterialState(stateIN);
    double pBARD = pq[0];
    double qBARD = pq[1];
    double pBARG = pq[2];
    double qBARG = pq[3];
    double pSTAR = pq[4];

    // Compute Combined P,Q
    double pBAR = stateIN.D * pBARG + (1.0 - stateIN.D) * pBARD;
    double qBAR = stateIN.D * qBARG + (1.0 - stateIN.D) * qBARD;

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
    if (is_compressible && pSTAR > 0){
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

MaterialTensor CompressibleDamageMechanicsSandstone::YieldStressFromMaterialState(MaterialState& stateIN){
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
    double pBARD = pq[0];
    double qBARD = pq[1];
    double pBARG = pq[2];
    double qBARG = pq[3];
    double pSTAR = pq[4];

    // Compute Combined P,Q
    double pBAR = stateIN.D * pBARG + (1.0 - stateIN.D) * pBARD;
    double qBAR = stateIN.D * qBARG + (1.0 - stateIN.D) * qBARD;

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
    if (is_compressible && pSTAR > 0){
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

double CompressibleDamageMechanicsSandstone::PorosityFromMaterialState(MaterialState &stateIN) {
    // Estimate Porosity from Material Deformation

    // Elastic Volume Ratio
    double Je  = std::sqrt(stateIN.Be.det());

    // Compute Porosity Using Bisection
    double J = 1.0;
    double a, phiSTemp;
    double phiA, phiB, phiC;
    double yA, yB, yC, dydp;

    // Bisection Parameters
    int n = 0;
    double A = (1.0/Je - 1.0);
    phiA = 0; yA = -1;
    phiB = 1; yB =  1;
    phiC = 1; yC =  1; dydp = 1;

    // Check for Monotonicity
    bool is_monotonic = (Je <= 1.0);

    // If Je > 1, Then Choose New phiB, phiC
    if (Je > 1.0){
        phiB = std::pow(-1.0 / (A * (b + 1.0)), 1.0/b);
        phiC = phiB;
    }

    // Bisection or Newton's Method
    if (is_compressible && std::isfinite(b)){
        // Use Compressible Model

        while (std::abs(yC) > AbsTOL && n < MaxIter){
            //increment n
            n++;

            if (use_newtons_method && is_monotonic){
                // Use Newton's Method

                //evaluate yC
                yC = phiC + A * std::pow(phiC, b + 1) - stateIN.rho / rho_0;

                //evaluate dy/dphi
                dydp = (b + 1) * A * std::pow(phiC, b) + 1.0;

                //update phiC
                phiC -= yC / dydp;

                //ensure phiC physical
                if (phiC > 1.0){
                    phiC = 1.0;
                } else if (phiC < 0.0){
                    phiC = 0.0;
                }

            } else {
                // Use Bisection

                //set phiC to average of phiA and phiB
                phiC = 0.5 * (phiA + phiB);

                //evaluate yC
                yC = phiC + A * std::pow(phiC, b + 1) - stateIN.rho / rho_0;

                // assign phiA, phiB
                if (yC * yA > 0) {
                    phiA = phiC;
                    yA = yC;
                } else {
                    phiB = phiC;
                    yB = yC;
                }
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


double CompressibleDamageMechanicsSandstone::StrainEnergyFromMaterialState(MaterialState &stateIN) {
    // Compute Contact Strain Energy (J/kg) from Material Deformation

    // Strain Invaraints
    MaterialTensor Ee = 0.5 * (stateIN.Be - MaterialTensor::Identity());
    MaterialTensor Be0 = stateIN.Be.deviator();
    MaterialTensor Ee0 = Ee.deviator();
    double eve = Ee.trace();
    double ese = std::sqrt(2.0 / 3.0 * Ee0.dot(Ee0));

    // Compute Specific Porous Rock Strain Energy
    double psiD = (1.0 / rho_0) * (Kd * eve * eve / 2.0 + 3.0 * Gd * ese * ese / 2.0);

    // Compute Specific Granular Strain Energy
    if (eve > 0){
        eve = 0;
    }
    double psiG = (1.0/ rho_0) * (1.0 - theta * stateIN.B) * pr * (-Kg * Kg * eve * eve * eve / 12.0
                                                                   - 3.0 * Gg * Kg * eve * ese * ese / 4.0
                                                                   + Gg * std::sqrt(3.0 * Kg * Gg) * ese * ese * ese / 4.0);

    return stateIN.D * psiG + (1.0 - stateIN.D) * psiD;
}


double CompressibleDamageMechanicsSandstone::ColdEnergyFromMaterialState(MaterialState &stateIN) {
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
    if (is_compressible && std::isfinite(b)){
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


double CompressibleDamageMechanicsSandstone::TemperatureFromMaterialState(MaterialState &stateIN) {
    // Compute Temperature from Material State

    // Compute Small Deformation Strain Energy
    double psiE = StrainEnergyFromMaterialState(stateIN);

    // Compute Cold Strain Energy
    double eC = ColdEnergyFromMaterialState(stateIN);

    // Compute Temperature or Energy
    if (is_adiabatic){
        // Use Thermal Version of Model
        stateIN.Ts = (stateIN.Es - psiE - eC)/cv;
        return stateIN.Ts;

    } else {
        // Use Isothermal Version of Model
        stateIN.Es = psiE + eC + cv * T0;
        stateIN.Ts = T0;
        return T0;
    }

}
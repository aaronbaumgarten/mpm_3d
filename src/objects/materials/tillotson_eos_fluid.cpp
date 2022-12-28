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
    if (fp64_props.size() < 15){
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
        std::cerr << "    eta0   [Pa*s]   Viscosity\n";
        std::cerr << "    PCav   [Pa]     Cavitation Pressure\n";

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
        eta0    = fp64_props[13];
        PCav    = fp64_props[14];

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
        std::cout << "    eta0   = " << eta0 << " [Pa*s]\n";
        std::cout << "    PCav   = " << PCav << " [Pa]\n";

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

        //tillotson EOS has E = Ec + cv * (T - T0)
        Ef = Eigen::VectorXd::Zero(body->points->x.size()); //cv * T0 * Eigen::VectorXd::Ones(body->points->x.size());
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

    //material state for internal variables
    MaterialState mat_state;

    //tensors for matrix calculation
    MaterialTensor S, Sn, L, D, W;

    //artificial viscosity
    double C0 = std::sqrt(A/r0);

    for (int i=0;i<body->points->x.size();i++){
        if (body->points->active[i] == 0){
            continue;
        }

        // Initialize Material State Structure
        mat_state = MaterialState();
        mat_state.rho  = body->points->m(i) / body->points->v(i);
        mat_state.Tf   = Tf(i);
        mat_state.Ef   = Ef(i);
        mat_state.S    = body->points->T[i];
        mat_state.p    = -body->points->T[i].trace()/3.0;

        // Strain-Rate
        L = body->points->L[i];

        // Stress Tensor
        S = body->points->T[i];
        Sn = S;

        // Deformation and Spin Rate
        D = 0.5*(L+L.transpose());
        W = 0.5*(L-L.transpose());

        // Estimate Pressure
        mat_state.p = PressureFromMaterialState(mat_state);

        // Estimate Material Stress
        mat_state.S = CauchyStressFromMaterialState(mat_state);

        // Check for Cavitation
        if (mat_state.p >= PCav) {

            // Add Artificial Viscosity in Compression
            if (use_artificial_viscosity && D.trace() < 0) {
                mat_state.S += r0 * C0 * h_i * D.trace() * MaterialTensor::Identity();
            }

            // Compute Temperature/Energy Increment
            mat_state.Ef += 0.5 * (Sn.dot(D) + mat_state.S.dot(D)) * job->dt / mat_state.rho;
            mat_state.Tf = TemperatureFromMaterialState(mat_state);

            // Final Stress Update
            body->points->T[i] = CauchyStressFromMaterialState(mat_state);

            // Final State Updtae
            Tf(i) = mat_state.Tf;
            Ef(i) = mat_state.Ef;
            rho(i) = mat_state.rho;
            J(i) = r0 / mat_state.rho;

            // Add Artificial Viscosity in Compression
            if (use_artificial_viscosity && D.trace() < 0) {
                body->points->T[i] += r0 * C0 * h_i * D.trace() * MaterialTensor::Identity();
            }

        } else {

            // Fluid has Cavitated

            // Set Stress to Cavitation Pressure
            body->points->T[i] = -PCav * MaterialTensor::Identity();

            // Set Temperature to T0
            Tf(i) = mat_state.Tf;
            Ef(i) = 0.0;
            rho(i) = mat_state.rho;
            J(i) = r0 / mat_state.rho;
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

    // Ensure No Cavitation
    if (pressureIN < PCav){
        // Do Nothing
        return;
    }

    // Assign Stress
    MaterialTensor tmp = body->points->T[idIN];
    body->points->T[idIN] = tmp - (1.0/3.0 * tmp.trace() + pressureIN)*MaterialTensor::Identity();

    // Approximate Density Calculation
    double mu = (-A + std::sqrt(A*A + 4.0*B*pressureIN)) / (2.0 * B);
    double eta = mu + 1;
    double r = r0 * eta;

    // Assign Density, Volume, and Mass of MPM Points
    body->points->v0(idIN) *= eta;
    body->points->m(idIN) *= eta;

    // Create MaterialState Structure
    MaterialState mat_state = MaterialState();
    mat_state.rho = r;
    mat_state.Tf = Tf(idIN);
    mat_state.Ef = ColdEnergyFromMaterialState(mat_state);

    // Compute Stress
    body->points->T[idIN] = CauchyStressFromMaterialState(mat_state);

    return;
}


/*----------------------------------------------------------------------------*/
//frame and state writing
void TillotsonEOSFluid::writeFrame(Job* job, Body* body, Serializer* serializer){

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
    serializer->writeScalarArray(Ef, "Ef");
    serializer->writeScalarArray(Tf, "Tf");
    serializer->writeScalarArray(rho, "rho");
    serializer->writeScalarArray(J, "J");

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
    rVec = std::vector<double>(nBins + 1);
    eCVec = std::vector<double>(nBins + 1);

    // Initialize Jmin, Jmax, and dJ
    rmax = 5.0 * r0;                     // 5x Density
    rmin = 0.5 * r0;                     // 0.5x Density
    dr   = (rmax - rmin) / nBins;        // Volume Ratio Increment

    // Correct rmin, rmax s.t. r = r0 gives eCVec = 0
    int zBin = (int)((r0 - rmin)/dr);
    rmin = r0 - dr * zBin;             // Corrected Minimum
    rmax = r0 + dr * (nBins - zBin);   // Corrected Maximum
    rVec[zBin] = r0;
    eCVec[zBin] = 0.0;

    // Intermediate Variables
    double rp, P;
    MaterialState stateCOLD = MaterialState();

    // Compute Cold Curve for rho > r0 (Compression)
    for (int n = zBin+1; n < nBins+1; n++){

        rVec[n] = rVec[n-1] + dr;               // Reference Density

        rp = 0.5 * (rVec[n] + rVec[n-1]);       // Density for Integration

        stateCOLD.rho = rp;                     // Cold Density State
        stateCOLD.Ef  = eCVec[n-1];             // Cold Energy State

        P = PressureFromMaterialState(stateCOLD); // Pressure Calculation

        eCVec[n] = eCVec[n-1] + dr * P / (rp*rp);   // Integrate Cold Energy
    }

    // Compute Cold Curve for rho < r0 (Expansion)
    for (int n = zBin-1; n >= 0; n--){

        rVec[n] = rVec[n+1] - dr;               // Reference Density

        rp = 0.5 * (rVec[n+1] + rVec[n]);       // Density for Integration

        stateCOLD.rho = rp;                     // Cold Density State
        stateCOLD.Ef  = eCVec[n+1];             // Cold Energy State

        P = PressureFromMaterialState(stateCOLD); // Pressure Calculation

        eCVec[n] = eCVec[n+1] - dr * P / (rp*rp);   // Integrate Cold Energy
    }

    // List Should Be Complete
    std::cout << "Cold Energy Curve Computed!" << std::endl;

    return;

}

double TillotsonEOSFluid::PressureFromMaterialState(MaterialState &stateIN) {

    // State Quantities of Interest
    double P, eta, mu, J;
    double P2, P3;

    // Density Ratio and Strain
    eta = stateIN.rho / r0;
    mu = eta - 1.0;
    J = r0 / stateIN.rho;

    // Initialize P
    P = 0;

    // Determine Pressure Regime
    if (stateIN.rho >= r0) {

        // Pressure Equation (1) from Brundage (2013)
        P = (a + b / (stateIN.Ef / (E0 * eta * eta) + 1.0)) *
                stateIN.rho * stateIN.Ef + A * mu + B * mu * mu;

    } else if (stateIN.rho < r0 &&
               stateIN.rho >= rIV &&
               stateIN.Ef <= EIV) {

        // Pressure Equation (2) from Brundage (2013)
        P = (a + b / (stateIN.Ef / (E0 * eta * eta) + 1.0)) *
            stateIN.rho * stateIN.Ef + A * mu + B * mu * mu;

    } else if (stateIN.rho <= r0 &&
               stateIN.Ef >= ECV) {

        // Pressure Equation (3) from Brundage (2013)
        P = a * stateIN.rho * stateIN.Ef
                + (b * stateIN.rho * stateIN.Ef / (stateIN.Ef / (E0 * eta * eta) + 1.0)
                + A * mu * std::exp(-beta * (J - 1.0)))
                * std::exp(-alfa * (J - 1.0) * (J - 1.0));

    } else if (stateIN.rho < r0 &&
               stateIN.rho >= rIV &&
               stateIN.Ef > EIV &&
               stateIN.Ef < ECV){

        // Pressure Equation (5) from Brundage (2013)

        P2 = (a + b / (stateIN.Ef / (E0 * eta * eta) + 1.0)) *
             stateIN.rho * stateIN.Ef + A * mu + B * mu * mu;

        P3 = a * stateIN.rho * stateIN.Ef
             + (b * stateIN.rho * stateIN.Ef / (stateIN.Ef / (E0 * eta * eta) + 1.0)
                + A * mu * std::exp(-beta * (J - 1.0)))
               * std::exp(-alfa * (J - 1.0) * (J - 1.0));

        P = ((stateIN.Ef - EIV)*P3 + (ECV - stateIN.Ef)*P2) / (ECV - EIV);

    } else if (stateIN.rho < rIV &&
               stateIN.Ef < ECV){

        // Pressure Equation (6) from Brundage (2013)
        P = (a + b / (stateIN.Ef / (E0 * eta * eta) + 1.0)) * stateIN.rho * stateIN.Ef
                + A * mu;

    } else {

        // Hmm... That isn't supposed to happen...
        P = 0.0;

    }

    // Check for Cavitation
    /*if (P < PCav){
        P = PCav;
    }*/

    return P;
}

MaterialTensor TillotsonEOSFluid::CauchyStressFromMaterialState(MaterialState& stateIN){
    // Compute Cauchy Stress from Material State and Deformation

    // Compute Pressure
    double P = PressureFromMaterialState(stateIN);

    // Compute S
    MaterialTensor S = MaterialTensor();
    S = eta0 * (stateIN.L + stateIN.L.transpose()) - P * MaterialTensor::Identity();

    return S;
}



double TillotsonEOSFluid::ColdEnergyFromMaterialState(MaterialState &stateIN) {
    // Compute Cold Energy (J/kg) from Material State
    double eC = 0;

    // Density from Material State
    double rho = stateIN.rho;

    // Cold Energy Computation
    double drp, rp, P;
    double rminus, rplus, eCminus, eCplus;
    int N = 25;
    int rBin;

    // Look Up or Compute eC
    if (rho > rmin && rho < rmax){
        // Use Look Up Table

        // Identify Bin for r
        rBin = (int) ((rho - rmin) / dr);

        // Identify Bin Variables
        rminus  = rVec[rBin];
        rplus   = rVec[rBin+1];
        eCminus = eCVec[rBin];
        eCplus  = eCVec[rBin+1];

        // Linear Interpolation w/in Bin
        eC = eCminus + (rho - rminus) * (eCplus - eCminus) / (rplus - rminus);

    } else {

        MaterialState stateTMP = MaterialState();

        // Compute Numerically
        for (int i=0; i<N; i++){
            drp = (rho - r0)/N;
            rp = r0 + (i + 0.5) * dr;

            stateTMP.rho = rp;
            stateTMP.Ef  = eC;

            P = PressureFromMaterialState(stateTMP);

            eC += dr * P / (rp * rp);
        }
    }

    return eC;
}


double TillotsonEOSFluid::TemperatureFromMaterialState(MaterialState &stateIN) {
    // Compute Temperature from Material State

    // Compute Cold Strain Energy
    double eC = ColdEnergyFromMaterialState(stateIN);

    // Compute Temperature or Energy
    if (is_adiabatic){
        // Use Thermal Version of Model
        stateIN.Tf = T0 + (stateIN.Ef - eC)/cv;
        return stateIN.Tf;

    } else {
        // Use Isothermal Version of Model
        stateIN.Ef = eC;
        stateIN.Tf = T0;
        return T0;
    }

}
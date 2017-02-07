//
// Created by aaron on 1/31/17.
// viscoelastic_polymer_liquid.cpp
//

/*
	standard include for files material source files
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include "particle.hpp"
#include "node.hpp"
#include "body.hpp"
#include "process.hpp"
#include "tensor.hpp"
#include <Eigen/Dense>

/*
	definitions
 */
#define xTOL 0.6597 //slope of inverse langevin reaches 10
#define INVERSE_LANGEVIN(x) (x*(3-x*x)/(1-x*x))// < 5 ? (x*(3-x*x)/(1-x*x)) : 5//(x/abs(x))*(xTOL*(3-xTOL*xTOL)/(1-xTOL*xTOL))+10*(x-xTOL)//A. Cohen approx.
#define THRESHOLD 0.0001 //weak threshold for error in newton approx.
#define TAUBAR_THRESHOLD 1e-11 // threshold for negligible taubar in calculation
#define FP_THRESHOLD 1e-11 //threshold for Fp newton-rhapson method

/*
	undefine by Sachith to free up potential variable names used elsewhere
*/
#undef EMOD
#undef NUMOD
#undef G
#undef K

/*
    material properties for model defined in configuration file
*/
double mus; /* solvent viscosity 3000 */ 
double GR; /* rubbery modulus 10000 */
double lambda_L; /* network locking strength 10 */
double K; /* bulk modulus for weak incomp. */
double cstC; /* 0 */
double eta0; /* 60000 */
double zeta; /* 4.7 */
double m; /* 0.65 */

/*
	declare functions
*/
extern "C" void material_init(Body *body);

extern "C" void calculate_stress_threaded(threadtask_t *task, Body *body, double dt);

extern "C" void calculate_stress(Body *body, double dt);

extern "C" void calculate_stress_implicit(Body *body, double dt);

void material_init(Body *body)
{
    int i, j;

    // store F (F(t=0) = 1)
    for (i = 0; i < body->p; i++) {
        for (j = 0; j < NDIM*NDIM; j++) {
            body->particles.F(i,j) = 0;
        }
        body->particles.F(i,XX) = 1;
        body->particles.F(i,YY) = 1;
        body->particles.F(i,ZZ) = 1;
        body->particles.Fp(i,XX) = 1;
        body->particles.Fp(i,YY) = 1;
        body->particles.Fp(i,ZZ) = 1;
    }

    if (body->material.num_fp64_props < 3) {
        // Bit of a hack, but it's okay for now. Just close your eyes and code it anyways.
        std::cout << body->material.num_fp64_props << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 3 properties defined (mus, GR, lambda_L, [K], [C, eta0, zeta, m]).\n",
                __FILE__, __func__);
        exit(EXIT_FAILURE);
    }

    // set material properties from config file
    if (body->material.num_fp64_props < 8) {
        cstC = 0;
        eta0 = 60000.0;
        zeta = 4.7;
        m = 0.65;
    } else {
        cstC = body->material.fp64_props[4];
        eta0 = body->material.fp64_props[5];
        zeta = body->material.fp64_props[6];
        m = body->material.fp64_props[7];
    }

    if (body->material.num_fp64_props < 4){
        K = 1e6;
    } else {
        K = body->material.fp64_props[3];
    }

    mus = body->material.fp64_props[0];
    GR = body->material.fp64_props[1];
    lambda_L = body->material.fp64_props[2];

    printf("Material properties (mus = %g Pas, GR = %g Pa, lambda_L = %g, K = %g, cstC = %g, eta0 = %g Pas, zeta = %g s, m = %g).\n",
           mus, GR, lambda_L, K, cstC, eta0, zeta, m);

    std::cout << "Done initializing material (" << body->id << ").\n";
    return;
}
/*----------------------------------------------------------------------------*/

/* Viscoelastic polymer fluid model. */
void calculate_stress(Body *body, double dt)
{
    threadtask_t t;
    calculate_stress_threaded(&t, body, dt);
    return;
}

/*----------------------------------------------------------------------------*/
void calculate_stress_threaded(threadtask_t *task, Body *body, double dt)
{
    int i, j;

    for (i = 0; i < body->p; i++) {
        if (body->particles.active[i] == 0) {
            continue;
        }

        // arrays for intermediate calculations
        Eigen::Matrix3d tmp;
        Eigen::Matrix3d tmp2;
        Eigen::Matrix3d tmp3;
        Eigen::VectorXd tmpVec(9);

        // deformation gradient
        tmpVec << body->particles.F.row(i).transpose();
        Eigen::Matrix<double,3,3,Eigen::RowMajor> F(tmpVec.data());

        // plastic distortion
        tmpVec << body->particles.Fp.row(i).transpose();
        Eigen::Matrix<double,3,3,Eigen::RowMajor> Fp(tmpVec.data());

        // implicit calculation of plastic deformation Fp
        Eigen::Matrix3d RFp = Eigen::Matrix3d::Identity();
        double magRFp = 0;
        Eigen::Matrix3d FpTAU = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d FpOLD = Fp; //store Fp at time t as FpOLD
        Eigen::Matrix3d FOLD = F; //store F at time t as FOLD

        //elastic distortion
        Eigen::Matrix3d Fe;
        double Je;
        int loopCount = 0;

        // loop declarations to avoid compile errors
        Eigen::Matrix3d Se, Sv, D, Be, Me, Dp, Fpdot;
        Eigen::Matrix<double,3,3,Eigen::RowMajor> L;
        double G, gammadotp, I1, I2, gammabar_dot, lambdabar_e, taubar;
        double gammap = body->particles.state(i,9);
        gammadotp = body->particles.state(i,10);
        double gammadotp2 = gammadotp;
        double R;
        double dR;
        double R2;
        double dR2;
        do {
            loopCount += 1;
            if (loopCount>100 && loopCount%100==0){
                printf("Number of cycles in material exceeded %d! R: %g > %g\n",loopCount,magRFp,FP_THRESHOLD);
                //loopCount = 0;
            }
            // elastic distortion, transpose, and inverse transpose
            //double Fe[9];
            //FeT.setZero();
            //FeTinv.setZero();

            //Fp^T * Fe^T = F^T
            Fe = (Fp.transpose().colPivHouseholderQr().solve(F.transpose())).transpose();
            Je = Fe.determinant();
            Fe *= 1.0/std::cbrt(Je);

            //inverses
            //tensor_transpose3(FeT,Fe);
            //tensor_inverse3(FeTinv,FeT);

            // velocity gradient
            tmpVec << body->particles.L.row(i).transpose();
            L = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor>>(tmpVec.data());

            // deformation strain rate
            D = 0.5*(L + L.transpose());

            // elastic left Cauchy-Green tensor, first and second invariant
            Be = Fe*Fe.transpose();
            I1 = Be.trace();
            I2 = 0.5*(I1*I1 - (Be*Be).trace());

            // scalar equivalent shearing rate (messy calculation)
            gammabar_dot = std::sqrt(2)*D.norm();

            // viscous stress
            Sv = 2*mus*(D - D.trace()/3.0 * Eigen::Matrix3d::Identity());

            // effective elastic stretch
            lambdabar_e = std::sqrt(I1/3.0);

            // modulus G
            //double G;
            G = GR*lambda_L/(3*lambdabar_e)*INVERSE_LANGEVIN(lambdabar_e/lambda_L);
            //printf("G = %g\tLe/Ll = %g\n",G,lambdabar_e/lambda_L);
            if (lambdabar_e > lambda_L){
                printf("Inverse Langevin function approximation not valid for values greater than 1!\n");
            }

            /*if (lambdabar_e/lambda_L > xTOL){
                printf("Approaching Locking @ t:%g\n",job->t);
            }*/

            // elastic stress (really messy calculation)
            //double Se[9];
            tmp = I1*Be - Be*Be;
            Se = G*(Be - Be.trace()/3.0 * Eigen::Matrix3d::Identity()) + cstC/I2 * (tmp - tmp.trace()/3.0 * Eigen::Matrix3d::Identity());

            // symetric Mandel stress
            tmp = Fe.colPivHouseholderQr().solve(Se.transpose()).transpose();
            Me = Fe.transpose()*tmp;

            // equivalent shear stress for plastic flow
            taubar = 1.0/std::sqrt(2) * (Me - Me.trace()/3.0 * Eigen::Matrix3d::Identity()).norm();

            // plastic shear rate by solving strength relation (bi-directional newton-rhapson method)
            if (std::abs(taubar) > TAUBAR_THRESHOLD){ //if taubar is zero, then gammadotp must be zero
                do{
                    R = std::pow(zeta*gammadotp,m) - eta0*gammadotp/taubar + 1; //solve remainder of gammadotp eq.
                    dR = m*std::pow(zeta,m)*pow(gammadotp,m-1) - eta0/taubar; //solve derivitive
                    R2 = std::pow(zeta*gammadotp2,m) - eta0*gammadotp/taubar + 1;
                    dR2 = m*std::pow(zeta,m)*pow(gammadotp2,m-1) - eta0/taubar;
                    //printf("N-R %g,%g,%g,%g\n",R,dR,gammadotp,taubar);
                    if (dR >= 0) {
                        dR = -1; //set dR to negative to avoid getting caught near 0
                    }
                    if (dR2 <= 0) {
                        dR2 = 1; //set dR to positive to avoid getting caught near 0
                    }
                    gammadotp -= R/dR;
                    gammadotp2 -= R2/dR2;
                    if (gammadotp < 0 || !std::isfinite(gammadotp)){
                        gammadotp = 0;
                    }
                    if (gammadotp2 < 0 || !std::isfinite(gammadotp2)){
                        gammadotp2 = 0;
                    }
                } while (R*R > THRESHOLD && R2*R2 > THRESHOLD);
            } else {
                gammadotp = 0;
            }

            if (gammadotp != 0 && R*R > R2*R2){
                gammadotp = gammadotp2; //if reverse N-R is better than forward N-R, replace gammadotp
            }

            //printf("done\n");

            // plastic flow rate
            if (std::abs(taubar) > TAUBAR_THRESHOLD) {
                Dp = gammadotp * (Me - Me.trace() / 3.0 * Eigen::Matrix3d::Identity()) / (2 * taubar);
            } else {
                Dp.setZero();
            }

            /*--------------------------------------------------------------------*/
            /*																	  */
            /*																	  */
            /*					BEGIN CALCULATIONS FOR t + dt					  */
            /*																	  */
            /*																	  */
            /*--------------------------------------------------------------------*/
            // BEGIN CALCULATIONS FOR t + dt

            // plastic deformation at t+tau
            FpTAU = FpOLD * (Eigen::Matrix3d::Identity() + dt*Dp);

            // deformation F
            F = FOLD * (Eigen::Matrix3d::Identity() + dt*L);

            // compare FpTAU and Fp then store FpTAU as Fp
            RFp = Fp - FpTAU;
            magRFp = RFp.norm();
            Fp = FpTAU;
        } while (magRFp > FP_THRESHOLD);

        /*if (!std::isfinite(Fp(0))){
            std::cout << "\n" << Fp << "\n" << FpOLD << "\n" << Dp << "\n----------\n";
        }*/

        //printf("done\n\n");

        // adding weak incompressibility****************************************************************
        Se += K*std::log(std::cbrt(Je*Je))*Eigen::Matrix3d::Identity();
        //**********************************************************************************************

        // Cauchy (plane) stress calculation and saving of values
        for (size_t pos=0;pos<NDIM*NDIM;pos++){
            body->particles.T(i,pos) = Sv(pos) + Se(pos); //OK because T is symmetric
            body->particles.F(i,pos) = F(pos); //OK because F is RowMajor
            body->particles.Fp(i,pos) = Fp(pos); //OK because F is RowMajor
        }
        body->particles.state(i,9) += gammadotp*dt;
        body->particles.state(i,10) = gammadotp;
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
void calculate_stress_implicit(Body *body, double dt)
{
    //same as above, just don't alter state variables
    int i, j;

    for (i = 0; i < body->p; i++) {
        if (body->particles.active[i] == 0) {
            continue;
        }

        // arrays for intermediate calculations
        Eigen::Matrix3d tmp;
        Eigen::Matrix3d tmp2;
        Eigen::Matrix3d tmp3;
        Eigen::VectorXd tmpVec(9);

        // deformation gradient
        tmpVec << body->particles.F.row(i).transpose();
        Eigen::Matrix<double,3,3,Eigen::RowMajor> F(tmpVec.data());

        // plastic distortion
        tmpVec << body->particles.Fp.row(i).transpose();
        Eigen::Matrix<double,3,3,Eigen::RowMajor> Fp(tmpVec.data());

        // implicit calculation of plastic deformation Fp
        Eigen::Matrix3d RFp = Eigen::Matrix3d::Identity();
        double magRFp = 0;
        Eigen::Matrix3d FpTAU = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d FpOLD = Fp; //store Fp at time t as FpOLD
        Eigen::Matrix3d FOLD = F; //store F at time t as FOLD

        //elastic distortion
        Eigen::Matrix3d Fe;
        double Je;
        int loopCount = 0;

        // loop declarations to avoid compile errors
        Eigen::Matrix3d Se, Sv, D, Be, Me, Dp, Fpdot;
        Eigen::Matrix<double,3,3,Eigen::RowMajor> L;
        double G, gammadotp, I1, I2, gammabar_dot, lambdabar_e, taubar;
        double gammap = body->particles.state(i,9);
        gammadotp = body->particles.state(i,10);
        double gammadotp2 = gammadotp;
        double R;
        double dR;
        double R2;
        double dR2;
        do {
            loopCount += 1;
            if (loopCount>100 && loopCount%100==0){
                printf("Number of cycles in material exceeded %d! R: %g > %g\n",loopCount,magRFp,FP_THRESHOLD);
                //loopCount = 0;
            }
            // elastic distortion, transpose, and inverse transpose
            //double Fe[9];
            //FeT.setZero();
            //FeTinv.setZero();

            //Fp^T * Fe^T = F^T
            Fe = (Fp.transpose().colPivHouseholderQr().solve(F.transpose())).transpose();
            Je = Fe.determinant();
            Fe *= 1.0/std::cbrt(Je);

            //inverses
            //tensor_transpose3(FeT,Fe);
            //tensor_inverse3(FeTinv,FeT);

            // velocity gradient
            tmpVec << body->particles.L.row(i).transpose();
            L = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor>>(tmpVec.data());

            // deformation strain rate
            D = 0.5*(L + L.transpose());

            // elastic left Cauchy-Green tensor, first and second invariant
            Be = Fe*Fe.transpose();
            I1 = Be.trace();
            I2 = 0.5*(I1*I1 - (Be*Be).trace());

            // scalar equivalent shearing rate (messy calculation)
            gammabar_dot = std::sqrt(2)*D.norm();

            // viscous stress
            Sv = 2*mus*(D - D.trace()/3.0 * Eigen::Matrix3d::Identity());

            // effective elastic stretch
            lambdabar_e = std::sqrt(I1/3.0);

            // modulus G
            //double G;
            G = GR*lambda_L/(3*lambdabar_e)*INVERSE_LANGEVIN(lambdabar_e/lambda_L);
            //printf("G = %g\tLe/Ll = %g\n",G,lambdabar_e/lambda_L);
            if (lambdabar_e > lambda_L){
                printf("Inverse Langevin function approximation not valid for values greater than 1!\n");
            }

            /*if (lambdabar_e/lambda_L > xTOL){
                printf("Approaching Locking @ t:%g\n",job->t);
            }*/

            // elastic stress (really messy calculation)
            //double Se[9];
            tmp = I1*Be - Be*Be;
            Se = G*(Be - Be.trace()/3.0 * Eigen::Matrix3d::Identity()) + cstC/I2 * (tmp - tmp.trace()/3.0 * Eigen::Matrix3d::Identity());

            // symetric Mandel stress
            tmp = Fe.colPivHouseholderQr().solve(Se.transpose()).transpose();
            Me = Fe.transpose()*tmp;

            // equivalent shear stress for plastic flow
            taubar = 1.0/std::sqrt(2) * (Me - Me.trace()/3.0 * Eigen::Matrix3d::Identity()).norm();

            // plastic shear rate by solving strength relation (bi-directional newton-rhapson method)
            if (std::abs(taubar) > TAUBAR_THRESHOLD){ //if taubar is zero, then gammadotp must be zero
                do{
                    R = std::pow(zeta*gammadotp,m) - eta0*gammadotp/taubar + 1; //solve remainder of gammadotp eq.
                    dR = m*std::pow(zeta,m)*pow(gammadotp,m-1) - eta0/taubar; //solve derivitive
                    R2 = std::pow(zeta*gammadotp2,m) - eta0*gammadotp/taubar + 1;
                    dR2 = m*std::pow(zeta,m)*pow(gammadotp2,m-1) - eta0/taubar;
                    //printf("N-R %g,%g,%g,%g\n",R,dR,gammadotp,taubar);
                    if (dR >= 0) {
                        dR = -1; //set dR to negative to avoid getting caught near 0
                    }
                    if (dR2 <= 0) {
                        dR2 = 1; //set dR to positive to avoid getting caught near 0
                    }
                    gammadotp -= R/dR;
                    gammadotp2 -= R2/dR2;
                    if (gammadotp < 0 || !std::isfinite(gammadotp)){
                        gammadotp = 0;
                    }
                    if (gammadotp2 < 0 || !std::isfinite(gammadotp2)){
                        gammadotp2 = 0;
                    }
                } while (R*R > THRESHOLD && R2*R2 > THRESHOLD);
            } else {
                gammadotp = 0;
            }

            if (gammadotp != 0 && R*R > R2*R2){
                gammadotp = gammadotp2; //if reverse N-R is better than forward N-R, replace gammadotp
            }

            //printf("done\n");

            // plastic flow rate
            if (std::abs(taubar) > TAUBAR_THRESHOLD) {
                Dp = gammadotp * (Me - Me.trace() / 3.0 * Eigen::Matrix3d::Identity()) / (2 * taubar);
            } else {
                Dp.setZero();
            }

            /*--------------------------------------------------------------------*/
            /*																	  */
            /*																	  */
            /*					BEGIN CALCULATIONS FOR t + dt					  */
            /*																	  */
            /*																	  */
            /*--------------------------------------------------------------------*/
            // BEGIN CALCULATIONS FOR t + dt

            // plastic deformation at t+tau
            FpTAU = FpOLD * (Eigen::Matrix3d::Identity() + dt*Dp);

            // deformation F
            F = FOLD * (Eigen::Matrix3d::Identity() + dt*L);

            // compare FpTAU and Fp then store FpTAU as Fp
            RFp = Fp - FpTAU;
            magRFp = RFp.norm();
            Fp = FpTAU;
        } while (magRFp > FP_THRESHOLD);


        /*if (!std::isfinite(Fp(0))){
            std::cout << "\n" << Fp << "\n" << FpOLD << "\n" << Dp << "\n----------\n";
        }*/

        //printf("done\n\n");

        // adding weak incompressibility****************************************************************
        Se += K*std::log(std::cbrt(Je*Je))*Eigen::Matrix3d::Identity();
        //**********************************************************************************************

        // Cauchy (plane) stress calculation and saving of values
        for (size_t pos=0;pos<NDIM*NDIM;pos++){
            body->particles.Ttrial(i,pos) = Sv(pos) + Se(pos); //OK because T is symmetric
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
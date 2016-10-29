//
// Created by aaron on 9/8/16.
// tensor.cpp
//

/**
    \file tensor.c
    \author Sachith Dunatunga
    \date 04.06.12

    mpm_2d -- An implementation of the Material Point Method in 2D.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "particle.hpp"

#include "tensor.hpp"
#define NDEBUG
#include <assert.h>

/*----------------------------------------------------------------------------*/
void zero_array(double *array, size_t numel)
{
#ifdef __STDC_IEC_559__
    /*
        We conform to IEEE754-1985, so we can safely use memset to set all
        zero bytes (0b -> 0.0).
    */
    memset(array, 0, numel * sizeof(double));
#else
    size_t i;
    for (i = 0; i < numel; i++) {
        array[i] = 0;
    }
#endif
    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** \brief Calculates the deviatoric and spherical part of a tensor.
    \param A_0 Where to store the deviator.
    \param A Tensor to take the deviator of.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_trace(double*  trA, const double*  A, size_t dim)
{
    size_t i;
    assert(dim == 2 || dim == 3);

/*    if (dim == 2) {*/
/*        *trA = A[XX] + A[YY];*/
/*    } else if (dim == 3) {*/
/*        *trA = A[XX] + A[YY] + A[ZZ];*/
/*    } else {*/
/*        *trA = 0;*/
/*    }*/
    *trA = 0;
    for (i = 0; i < dim; i++) {
        *trA += A[i * dim + i];
    }

    return;
}

void tensor_trace3(double*  trA, const double*  A)
{
    size_t i;

    *trA = 0;
    for (i = 0; i < 3; i++) {
        *trA += A[i * 3 + i];
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** \brief Adds two tensors together and stores the result in a third.
    \param C Where to store the result (tensor).
    \param A First tensor (tensor).
    \param B Second tensor (tensor).
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.

    This function adds the B tensor to the A tensor and stores the result in
    a tensor C. You may use the same tensor for both A and C to overwrite the
    value of A.
*/
void tensor_add(double *C, const double *A, const double *B, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(C != NULL);
    assert(A != NULL);
    assert(B != NULL);

    if (C != A) {
        tensor_copy(C, A, dim);
    }

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            C[i*dim + j] += B[i*dim + j];
        }
    }

    return;
}

void tensor_add3(double *C, const double *A, const double *B)
{
    size_t i, j;
    assert(C != NULL);
    assert(A != NULL);
    assert(B != NULL);

    if (C != A) {
        tensor_copy3(C, A);
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            C[i*3 + j] += B[i*3 + j];
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** \brief Copies a tensor.
    \param C Where to store the result (tensor).
    \param A Tensor to copy (tensor).
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_copy(double*  C, const double*  A, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(C != NULL);
    assert(A != NULL);

    if (C == A) {
        return;
    }

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            C[i*dim + j] = A[i*dim + j];
        }
    }

    return;
}

void tensor_copy3(double*  C, const double*  A)
{
    size_t i, j;
    assert(C != NULL);
    assert(A != NULL);

    if (C == A) {
        return;
    }

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            C[i*3 + j] = A[i*3 + j];
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** \brief Calculates the deviatoric and spherical part of a tensor.
    \param A_0 Where to store the deviator (tensor).
    \param cI Where to store the spherical part of A (scalar).
    \param A Tensor to take the deviator of.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_decompose(double*  A_0, double*  c, const double*  A, size_t dim)
{
    size_t i;
    double trA;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);
    assert(A_0 != NULL);
    assert(c != NULL);
    assert(A != A_0);


    tensor_trace(&trA, A, dim);
    *c = trA / dim;

    tensor_copy(A_0, A, dim);
    for (i = 0; i < dim; i++) {
        A_0[i*dim + i] = A_0[i*dim + i] - (*c);
    }

    return;
}

void tensor_decompose3(double*  A_0, double*  c, const double*  A)
{
    size_t i;
    double trA;
    assert(A != NULL);
    assert(A_0 != NULL);
    assert(c != NULL);
    assert(A != A_0);


    tensor_trace3(&trA, A);
    *c = trA / 3.0;

    tensor_copy3(A_0, A);
    for (i = 0; i < 3; i++) {
        A_0[i*3 + i] = A_0[i*3 + i] - (*c);
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Multiples two tensors and stores the result in a third.
    \param C where to store result tensor.
    \param A left tensor.
    \param B right tensor.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_multiply(double*  C, const double*  A, const double*  B, size_t dim)
{
    size_t i, j, k;
    assert(dim == 2 || dim == 3);
    assert(C != NULL);
    assert(A != NULL);
    assert(B != NULL);
    assert(A != C);

    if (dim == 3) {
        tensor_multiply3_helper(C, A, B);
    } else {
        for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
                C[i*dim + j] = 0;
                for (k = 0; k < dim; k++) {
                    C[i*dim + j] += A[i*dim + k]*B[k*dim + j];
                }
            }
        }
    }

    return;
}

void tensor_multiply3(double*  C, const double*  A, const double*  B)
{
    assert(C != NULL);
    assert(A != NULL);
    assert(B != NULL);
    assert(A != C);

    tensor_multiply3_helper(C, A, B);

    return;
}

void tensor_multiply3_helper(double*  C, const double*  A, const double*  B)
{
    size_t i, j;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            C[i*3+j] = A[3*i+0] * B[3*0+j]
                       + A[3*i+1] * B[3*1+j]
                       + A[3*i+2] * B[3*2+j];
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Scales tensor by a constant.
    \param A tensor to scale.
    \param c scaling constant.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_scale(double *A, double c, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            A[i*dim + j] *= c;
        }
    }

    return;
}

void tensor_scale3(double *A, double c)
{
    size_t i, j;
    assert(A != NULL);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            A[i*3 + j] *= c;
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Gets symmetric part of tensor.
    \param C where to store sym(A).
    \param A tensor to symmetrize.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_sym(double *  C, const double *  A, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);
    assert(C != NULL);
    assert(A != C);

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            C[i*dim + j] = 0.5 *(A[i*dim + j] + A[j*dim + i]);
        }
    }

    return;
}

void tensor_sym3(double *  C, const double *  A)
{
    size_t i, j;
    assert(A != NULL);
    assert(C != NULL);
    assert(A != C);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            C[i*3 + j] = 0.5 *(A[i*3 + j] + A[j*3 + i]);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Gets skew part of tensor.
    \param C where to store skw(A).
    \param A tensor to skew.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_skw(double *  C, const double *  A, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);
    assert(C != NULL);
    assert(A != C);

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            C[i*dim + j] = 0.5 *(A[i*dim + j] - A[j*dim + i]);
        }
    }

    return;
}

void tensor_skw3(double *  C, const double *  A)
{
    size_t i, j;
    assert(A != NULL);
    assert(C != NULL);
    assert(A != C);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            C[i*3 + j] = 0.5 *(A[i*3 + j] - A[j*3 + i]);
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Tensor contraction of A and B.
    \param C where to store skw(A).
    \param A tensor to skew.
    \param dim Is this a 2D or 3D tensor? Pass in either 2 or 3.
*/
void tensor_contraction(double *c, double *A, double *B, size_t dim)
{
    size_t i, j;
    assert(dim == 2 || dim == 3);
    assert(A != NULL);
    assert(B != NULL);
    assert(c != NULL);

    *c = 0;
    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            (*c) += A[i*dim + j] * B[i*dim + j];
        }
    }

    return;
}

void tensor_contraction3(double *c, const double *A, const double *B)
{
    size_t i, j;
    assert(A != NULL);
    assert(B != NULL);
    assert(c != NULL);

    *c = 0;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            (*c) += A[i*3 + j] * B[i*3 + j];
        }
    }

    return;
}
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief Inverse of tensor A.
    \param Ainv where to store inv(A).
    \param A tensor to invert.
*/
void tensor_inverse3(double *  Ainv, const double *  A)
{
    const size_t dim = 3;
    const double det = A[0*dim + 0] * (A[1*dim + 1] * A[2*dim + 2] - A[2*dim + 1] * A[1*dim + 2]) -
                       A[0*dim + 1] * (A[1*dim + 0] * A[2*dim + 2] - A[1*dim + 2] * A[2*dim + 0]) +
                       A[0*dim + 2] * (A[1*dim + 0] * A[2*dim + 1] - A[1*dim + 1] * A[2*dim + 0]);

    double invdet = 1 / det;

    Ainv[0*dim + 0] = (A[1*dim + 1] * A[2*dim + 2] - A[2*dim + 1] * A[1*dim + 2]) * invdet;
    Ainv[0*dim + 1] = (A[0*dim + 2] * A[2*dim + 1] - A[0*dim + 1] * A[2*dim + 2]) * invdet;
    Ainv[0*dim + 2] = (A[0*dim + 1] * A[1*dim + 2] - A[0*dim + 2] * A[1*dim + 1]) * invdet;
    Ainv[1*dim + 0] = (A[1*dim + 2] * A[2*dim + 0] - A[1*dim + 0] * A[2*dim + 2]) * invdet;
    Ainv[1*dim + 1] = (A[0*dim + 0] * A[2*dim + 2] - A[0*dim + 2] * A[2*dim + 0]) * invdet;
    Ainv[1*dim + 2] = (A[1*dim + 0] * A[0*dim + 2] - A[0*dim + 0] * A[1*dim + 2]) * invdet;
    Ainv[2*dim + 0] = (A[1*dim + 0] * A[2*dim + 1] - A[2*dim + 0] * A[1*dim + 1]) * invdet;
    Ainv[2*dim + 1] = (A[2*dim + 0] * A[0*dim + 1] - A[0*dim + 0] * A[2*dim + 1]) * invdet;
    Ainv[2*dim + 2] = (A[0*dim + 0] * A[1*dim + 1] - A[1*dim + 0] * A[0*dim + 1]) * invdet;

    return;
}
/*----------------------------------------------------------------------------*/


// Slightly modified version of  Stan Melax's code for 3x3 matrix diagonalization (Thanks Stan!)
// source: http://www.melax.com/diag.html?attredirects=0
// (answer from stackoverflow: http://stackoverflow.com/a/26949066)
void tensor_symeig3(double *  P, double *  D, const double *  A)
{
    // A must be a symmetric matrix.
    // returns P and D such that
    // Diagonal matrix D = PT * A * P;  and  A = P*D*PT
    const int maxsteps=24;  // certainly wont need that many.
    const size_t dim = 3;
    int k0, k1, k2;
    double o[3], m[3];
    double q [4] = {0.0,0.0,0.0,1.0};
    double jr[4];
    double sqw, sqx, sqy, sqz;
    double tmp1, tmp2, mq;
    double AP[9];
    double thet, sgn, t, c;
    for(int i=0;i < maxsteps;++i)
    {
        // quat to matrix
        sqx      = q[0]*q[0];
        sqy      = q[1]*q[1];
        sqz      = q[2]*q[2];
        sqw      = q[3]*q[3];
        P[0*dim + 0]  = ( sqx - sqy - sqz + sqw);
        P[1*dim + 1]  = (-sqx + sqy - sqz + sqw);
        P[2*dim + 2]  = (-sqx - sqy + sqz + sqw);
        tmp1     = q[0]*q[1];
        tmp2     = q[2]*q[3];
        P[1*dim + 0]  = 2.0 * (tmp1 + tmp2);
        P[0*dim + 1]  = 2.0 * (tmp1 - tmp2);
        tmp1     = q[0]*q[2];
        tmp2     = q[1]*q[3];
        P[2*dim + 0]  = 2.0 * (tmp1 - tmp2);
        P[0*dim + 2]  = 2.0 * (tmp1 + tmp2);
        tmp1     = q[1]*q[2];
        tmp2     = q[0]*q[3];
        P[2*dim + 1]  = 2.0 * (tmp1 + tmp2);
        P[1*dim + 2]  = 2.0 * (tmp1 - tmp2);

        // AP = A * P
        AP[0*dim + 0] = P[0*dim + 0]*A[0*dim + 0]+P[1*dim + 0]*A[0*dim + 1]+P[2*dim + 0]*A[0*dim + 2];
        AP[0*dim + 1] = P[0*dim + 1]*A[0*dim + 0]+P[1*dim + 1]*A[0*dim + 1]+P[2*dim + 1]*A[0*dim + 2];
        AP[0*dim + 2] = P[0*dim + 2]*A[0*dim + 0]+P[1*dim + 2]*A[0*dim + 1]+P[2*dim + 2]*A[0*dim + 2];
        AP[1*dim + 0] = P[0*dim + 0]*A[0*dim + 1]+P[1*dim + 0]*A[1*dim + 1]+P[2*dim + 0]*A[1*dim + 2];
        AP[1*dim + 1] = P[0*dim + 1]*A[0*dim + 1]+P[1*dim + 1]*A[1*dim + 1]+P[2*dim + 1]*A[1*dim + 2];
        AP[1*dim + 2] = P[0*dim + 2]*A[0*dim + 1]+P[1*dim + 2]*A[1*dim + 1]+P[2*dim + 2]*A[1*dim + 2];
        AP[2*dim + 0] = P[0*dim + 0]*A[0*dim + 2]+P[1*dim + 0]*A[1*dim + 2]+P[2*dim + 0]*A[2*dim + 2];
        AP[2*dim + 1] = P[0*dim + 1]*A[0*dim + 2]+P[1*dim + 1]*A[1*dim + 2]+P[2*dim + 1]*A[2*dim + 2];
        AP[2*dim + 2] = P[0*dim + 2]*A[0*dim + 2]+P[1*dim + 2]*A[1*dim + 2]+P[2*dim + 2]*A[2*dim + 2];
        // D = Pt * AP
        D[0*dim + 0] = AP[0*dim + 0]*P[0*dim + 0]+AP[1*dim + 0]*P[1*dim + 0]+AP[2*dim + 0]*P[2*dim + 0];
        D[0*dim + 1] = AP[0*dim + 0]*P[0*dim + 1]+AP[1*dim + 0]*P[1*dim + 1]+AP[2*dim + 0]*P[2*dim + 1];
        D[0*dim + 2] = AP[0*dim + 0]*P[0*dim + 2]+AP[1*dim + 0]*P[1*dim + 2]+AP[2*dim + 0]*P[2*dim + 2];
        D[1*dim + 0] = AP[0*dim + 1]*P[0*dim + 0]+AP[1*dim + 1]*P[1*dim + 0]+AP[2*dim + 1]*P[2*dim + 0];
        D[1*dim + 1] = AP[0*dim + 1]*P[0*dim + 1]+AP[1*dim + 1]*P[1*dim + 1]+AP[2*dim + 1]*P[2*dim + 1];
        D[1*dim + 2] = AP[0*dim + 1]*P[0*dim + 2]+AP[1*dim + 1]*P[1*dim + 2]+AP[2*dim + 1]*P[2*dim + 2];
        D[2*dim + 0] = AP[0*dim + 2]*P[0*dim + 0]+AP[1*dim + 2]*P[1*dim + 0]+AP[2*dim + 2]*P[2*dim + 0];
        D[2*dim + 1] = AP[0*dim + 2]*P[0*dim + 1]+AP[1*dim + 2]*P[1*dim + 1]+AP[2*dim + 2]*P[2*dim + 1];
        D[2*dim + 2] = AP[0*dim + 2]*P[0*dim + 2]+AP[1*dim + 2]*P[1*dim + 2]+AP[2*dim + 2]*P[2*dim + 2];
        o[0]    = D[1*dim + 2];
        o[1]    = D[0*dim + 2];
        o[2]    = D[0*dim + 1];
        m[0]    = fabs(o[0]);
        m[1]    = fabs(o[1]);
        m[2]    = fabs(o[2]);

        k0      = (m[0] > m[1] && m[0] > m[2])?0: (m[1] > m[2])? 1 : 2; // index of largest element of offdiag
        k1      = (k0+1)%3;
        k2      = (k0+2)%3;
        if (o[k0]==0.0)
        {
            break;  // diagonal already
        }
        thet    = (D[k2*dim + k2]-D[k1*dim + k1])/(2.0*o[k0]);
        sgn     = (thet > 0.0)?1.0:-1.0;
        thet   *= sgn; // make it positive
        t       = sgn /(thet +((thet < 1.E6)?sqrt(thet*thet+1.0):thet)) ; // sign(T)/(|T|+sqrt(T^2+1))
        c       = 1.0/sqrt(t*t+1.0); //  c= 1/(t^2+1) , t=s/c
        if(c==1.0)
        {
            break;  // no room for improvement - reached machine precision.
        }
        jr[0 ]  = jr[1] = jr[2] = jr[3] = 0.0;
        jr[k0]  = sgn*sqrt((1.0-c)/2.0);  // using 1/2 angle identity sin(a/2) = sqrt((1-cos(a))/2)
        jr[k0] *= -1.0; // since our quat-to-matrix convention was for v*M instead of M*v
        jr[3 ]  = sqrt(1.0f - jr[k0] * jr[k0]);
        if(jr[3]==1.0)
        {
            break; // reached limits of floating point precision
        }
        q[0]    = (q[3]*jr[0] + q[0]*jr[3] + q[1]*jr[2] - q[2]*jr[1]);
        q[1]    = (q[3]*jr[1] - q[0]*jr[2] + q[1]*jr[3] + q[2]*jr[0]);
        q[2]    = (q[3]*jr[2] + q[0]*jr[1] - q[1]*jr[0] + q[2]*jr[3]);
        q[3]    = (q[3]*jr[3] - q[0]*jr[0] - q[1]*jr[1] - q[2]*jr[2]);
        mq      = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
        q[0]   /= mq;
        q[1]   /= mq;
        q[2]   /= mq;
        q[3]   /= mq;
    }
    return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/** \brief determinant of tensor A.
    \param detA where to store det(A).
    \param A tensor to determine.
*/
void tensor_det3(double *  detA, const double *  A)
{
    *detA = A[0]*A[4]*A[8] + A[1]*A[5]*A[6] + A[2]*A[3]*A[7] - A[2]*A[4]*A[6] - A[1]*A[3]*A[8] - A[0]*A[5]*A[7];
    return;
}
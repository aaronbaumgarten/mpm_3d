//
// Created by aaron on 4/22/18.
// algebra.cpp
//

#include <math.h>
#include "mpmtensor.hpp"
#include <iostream>

/*------------------------------------------------------------------------*/
//construct MaterialTensor with pointer to data (not safe, but who cares right?)
MaterialTensor::MaterialTensor(double* otherdata){
    for (int i=0; i<TENSOR_MAX_LENGTH; i++){
        data_ptr[i] = otherdata[i];
    }
}


//construct MaterialTensor from a KinematicTensor
MaterialTensor::MaterialTensor(KinematicTensor& other){
    for(int i=0; i<TENSOR_MAX_LENGTH; i++){
        data_ptr[i] = other(i);
    }
}

/*
//construct MaterialTensor from another MaterialTensor
MaterialTensor::MaterialTensor(const MaterialTensor& other){
    for(int i=0; i<TENSOR_MAX_LENGTH; i++){
        data_ptr[i] = other(i);
    }
}
 */

/*------------------------------------------------------------------------*/
//define -A
MaterialTensor MaterialTensor::operator-(){
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp[i] = -data_ptr[i];
    }
    return MaterialTensor(tmp.data());
}

//define A *= s for integer values
MaterialTensor& MaterialTensor::operator*= (const int &rhs){
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = data_ptr[i]*rhs;
    }
    return *this;
}

//define A *= s for double values
MaterialTensor& MaterialTensor::operator*= (const double &rhs){
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = data_ptr[i]*rhs;
    }
    return *this;
}

//define A /= s for integer values
MaterialTensor& MaterialTensor::operator/= (const int &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = data_ptr[i]*tmp;
    }
    return *this;
}

//define A /= s for double values
MaterialTensor& MaterialTensor::operator/= (const double &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = data_ptr[i]*tmp;
    }
    return *this;
}

/*------------------------------------------------------------------------*/
//set A to zero tensor
void MaterialTensor::setZero(){
    for (int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = 0;
    }
    return;
}

//set A to identity tensor
void MaterialTensor::setIdentity(){
    data_ptr[XX] = 1; data_ptr[XY] = 0; data_ptr[XZ] = 0;
    data_ptr[YX] = 0; data_ptr[YY] = 1; data_ptr[YZ] = 0;
    data_ptr[ZX] = 0; data_ptr[ZY] = 0; data_ptr[ZZ] = 1;
    return;
}

/*------------------------------------------------------------------------*/
//tensor inner products
double MaterialTensor::dot(const MaterialTensor &rhs){
    double tmp=0;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp += data_ptr[i]*rhs[i];
    }
    return tmp;
}

double MaterialTensor::dot(const KinematicTensor &rhs){
    double tmp=0;
    for(int i=0;i<rhs.DIM;i++){
        for(int j=0;j<rhs.DIM;j++) {
            tmp += data_ptr[TENSOR_MAX_DIM*i+j] * rhs(i,j);
        }
    }
    return tmp;
}

//trace of MaterialTensor
double MaterialTensor::trace(){
    return data_ptr[XX]+data_ptr[YY]+data_ptr[ZZ];
}

//tensor norm given by sqrt(A:A)
double MaterialTensor::norm(){
    double tmp = 0;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp += data_ptr[i]*data_ptr[i];
    }
    return std::sqrt(tmp);
}

//determinant of MaterialTensor
double MaterialTensor::det(){
    double a = data_ptr[XX]*(data_ptr[YY]*data_ptr[ZZ] - data_ptr[ZY]*data_ptr[YZ]);
    double b = data_ptr[XY]*(data_ptr[YX]*data_ptr[ZZ] - data_ptr[ZX]*data_ptr[YZ]);
    double c = data_ptr[XZ]*(data_ptr[YX]*data_ptr[ZY] - data_ptr[ZX]*data_ptr[YY]);
    return a-b+c;
}

/*------------------------------------------------------------------------*/
//tensor inverse of MaterialTensor (not safe)
MaterialTensor MaterialTensor::inverse(){
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    double detA = det();
    tmp[XX] = 1/detA * (data_ptr[YY]*data_ptr[ZZ] - data_ptr[ZY]*data_ptr[YZ]);
    tmp[XY] = 1/detA * (data_ptr[XZ]*data_ptr[ZY] - data_ptr[ZZ]*data_ptr[XY]);
    tmp[XZ] = 1/detA * (data_ptr[XY]*data_ptr[YZ] - data_ptr[YY]*data_ptr[XZ]);
    tmp[YX] = 1/detA * (data_ptr[YZ]*data_ptr[ZX] - data_ptr[ZZ]*data_ptr[YX]);
    tmp[YY] = 1/detA * (data_ptr[XX]*data_ptr[ZZ] - data_ptr[ZX]*data_ptr[XZ]);
    tmp[YZ] = 1/detA * (data_ptr[XZ]*data_ptr[YX] - data_ptr[YZ]*data_ptr[XX]);
    tmp[ZX] = 1/detA * (data_ptr[YX]*data_ptr[ZY] - data_ptr[ZX]*data_ptr[YY]);
    tmp[ZY] = 1/detA * (data_ptr[XY]*data_ptr[ZX] - data_ptr[ZY]*data_ptr[ZX]);
    tmp[ZZ] = 1/detA * (data_ptr[XX]*data_ptr[YY] - data_ptr[YX]*data_ptr[XY]);
    for (int i=0;i<TENSOR_MAX_LENGTH;i++){
        if (!std::isfinite(tmp[i])){
            std::cerr << "WARNING: MaterialTensor::inverse() has resulted in non-finite value. So that's bad." << std::endl;
        }
    }
    return MaterialTensor(tmp.data());
}

//tensor deviator (A_0 = A - 1/3 tr(A) I)
MaterialTensor MaterialTensor::deviator(){
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    double trA = trace();
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp[i] = data_ptr[i];
    }
    tmp[XX] -= trA/3.0; tmp[YY] -= trA/3.0; tmp[ZZ] -= trA/3.0;
    return MaterialTensor(tmp.data());
}

//transpose (A^T)
MaterialTensor MaterialTensor::transpose(){
    std::array<double,9> tmp;
    tmp[XX] = data_ptr[XX]; tmp[XY] = data_ptr[YX]; tmp[XZ] = data_ptr[ZX];
    tmp[YX] = data_ptr[XY]; tmp[YY] = data_ptr[YY]; tmp[YZ] = data_ptr[ZY];
    tmp[ZX] = data_ptr[XZ]; tmp[ZY] = data_ptr[YZ]; tmp[ZZ] = data_ptr[ZZ];
    return MaterialTensor(tmp.data());
}

//symmetric and skew decompositions of A
MaterialTensor MaterialTensor::sym(){
    MaterialTensor tmp(this->data());
    return 0.5*(tmp + tmp.transpose());
}
MaterialTensor MaterialTensor::skw(){
    MaterialTensor tmp(this->data());
    return 0.5*(tmp - tmp.transpose());
}


/*------------------------------------------------------------------------*/
//construct KinematicTensor from data pointer and input tensor_type
KinematicTensor::KinematicTensor(double* otherdata, int input){
    assignTensorType(input);
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++){
            data_ptr[TENSOR_MAX_DIM * i + j] = otherdata[TENSOR_MAX_DIM * i + j];
        }
    }
    for (int i=DIM;i<3;i++) {
        for (int j = DIM; j < TENSOR_MAX_DIM; j++) {
            if (otherdata[3*i+j] != 0){
                std::cerr << "WARNING: ignoring non-zero entries in input data to KinematicTensor." << std::endl;
            }
            data_ptr[TENSOR_MAX_DIM * i + j] = 0;
        }
    }
}

/*
//copy KinematicTensor from KinematicTensor
KinematicTensor::KinematicTensor(const KinematicTensor& other){
    DIM = other.DIM;
    TENSOR_TYPE = other.TENSOR_TYPE;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = other[i];
    }
}

KinematicTensor::KinematicTensor(const KinematicTensor& other, int input){
    assignTensorType(input);
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++){
            data_ptr[TENSOR_MAX_DIM * i + j] = other[TENSOR_MAX_DIM * i + j];
        }
    }
    for (int i=DIM;i<3;i++) {
        for (int j = DIM; j < TENSOR_MAX_DIM; j++) {
            if (other[TENSOR_MAX_DIM * i + j] != 0){
                std::cerr << "WARNING: ignoring non-zero entries in input data to KinematicTensor." << std::endl;
            }
            data_ptr[TENSOR_MAX_DIM * i + j] = 0;
        }
    }
}
*/

//copy KinematicTensor from MaterialTensor
KinematicTensor::KinematicTensor(MaterialTensor& other, int input){
    assignTensorType(input);
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++){
            data_ptr[TENSOR_MAX_DIM * i + j] = other[TENSOR_MAX_DIM * i + j];
        }
    }
    for (int i=DIM;i<3;i++) {
        for (int j = DIM; j < TENSOR_MAX_DIM; j++) {
            if (other[TENSOR_MAX_DIM * i + j] != 0){
                std::cerr << "WARNING: ignoring non-zero entries in input data to KinematicTensor." << std::endl;
            }
            data_ptr[TENSOR_MAX_DIM * i + j] = 0;
        }
    }
}

/*------------------------------------------------------------------------*/
//define -A
KinematicTensor KinematicTensor::operator-(){
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp[i] = -data_ptr[i];
    }
    return KinematicTensor(tmp.data(),TENSOR_TYPE);
}

//define A *= s for integer values
KinematicTensor& KinematicTensor::operator*= (const int &rhs){
    if (DIM == 1){
        data_ptr[XX] *= rhs;
    } else if (DIM ==2){
        data_ptr[XX] *= rhs; data_ptr[XY] *= rhs;
        data_ptr[YX] *= rhs; data_ptr[YY] *= rhs;
    } else {
        for (int i = 0; i < TENSOR_MAX_LENGTH; i++) {
            data_ptr[i] = data_ptr[i] * rhs;
        }
    }
    return *this;
}

//define A *= s for double values
KinematicTensor& KinematicTensor::operator*= (const double &rhs){
    if (DIM == 1){
        data_ptr[XX] *= rhs;
    } else if (DIM ==2){
        data_ptr[XX] *= rhs; data_ptr[XY] *= rhs;
        data_ptr[YX] *= rhs; data_ptr[YY] *= rhs;
    } else {
        for (int i = 0; i < TENSOR_MAX_LENGTH; i++) {
            data_ptr[i] = data_ptr[i] * rhs;
        }
    }
    return *this;
}

//define A /= s for integer values
KinematicTensor& KinematicTensor::operator/= (const int &rhs){
    double tmp = 1.0/rhs;
    if (DIM == 1){
        data_ptr[XX] *= tmp;
    } else if (DIM ==2){
        data_ptr[XX] *= tmp; data_ptr[XY] *= tmp;
        data_ptr[YX] *= tmp; data_ptr[YY] *= tmp;
    } else {
        for (int i = 0; i < TENSOR_MAX_LENGTH; i++) {
            data_ptr[i] = data_ptr[i] * tmp;
        }
    }
    return *this;
}

//define A /= s for double values
KinematicTensor& KinematicTensor::operator/= (const double &rhs){
    double tmp = 1.0/rhs;
    if (DIM == 1){
        data_ptr[XX] *= tmp;
    } else if (DIM ==2){
        data_ptr[XX] *= tmp; data_ptr[XY] *= tmp;
        data_ptr[YX] *= tmp; data_ptr[YY] *= tmp;
    } else {
        for (int i = 0; i < TENSOR_MAX_LENGTH; i++) {
            data_ptr[i] = data_ptr[i] * tmp;
        }
    }
    return *this;
}

/*------------------------------------------------------------------------*/
//set to standard tensors
void KinematicTensor::setZero(){
    for (int i = 0; i < TENSOR_MAX_LENGTH; i++) {
        data_ptr[i] = 0;
    }
    return;
}

void KinematicTensor::setIdentity(){
    if (DIM == 1){
        data_ptr[XX] = 1;
    } else if (DIM ==2){
        data_ptr[XX] = 1; data_ptr[XY] = 0;
        data_ptr[YX] = 0; data_ptr[YY] = 1;
    } else {
        data_ptr[XX] = 1; data_ptr[XY] = 0; data_ptr[XZ] = 0;
        data_ptr[YX] = 0; data_ptr[YY] = 1; data_ptr[YZ] = 0;
        data_ptr[ZX] = 0; data_ptr[ZY] = 0; data_ptr[ZZ] = 1;
    }
    return;
}

/*------------------------------------------------------------------------*/
//tensor inner products
double KinematicTensor::dot(const MaterialTensor &rhs){
    double tmp=0;
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++) {
            tmp += data_ptr[TENSOR_MAX_DIM*i+j] * rhs(i,j);
        }
    }
    return tmp;
}

double KinematicTensor::dot(const KinematicTensor &rhs){
    double tmp=0;
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++) {
            tmp += data_ptr[TENSOR_MAX_DIM*i+j] * rhs(i,j);
        }
    }
    return tmp;
}

//trace of KinematicTensor
double KinematicTensor::trace(){
    return data_ptr[XX]+data_ptr[YY]+data_ptr[ZZ];
}

//tensor norm given by sqrt(A:A)
double KinematicTensor::norm(){
    double tmp = 0;
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++) {
            tmp += data_ptr[TENSOR_MAX_DIM * i + j] * data_ptr[TENSOR_MAX_DIM * i + j];
        }
    }
    return std::sqrt(tmp);
}

//determinant of KinematicTensor
double KinematicTensor::det(){
    if (DIM == 1){
        return data_ptr[XX];
    } else if (DIM == 2){
        return (data_ptr[XX]*data_ptr[YY] - data_ptr[XY]*data_ptr[YX]);
    } else {
        double a = data_ptr[XX] * (data_ptr[YY] * data_ptr[ZZ] - data_ptr[ZY] * data_ptr[YZ]);
        double b = data_ptr[XY] * (data_ptr[YX] * data_ptr[ZZ] - data_ptr[ZX] * data_ptr[YZ]);
        double c = data_ptr[XZ] * (data_ptr[YX] * data_ptr[ZY] - data_ptr[ZX] * data_ptr[YY]);
        return a-b+c;
    }
}

/*------------------------------------------------------------------------*/
//tensor inverse of KinematicTensor (not safe)
KinematicTensor KinematicTensor::inverse(){
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    double detA = det();
    if (DIM==1){
        tmp[XX] = 1/detA;
        for (int i=1;i<TENSOR_MAX_LENGTH;i++){
            tmp[i] = 0;
        }
    } else if (DIM==2){
        tmp[XX] = 1/detA * data_ptr[YY]; tmp[XY] = -1/detA * data_ptr[XY]; tmp[XZ] = 0;
        tmp[YX] = -1/detA * data_ptr[YX]; tmp[YY] = 1/detA * data_ptr[XX]; tmp[YZ] = 0;
        tmp[ZX] = 0; tmp[ZY] = 0; tmp[ZZ] = 0;
    } else {
        tmp[XX] = 1 / detA * (data_ptr[YY] * data_ptr[ZZ] - data_ptr[ZY] * data_ptr[YZ]);
        tmp[XY] = 1 / detA * (data_ptr[XZ] * data_ptr[ZY] - data_ptr[ZZ] * data_ptr[XY]);
        tmp[XZ] = 1 / detA * (data_ptr[XY] * data_ptr[YZ] - data_ptr[YY] * data_ptr[XZ]);
        tmp[YX] = 1 / detA * (data_ptr[YZ] * data_ptr[ZX] - data_ptr[ZZ] * data_ptr[YX]);
        tmp[YY] = 1 / detA * (data_ptr[XX] * data_ptr[ZZ] - data_ptr[ZX] * data_ptr[XZ]);
        tmp[YZ] = 1 / detA * (data_ptr[XZ] * data_ptr[YX] - data_ptr[YZ] * data_ptr[XX]);
        tmp[ZX] = 1 / detA * (data_ptr[YX] * data_ptr[ZY] - data_ptr[ZX] * data_ptr[YY]);
        tmp[ZY] = 1 / detA * (data_ptr[XY] * data_ptr[ZX] - data_ptr[ZY] * data_ptr[ZX]);
        tmp[ZZ] = 1 / detA * (data_ptr[XX] * data_ptr[YY] - data_ptr[YX] * data_ptr[XY]);
    }
    for (int i=0;i<TENSOR_MAX_LENGTH;i++){
        if (!std::isfinite(tmp[i])){
            std::cerr << "WARNING: MaterialTensor::inverse() has resulted in non-finite value. So that's bad." << std::endl;
        }
    }
    return KinematicTensor(tmp.data(),TENSOR_TYPE);
}

//tensor deviator (A_0 = A - 1/3 tr(A) I)
MaterialTensor KinematicTensor::deviator(){
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    double trA = trace();
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp[i] = data_ptr[i];
    }
    tmp[XX] -= trA/3.0; tmp[YY] -= trA/3.0; tmp[ZZ] -= trA/3.0;
    return MaterialTensor(tmp.data());
}

//transpose (A^T)
KinematicTensor KinematicTensor::transpose(){
    std::array<double,9> tmp;
    if (DIM == 1){
        return KinematicTensor(data_ptr,TENSOR_TYPE);
    } if (DIM == 2){
        tmp[XX] = data_ptr[XX]; tmp[XY] = data_ptr[YX]; tmp[XZ] = 0;
        tmp[YX] = data_ptr[XY]; tmp[YY] = data_ptr[YY]; tmp[YZ] = 0;
        tmp[ZX] = 0; tmp[ZY] = 0; tmp[ZZ] = 0;
    }
    tmp[XX] = data_ptr[XX]; tmp[XY] = data_ptr[YX]; tmp[XZ] = data_ptr[ZX];
    tmp[YX] = data_ptr[XY]; tmp[YY] = data_ptr[YY]; tmp[YZ] = data_ptr[ZY];
    tmp[ZX] = data_ptr[XZ]; tmp[ZY] = data_ptr[YZ]; tmp[ZZ] = data_ptr[ZZ];
    return KinematicTensor(tmp.data(), TENSOR_TYPE);
}

//symmetric and skew decompositions of A
KinematicTensor KinematicTensor::sym(){
    KinematicTensor tmp(this->data(),this->TENSOR_TYPE);
    return 0.5*(tmp + tmp.transpose());
}
KinematicTensor KinematicTensor::skw(){
    KinematicTensor tmp(this->data(),this->TENSOR_TYPE);
    return 0.5*(tmp - tmp.transpose());
}

/*----------------------------------------------------------------------------*/
//Tensor Addition
MaterialTensor operator+ (const MaterialTensor& lhs, const MaterialTensor& rhs){
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    return MaterialTensor(tmp.data());
}

MaterialTensor operator+ (const MaterialTensor& lhs, const KinematicTensor& rhs){
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    return MaterialTensor(tmp.data());
}

MaterialTensor operator+ (const KinematicTensor& lhs, const MaterialTensor& rhs){
    return rhs+lhs;
}

KinematicTensor operator+ (const KinematicTensor& lhs, const KinematicTensor& rhs){
    assert(lhs.TENSOR_TYPE == rhs.TENSOR_TYPE && "Addition failed.");
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    return KinematicTensor(tmp.data(),lhs.TENSOR_TYPE);
}

/*----------------------------------------------------------------------------*/
//Tensor Subtraction
MaterialTensor operator- (const MaterialTensor& lhs, const MaterialTensor& rhs){
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialTensor(tmp.data());
}

MaterialTensor operator- (const MaterialTensor& lhs, const KinematicTensor& rhs){
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialTensor(tmp.data());
}

MaterialTensor operator- (const KinematicTensor& lhs, const MaterialTensor& rhs){
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialTensor(tmp.data());
}

KinematicTensor operator- (const KinematicTensor& lhs, const KinematicTensor& rhs){
    assert(lhs.TENSOR_TYPE == rhs.TENSOR_TYPE && "Addition failed.");
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return KinematicTensor(tmp.data(), lhs.TENSOR_TYPE);
}

/*----------------------------------------------------------------------------*/
//Tensor Multiplication (returning MaterialTensor)
MaterialTensor operator* (const double& lhs, const MaterialTensor& rhs){
    MaterialTensor tmp = MaterialTensor(rhs);
    return tmp*=lhs;
}

MaterialTensor operator* (const int& lhs, const MaterialTensor& rhs){
    MaterialTensor tmp = MaterialTensor(rhs);
    return tmp*=lhs;
}

MaterialTensor operator* (const MaterialTensor& lhs, const double& rhs){
    MaterialTensor tmp = MaterialTensor(lhs);
    return tmp*=rhs;
}

MaterialTensor operator* (const MaterialTensor& lhs, const int& rhs){
    MaterialTensor tmp = MaterialTensor(lhs);
    return tmp*=rhs;
}

MaterialTensor operator/ (const MaterialTensor& lhs, const double& rhs){
    MaterialTensor tmp = MaterialTensor(lhs);
    return tmp/=rhs;
}

MaterialTensor operator/ (const MaterialTensor& lhs, const int& rhs){
    MaterialTensor tmp = MaterialTensor(lhs);
    return tmp/=rhs;
}

MaterialTensor operator* (const MaterialTensor& lhs, const MaterialTensor& rhs){
    std::array<double,MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_DIM;i++){
        for(int j=0;j<MPMTensor::TENSOR_MAX_DIM;j++) {
            tmp[MPMTensor::TENSOR_MAX_DIM*i+j] = lhs[MPMTensor::TENSOR_MAX_DIM * i + 0]*rhs[MPMTensor::TENSOR_MAX_DIM * 0 + j] +
                                      lhs[MPMTensor::TENSOR_MAX_DIM * i + 1]*rhs[MPMTensor::TENSOR_MAX_DIM * 1 + j] +
                                      lhs[MPMTensor::TENSOR_MAX_DIM * i + 2]*rhs[MPMTensor::TENSOR_MAX_DIM * 2 + j];
        }
    }
    return MaterialTensor(tmp.data());
}

MaterialTensor operator* (const MaterialTensor& lhs, const KinematicTensor& rhs){
    std::array<double,MPMTensor::TENSOR_MAX_LENGTH> tmp;
    if (rhs.DIM == MPMTensor::TENSOR_MAX_DIM){
        for(int i=0;i<MPMTensor::TENSOR_MAX_DIM;i++){
            for(int j=0;j<MPMTensor::TENSOR_MAX_DIM;j++) {
                tmp[MPMTensor::TENSOR_MAX_DIM*i+j] = lhs[MPMTensor::TENSOR_MAX_DIM * i + 0]*rhs[MPMTensor::TENSOR_MAX_DIM * 0 + j] +
                                          lhs[MPMTensor::TENSOR_MAX_DIM * i + 1]*rhs[MPMTensor::TENSOR_MAX_DIM * 1 + j] +
                                          lhs[MPMTensor::TENSOR_MAX_DIM * i + 2]*rhs[MPMTensor::TENSOR_MAX_DIM * 2 + j];
            }
        }
    } else {
        for (int i = 0; i < MPMTensor::TENSOR_MAX_DIM; i++) {
            for (int j = 0; j < rhs.DIM; j++) {
                tmp[MPMTensor::TENSOR_MAX_DIM * i + j] = 0;
                for (int k = 0; k < rhs.DIM; k++) {
                    tmp[MPMTensor::TENSOR_MAX_DIM * i + j] += lhs[MPMTensor::TENSOR_MAX_DIM * i + k] * rhs[3 * k + j];
                }
            }
        }
        for (int j=rhs.DIM;j<MPMTensor::TENSOR_MAX_DIM;j++){
            tmp[MPMTensor::TENSOR_MAX_DIM * 0 + j] = 0;
            tmp[MPMTensor::TENSOR_MAX_DIM * 1 + j] = 0;
            tmp[MPMTensor::TENSOR_MAX_DIM * 2 + j] = 0;
        }
    }
    return MaterialTensor(tmp.data());
}

MaterialTensor operator* (const KinematicTensor& lhs, const MaterialTensor& rhs){
    std::array<double,MPMTensor::TENSOR_MAX_LENGTH> tmp;
    if (lhs.DIM == MPMTensor::TENSOR_MAX_DIM){
        for(int i=0;i<MPMTensor::TENSOR_MAX_DIM;i++){
            for(int j=0;j<MPMTensor::TENSOR_MAX_DIM;j++) {
                tmp[MPMTensor::TENSOR_MAX_DIM*i+j] = lhs[MPMTensor::TENSOR_MAX_DIM*i+0]*rhs[MPMTensor::TENSOR_MAX_DIM *0+j] +
                                          lhs[MPMTensor::TENSOR_MAX_DIM*i+1]*rhs[MPMTensor::TENSOR_MAX_DIM*1+j] +
                                          lhs[MPMTensor::TENSOR_MAX_DIM*i+2]*rhs[MPMTensor::TENSOR_MAX_DIM*2+j];
            }
        }
    } else {
        for (int i = 0; i < lhs.DIM; i++) {
            for (int j = 0; j < MPMTensor::TENSOR_MAX_DIM; j++) {
                tmp[MPMTensor::TENSOR_MAX_DIM * i + j] = 0;
                for (int k = 0; k < lhs.DIM; k++) {
                    tmp[MPMTensor::TENSOR_MAX_DIM * i + j] += lhs[MPMTensor::TENSOR_MAX_DIM * i + k] * rhs[MPMTensor::TENSOR_MAX_DIM * k + j];
                }
            }
        }
        for (int i=lhs.DIM;i<MPMTensor::TENSOR_MAX_DIM;i++){
            tmp[MPMTensor::TENSOR_MAX_DIM*i + 0] = 0;
            tmp[MPMTensor::TENSOR_MAX_DIM*i + 1] = 0;
            tmp[MPMTensor::TENSOR_MAX_DIM*i + 2] = 0;
        }
    }
    return MaterialTensor(tmp.data());
}


/*----------------------------------------------------------------------------*/
//Tensor Multiplication (returning KinematicTensor)
KinematicTensor operator* (const double& lhs, const KinematicTensor& rhs){
    KinematicTensor tmp = KinematicTensor(rhs);
    return tmp*=lhs;
}

KinematicTensor operator* (const int& lhs, const KinematicTensor& rhs){
    KinematicTensor tmp = KinematicTensor(rhs);
    return tmp*=lhs;
}

KinematicTensor operator* (const KinematicTensor& lhs, const double& rhs){
    return rhs*lhs;
}

KinematicTensor operator* (const KinematicTensor& lhs, const int& rhs){
    return rhs*lhs;
}

KinematicTensor operator/ (const KinematicTensor& lhs, const double& rhs){
    KinematicTensor tmp = KinematicTensor(lhs);
    return tmp/=rhs;
}

KinematicTensor operator/ (const KinematicTensor& lhs, const int& rhs){
    KinematicTensor tmp = KinematicTensor(lhs);
    return tmp/=rhs;
}

KinematicTensor operator* (const KinematicTensor& lhs, const KinematicTensor& rhs){
    assert(lhs.TENSOR_TYPE == rhs.TENSOR_TYPE && "Multiplication failed.");
    std::array<double,MPMTensor::TENSOR_MAX_LENGTH> tmp;
    if (lhs.DIM == MPMTensor::TENSOR_MAX_DIM){
        for(int i=0;i<MPMTensor::TENSOR_MAX_DIM;i++){
            for(int j=0;j<MPMTensor::TENSOR_MAX_DIM;j++) {
                tmp[MPMTensor::TENSOR_MAX_DIM*i+j] = lhs[MPMTensor::TENSOR_MAX_DIM*i+0]*rhs[MPMTensor::TENSOR_MAX_DIM *0+j] +
                                          lhs[MPMTensor::TENSOR_MAX_DIM*i+1]*rhs[MPMTensor::TENSOR_MAX_DIM*1+j] +
                                          lhs[MPMTensor::TENSOR_MAX_DIM*i+2]*rhs[MPMTensor::TENSOR_MAX_DIM*2+j];
            }
        }
    } else {
        for (int i = 0; i < lhs.DIM; i++) {
            for (int j = 0; j < MPMTensor::TENSOR_MAX_DIM; j++) {
                tmp[MPMTensor::TENSOR_MAX_DIM * i + j] = 0;
                for (int k = 0; k < lhs.DIM; k++) {
                    tmp[MPMTensor::TENSOR_MAX_DIM * i + j] += lhs[MPMTensor::TENSOR_MAX_DIM * i + k] * rhs[MPMTensor::TENSOR_MAX_DIM * k + j];
                }
            }
        }
        for (int i=lhs.DIM;i<MPMTensor::TENSOR_MAX_DIM;i++){
            tmp[MPMTensor::TENSOR_MAX_DIM*i + 0] = 0;
            tmp[MPMTensor::TENSOR_MAX_DIM*i + 1] = 0;
            tmp[MPMTensor::TENSOR_MAX_DIM*i + 2] = 0;
        }
    }
    return KinematicTensor(tmp.data(),lhs.TENSOR_TYPE);
}



/*----------------------------------------------------------------------------*/
//Tensor Maps
MaterialTensor::Map& MaterialTensor::Map::operator= (const MaterialTensor &other){
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){data_ptr[i]=other(i);};
    return *this;
}

MaterialTensor::Map& MaterialTensor::Map::operator= (const KinematicTensor &other){
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){data_ptr[i]=other(i);};
    return *this;
}

KinematicTensor::Map& KinematicTensor::Map::operator= (const KinematicTensor &other){
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){data_ptr[i]=other(i);};
    return *this;
}
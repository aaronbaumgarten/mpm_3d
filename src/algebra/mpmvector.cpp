//
// Created by aaron on 4/23/18.
// mpmvector.cpp
//

#include "mpm_tensor.hpp"
#include "mpm_vector.hpp"
#include <math.h>
#include <iostream>

/*----------------------------------------------------------------------------*/
//construct MaterialVector with pointer to data (not safe, but who cares right?)
MaterialVector::MaterialVector(double* otherdata){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] = otherdata[i];
    }
}

//construct MaterialVector from KinematicVector
MaterialVector::MaterialVector(KinematicVector& other){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] = other[i];
    }
}

/*----------------------------------------------------------------------------*/
//define -v
MaterialVector MaterialVector::operator-(){
    std::array<double, VECTOR_MAX_DIM> tmp;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        tmp[i] = -data_ptr[i];
    }
    return MaterialVector(tmp.data());
}

//define v *= s for integer scalar value
MaterialVector& MaterialVector::operator*= (const int &rhs){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] *= rhs;
    }
    return *this;
}

//define v *= s for double scalar
MaterialVector& MaterialVector::operator*= (const double &rhs){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] *= rhs;
    }
    return *this;
}

//define v /= s for integer scalar value
MaterialVector& MaterialVector::operator/= (const int &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] *= tmp;
    }
    return *this;
}

//define v /= s for double scalar
MaterialVector& MaterialVector::operator/= (const double &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] *= tmp;
    }
    return *this;
}

/*----------------------------------------------------------------------------*/
//set MaterialVector to standard vectors
void MaterialVector::setZero(){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] = 0;
    }
    return;
}
void MaterialVector::setOnes(){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] = 1;
    }
    return;
}

/*----------------------------------------------------------------------------*/
//vector inner product
double MaterialVector::dot(const MaterialVector &rhs){
    double tmp = 0;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        tmp += data_ptr[i]*rhs[i];
    }
    return tmp;
}

double MaterialVector::dot(const KinematicVector &rhs){
    double tmp = 0;
    for(int i=0;i<rhs.DIM;i++){
        tmp += data_ptr[i]*rhs[i];
    }
    return tmp;
}

//vector cross product
MaterialVector MaterialVector::cross(const MaterialVector &rhs){
    std::array<double, VECTOR_MAX_DIM> tmp;
    tmp[X] = data_ptr[Y]*rhs[Z] - data_ptr[Z]*rhs[Y];
    tmp[Y] = data_ptr[Z]*rhs[X] - data_ptr[X]*rhs[Z];
    tmp[Z] = data_ptr[X]*rhs[Y] - data_ptr[Y]*rhs[X];
    return MaterialVector(tmp.data());
}

MaterialVector MaterialVector::cross(const KinematicVector &rhs){
    std::array<double, VECTOR_MAX_DIM> tmp;
    tmp[X] = data_ptr[Y]*rhs[Z] - data_ptr[Z]*rhs[Y];
    tmp[Y] = data_ptr[Z]*rhs[X] - data_ptr[X]*rhs[Z];
    tmp[Z] = data_ptr[X]*rhs[Y] - data_ptr[Y]*rhs[X];
    return MaterialVector(tmp.data());
}

//vector tensor product
MaterialTensor MaterialVector::tensor(const MaterialVector& rhs){
    std::array<double,MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        for(int j=0;j<VECTOR_MAX_DIM;j++){
            tmp[VECTOR_MAX_DIM * i + j] = data_ptr[i]*rhs[j];
        }
    }
    return MaterialTensor(tmp.data());
}
MaterialTensor MaterialVector::tensor(const KinematicVector& rhs){
    std::array<double,MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        for(int j=0;j<rhs.DIM;j++){
            tmp[VECTOR_MAX_DIM * i + j] = data_ptr[i]*rhs[j];
        }
        for(int j=rhs.DIM;j<VECTOR_MAX_DIM;j++){
            tmp[VECTOR_MAX_DIM * i + j] = 0;
        }
    }
    return MaterialTensor(tmp.data());
}

//2-norm of MaterialVector
double MaterialVector::norm(){
    double tmp = 0;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        tmp += data_ptr[i]*data_ptr[i];
    }
    return std::sqrt(tmp);
}


/*----------------------------------------------------------------------------*/
//construct KinematicVector from pointer to data and input type
KinematicVector::KinematicVector(double* otherdata, int input){
    assignVectorType(input);
    for(int i=0;i<DIM;i++){
        data_ptr[i] = otherdata[i];
    }
    for(int i=DIM;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] = 0;
    }
}

//construct KinematicVector from MaterialVector
KinematicVector::KinematicVector(MaterialVector& other, int input){
    assignVectorType(input);
    for(int i=0;i<DIM;i++){
        data_ptr[i] = other[i];
    }
    for(int i=DIM;i<VECTOR_MAX_DIM;i++){
        if (other[i] != 0){
            std::cerr << "WARNING: ignoring non-zero entries in input data to KinematicVector." << std::endl;
        }
        data_ptr[i] = 0;
    }
}

/*----------------------------------------------------------------------------*/
//define -v
KinematicVector KinematicVector::operator-(){
    std::array<double, VECTOR_MAX_DIM> tmp;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        tmp[i] = -data_ptr[i];
    }
    return KinematicVector(tmp.data(),VECTOR_TYPE);
}

//define v *= s for integer scalar value
KinematicVector& KinematicVector::operator*= (const int &rhs){
    for(int i=0;i<DIM;i++){
        data_ptr[i] *= rhs;
    }
    return *this;
}

//define v *= s for double scalar
KinematicVector& KinematicVector::operator*= (const double &rhs){
    for(int i=0;i<DIM;i++){
        data_ptr[i] *= rhs;
    }
    return *this;
}

//define v /= s for integer scalar value
KinematicVector& KinematicVector::operator/= (const int &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<DIM;i++){
        data_ptr[i] *= tmp;
    }
    return *this;
}

//define v /= s for double scalar
KinematicVector& KinematicVector::operator/= (const double &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<DIM;i++){
        data_ptr[i] *= tmp;
    }
    return *this;
}

/*----------------------------------------------------------------------------*/
//set KinematicVector to standard vectors
void KinematicVector::setZero(){
    for(int i=0;i<DIM;i++){
        data_ptr[i] = 0;
    }
    return;
}

void KinematicVector::setOnes(){
    for(int i=0;i<DIM;i++){
        data_ptr[i] = 1;
    }
    return;
}

/*----------------------------------------------------------------------------*/
//vector contraction
double KinematicVector::dot(const MaterialVector &rhs){
    double tmp = 0;
    for(int i=0;i<DIM;i++){
        tmp += data_ptr[i]*rhs[i];
    }
    return tmp;
}

double KinematicVector::dot(const KinematicVector &rhs){
    double tmp = 0;
    for(int i=0;i<DIM;i++){
        tmp += data_ptr[i]*rhs[i];
    }
    return tmp;
}

//vector cross product (always return MaterialVector)
MaterialVector KinematicVector::cross(const MaterialVector &rhs){
    std::array<double, VECTOR_MAX_DIM> tmp;
    tmp[X] = data_ptr[Y]*rhs[Z] - data_ptr[Z]*rhs[Y];
    tmp[Y] = data_ptr[Z]*rhs[X] - data_ptr[X]*rhs[Z];
    tmp[Z] = data_ptr[X]*rhs[Y] - data_ptr[Y]*rhs[X];
    return MaterialVector(tmp.data());
}

MaterialVector KinematicVector::cross(const KinematicVector &rhs){
    std::array<double, VECTOR_MAX_DIM> tmp;
    tmp[X] = data_ptr[Y]*rhs[Z] - data_ptr[Z]*rhs[Y];
    tmp[Y] = data_ptr[Z]*rhs[X] - data_ptr[X]*rhs[Z];
    tmp[Z] = data_ptr[X]*rhs[Y] - data_ptr[Y]*rhs[X];
    return MaterialVector(tmp.data());
}

//vector tensor product
MaterialTensor KinematicVector::tensor(const MaterialVector& rhs){
    std::array<double,MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<DIM;i++){
        for(int j=0;j<VECTOR_MAX_DIM;j++){
            tmp[VECTOR_MAX_DIM * i + j] = data_ptr[i]*rhs[j];
        }
    }
    for (int i=DIM;i<VECTOR_MAX_DIM;i++){
        for (int j=0;j<VECTOR_MAX_DIM;j++){
            tmp[VECTOR_MAX_DIM * i + j] = 0;
        }
    }
    return MaterialTensor(tmp.data());
}

KinematicTensor KinematicVector::tensor(const KinematicVector& rhs){
    assert(VECTOR_TYPE == rhs.VECTOR_TYPE && "Tensor product failed.");
    std::array<double,MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<DIM;i++){
        for(int j=0;j<rhs.DIM;j++){
            tmp[VECTOR_MAX_DIM * i + j] = data_ptr[i]*rhs[j];
        }
        for(int j=rhs.DIM;j<VECTOR_MAX_DIM;j++){
            tmp[VECTOR_MAX_DIM * i + j] = 0;
        }
    }
    for (int i=DIM;i<VECTOR_MAX_DIM;i++){
        for (int j=0;j<VECTOR_MAX_DIM;j++){
            tmp[VECTOR_MAX_DIM * i + j] = 0;
        }
    }
    return KinematicTensor(tmp.data(),VECTOR_TYPE);
}

//vector 2-norm
double KinematicVector::norm(){
    double tmp = 0;
    for(int i=0;i<DIM;i++){
        tmp += data_ptr[i]*data_ptr[i];
    }
    return tmp;
}


/*----------------------------------------------------------------------------*/
//Vector Addition
MaterialVector operator+ (const MaterialVector& rhs, const MaterialVector& lhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    return MaterialVector(tmp.data());
}

MaterialVector operator+ (const MaterialVector& rhs, const KinematicVector& lhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    return MaterialVector(tmp.data());
}

MaterialVector operator+ (const KinematicVector& lhs, const MaterialVector& rhs){
    return (rhs+lhs);
}

KinematicVector operator+ (const KinematicVector& lhs, const KinematicVector& rhs){
    assert(lhs.VECTOR_TYPE == rhs.VECTOR_TYPE && "Addition failed.");
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<lhs.DIM;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    for(int i=lhs.DIM;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = 0;
    }
    return KinematicVector(tmp.data(), lhs.VECTOR_TYPE);
}

/*----------------------------------------------------------------------------*/
//Vector Subtraction
MaterialVector operator- (const MaterialVector& lhs, const MaterialVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialVector(tmp.data());
}

MaterialVector operator- (const MaterialVector& lhs, const KinematicVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialVector(tmp.data());
}

MaterialVector operator- (const KinematicVector& lhs, const MaterialVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialVector(tmp.data());
}

KinematicVector operator- (const KinematicVector& lhs, const KinematicVector& rhs){
    assert(lhs.VECTOR_TYPE == rhs.VECTOR_TYPE && "Subtraction failed.");
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<lhs.DIM;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    for(int i=lhs.DIM;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = 0;
    }
    return KinematicVector(tmp.data(), lhs.VECTOR_TYPE);
}


/*----------------------------------------------------------------------------*/
//Multiplication (returning MaterialVector)
MaterialVector operator* (const double& lhs, const MaterialVector& rhs){
    MaterialVector tmp = MaterialVector(rhs);
    return tmp*=lhs;
}

MaterialVector operator* (const int& lhs, const MaterialVector& rhs){
    MaterialVector tmp = MaterialVector(rhs);
    return tmp*=lhs;
}

MaterialVector operator* (const MaterialVector& lhs, const double& rhs){
    MaterialVector tmp = MaterialVector(lhs);
    return tmp*=rhs;
}

MaterialVector operator* (const MaterialVector& lhs, const int& rhs){
    MaterialVector tmp = MaterialVector(lhs);
    return tmp*=rhs;
}

MaterialVector operator* (const MaterialTensor& lhs, const MaterialVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    tmp[MPMVector::X] = lhs[MPMTensor::XX]*rhs[MPMVector::X] + lhs[MPMTensor::XY]*rhs[MPMVector::Y] + lhs[MPMTensor::XZ]*rhs[MPMVector::Z];
    tmp[MPMVector::Y] = lhs[MPMTensor::YX]*rhs[MPMVector::X] + lhs[MPMTensor::YY]*rhs[MPMVector::Y] + lhs[MPMTensor::YZ]*rhs[MPMVector::Z];
    tmp[MPMVector::Z] = lhs[MPMTensor::ZX]*rhs[MPMVector::X] + lhs[MPMTensor::ZY]*rhs[MPMVector::Y] + lhs[MPMTensor::ZZ]*rhs[MPMVector::Z];
    return MaterialVector(tmp.data());
}

MaterialVector operator* (const MaterialTensor& lhs, const KinematicVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp = { {0} };
    for (int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        for (int j=0;j<rhs.DIM;j++){
            tmp[i] += lhs(i,j) * rhs(j);
        }
    }
    return MaterialVector(tmp.data());
}

MaterialVector operator* (const KinematicTensor& lhs, const MaterialVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp = { {0} };
    for (int i=0;i<lhs.DIM;i++){
        for (int j=0;j<MPMVector::VECTOR_MAX_DIM;j++){
            tmp[i] += lhs(i,j) * rhs(j);
        }
    }
    return MaterialVector(tmp.data());
}

/*----------------------------------------------------------------------------*/
//Division (returning MaterialVector)
MaterialVector operator/ (const MaterialVector& lhs, const double& rhs){
    MaterialVector tmp = MaterialVector(lhs);
    return tmp/=rhs;
}

MaterialVector operator/ (const MaterialVector& lhs, const int& rhs){
    MaterialVector tmp = MaterialVector(lhs);
    return tmp/=rhs;
}

/*----------------------------------------------------------------------------*/
//Multiplication (returning KinematicTensor)
KinematicVector operator* (const double& lhs, const KinematicVector& rhs){
    KinematicVector tmp = KinematicVector(rhs);
    return tmp*=lhs;
}

KinematicVector operator* (const int& lhs, const KinematicVector& rhs){
    KinematicVector tmp = KinematicVector(rhs);
    return tmp*=lhs;
}

KinematicVector operator* (const KinematicVector& lhs, const double& rhs){
    KinematicVector tmp = KinematicVector(lhs);
    return tmp*=rhs;
}

KinematicVector operator* (const KinematicVector& lhs, const int& rhs){
    KinematicVector tmp = KinematicVector(lhs);
    return tmp*=rhs;
}

KinematicVector operator* (const KinematicTensor& lhs, const KinematicVector& rhs){
    assert(lhs.TENSOR_TYPE == rhs.VECTOR_TYPE && "Multiplication failed.");
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp = { {0} };
    for (int i=0;i<lhs.DIM;i++){
        for (int j=0;j<rhs.DIM;j++){
            tmp[i] += lhs(i,j) * rhs(j);
        }
    }
    return KinematicVector(tmp.data(), rhs.VECTOR_TYPE);
}

KinematicVector operator/ (const KinematicVector& lhs, const double& rhs){
    KinematicVector tmp = KinematicVector(lhs);
    return tmp/=rhs;
}
KinematicVector operator/ (const KinematicVector& lhs, const int& rhs){
    KinematicVector tmp = KinematicVector(lhs);
    return tmp/=rhs;
}


/*----------------------------------------------------------------------------*/
//Vector Maps
MaterialVector::Map& MaterialVector::Map::operator= (const MaterialVector &other){
    for(int i=0;i<VECTOR_MAX_DIM;i++){data_ptr[i]=other(i);};
    return *this;
}

MaterialVector::Map& MaterialVector::Map::operator= (const KinematicVector &other){
    for(int i=0;i<VECTOR_MAX_DIM;i++){data_ptr[i]=other(i);};
    return *this;
}

KinematicVector::Map& KinematicVector::Map::operator= (const KinematicVector &other){
    for(int i=0;i<VECTOR_MAX_DIM;i++){data_ptr[i]=other(i);};
    return *this;
}
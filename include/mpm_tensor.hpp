//
// Created by aaron on 4/22/18.
// algebra.hpp
//

#ifndef MPM_V2_TENSOR_HPP
#define MPM_V2_TENSOR_HPP

#include <stdlib.h>
#include <array>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <array>
#include "factorial.hpp"

/*----------------------------------------------------------------------------*/
//Eigen::Map object which interfaces with KinematicTensor object.

typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>,0,Eigen::Stride<3,1>> EIGEN_MAP_OF_KINEMATIC_TENSOR;

/*----------------------------------------------------------------------------*/
//Eigen::Map object which interfaces with MaterialTensor object.

typedef Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor>> EIGEN_MAP_OF_MATERIAL_TENSOR;

/*----------------------------------------------------------------------------*/
//base class for KinematicTensor and MaterialTensor.
class MPMTensor {
public:
    MPMTensor(){
        buffer = { {0} };
        data_ptr = buffer.data();
    }

    MPMTensor(const MPMTensor& other){
        for (int i=0;i<TENSOR_MAX_LENGTH;i++){
            buffer[i] = other.data_ptr[i];
        }
        data_ptr = buffer.data(); //make sure not to copy stale pointer
    }

    MPMTensor(const MPMTensor&& other){
        for (int i=0;i<TENSOR_MAX_LENGTH;i++){
            buffer[i] = other.data_ptr[i];
        }
        data_ptr = buffer.data(); //make sure not to copy stale pointer
    }

    MPMTensor& operator= (const MPMTensor& other){
        for (int i=0;i<TENSOR_MAX_LENGTH;i++){
            buffer[i] = other.data_ptr[i];
        }
        data_ptr = buffer.data(); //make sure not to copy stale pointer
        return *this;
    }

    MPMTensor& operator= (const MPMTensor&& other){
        for (int i=0;i<TENSOR_MAX_LENGTH;i++){
            buffer[i] = other.data_ptr[i];
        }
        data_ptr = buffer.data(); //make sure not to copy stale pointer
        return *this;
    }

    //set standard tensor length
    static const int TENSOR_MAX_DIM = 3, TENSOR_MAX_LENGTH = TENSOR_MAX_DIM*TENSOR_MAX_DIM;

    //set standard tensor access variables
    static const int XX = 0, XY = 1, XZ = 2,
                     YX = 3, YY = 4, YZ = 5,
                     ZX = 6, ZY = 7, ZZ = 8;

    //set standard vector access variables
    static const int X = 0,  Y = 1,  Z = 2;

    /*------------------------------------------------------------------------*/
    //functions for accessing tensor data
    double& operator() (int i, int j) const {
        return data_ptr[TENSOR_MAX_DIM*i + j];
    }

    double& operator() (int i) const {
        return data_ptr[i];
    }

    double& operator[] (int i) const {
        return data_ptr[i];
    }

    /*------------------------------------------------------------------------*/
    //return pointer to data
    double* data() const {
        return data_ptr;
    }

    /*------------------------------------------------------------------------*/
    //return size of tensor
    double size() const {
        return TENSOR_MAX_LENGTH;
    }

    /*------------------------------------------------------------------------*/
    //return rows and columns
    virtual double rows() const {
        return TENSOR_MAX_DIM;
    }

    virtual double cols() const {
        return TENSOR_MAX_DIM;
    }

    /*------------------------------------------------------------------------*/
private:
    //standard buffer for data
    std::array<double, TENSOR_MAX_LENGTH> buffer;

    /*------------------------------------------------------------------------*/
protected:
    //pointer to data (by default points to standard buffer)
    double* data_ptr;
};


/*----------------------------------------------------------------------------*/
class KinematicTensor;

/*----------------------------------------------------------------------------*/
class MaterialTensor : public MPMTensor{
public:
    //default constructor
    MaterialTensor(){}

    //copy constructor
    MaterialTensor(const MaterialTensor& other) : MPMTensor(other){};

    //construct tensor from data pointer
    MaterialTensor(const double* otherdata);

    //copy from KinematicTensor
    MaterialTensor(const KinematicTensor&);

    //forward declare Map class
    class Map;

    //construct from Eigen::Matrix
    template <typename OtherDerived>
    MaterialTensor(Eigen::MatrixBase<OtherDerived> &other){
        assert(other.rows() == TENSOR_MAX_DIM && other.cols() == TENSOR_MAX_DIM);
        for(int i=0;i<TENSOR_MAX_DIM;i++){
            for(int j=0;j<TENSOR_MAX_DIM;j++){
                data_ptr[TENSOR_MAX_DIM*i+j] = other(i,j);
            }
        }
    }

    /*------------------------------------------------------------------------*/
    //define operators
    operator EIGEN_MAP_OF_MATERIAL_TENSOR () const {return EIGEN_MAP_OF_MATERIAL_TENSOR(data_ptr);}
    MaterialTensor operator-();
    inline MaterialTensor& operator+= (const MaterialTensor &rhs);
    inline MaterialTensor& operator+= (const KinematicTensor &rhs);
    inline MaterialTensor& operator-= (const MaterialTensor &rhs);
    inline MaterialTensor& operator-= (const KinematicTensor &rhs);
    MaterialTensor& operator*= (const int &rhs);
    MaterialTensor& operator*= (const double &rhs);
    MaterialTensor& operator/= (const int &rhs);
    MaterialTensor& operator/= (const double &rhs);

    /*------------------------------------------------------------------------*/
    //set to standard tensors
    void setZero();
    void setIdentity();
    static MaterialTensor Identity();

    /*------------------------------------------------------------------------*/
    //tensor contractions
    double dot(const MaterialTensor &rhs) const;
    double dot(const KinematicTensor &rhs) const;

    //tensor trace
    double trace() const;

    //tensor frobenius norm
    double norm() const;

    //tensor determinant
    double det() const;

    /*------------------------------------------------------------------------*/
    //inverse
    MaterialTensor inverse() const;

    //deviator
    MaterialTensor deviator() const;

    //transpose
    MaterialTensor transpose() const;

    //symetric and skew decomposition
    MaterialTensor sym() const;
    MaterialTensor skw() const;

    //tensor exponent
    MaterialTensor exp(int order = 5) const;
};

/*----------------------------------------------------------------------------*/
class MaterialTensor::Map : public MaterialTensor{
public:
    //point to other data (not safe)
    Map(const double* dataIN){ data_ptr = (double*)(dataIN); }

    //assignment operator
    Map& operator= (const MaterialTensor &other);
    Map& operator= (const KinematicTensor &other);
};

/*----------------------------------------------------------------------------*/
class KinematicTensor : public MPMTensor{
public:
    //KinematicTensor types
    static const int TENSOR_1D = 1, TENSOR_2D = 2, TENSOR_3D = 3, TENSOR_AXISYM = 4;
    static const int TENSOR_2D_OOP = 5;

    //default dimension
    int DIM = TENSOR_MAX_DIM;
    int TENSOR_TYPE = TENSOR_3D;

    void assignTensorType(int input){
        TENSOR_TYPE = input;
        if (TENSOR_TYPE == TENSOR_1D){
            DIM = 1;
        } else if (TENSOR_TYPE == TENSOR_2D){
            DIM = 2;
        } else if (TENSOR_TYPE == TENSOR_3D){
            DIM = 3;
        } else if (TENSOR_TYPE == TENSOR_2D_OOP){
            DIM = 3;
        } else if (TENSOR_TYPE == TENSOR_AXISYM){
            DIM = 3;
        } else {
            std::cerr << "KinematicTensor doesn't have defined type for input " << TENSOR_TYPE << "." << std::endl;
        }
    }

    /*------------------------------------------------------------------------*/
    //return rows and columns
    double rows() const {
        return DIM;
    }

    double cols() const {
        return DIM;
    }

    //default constructor
    KinematicTensor() : MPMTensor(){
        DIM = TENSOR_MAX_DIM;
        TENSOR_TYPE = TENSOR_3D;
    };

    //copy constructor
    KinematicTensor(const KinematicTensor& other) : MPMTensor(other){
        DIM = other.DIM;
        TENSOR_TYPE = other.TENSOR_TYPE;
    }

    KinematicTensor(int input){
        assignTensorType(input);
    }

    //construct tensor from data pointer
    KinematicTensor(const double* otherdata, int input);

    //copy from KinematicTensor
    //KinematicTensor(KinematicTensor&);
    //KinematicTensor(const KinematicTensor&, int input);

    //copy from MaterialTensor
    KinematicTensor(const MaterialTensor&, int input);

    //forward declare Map class
    class Map;

    //construct from Eigen::Matrix
    template <typename OtherDerived>
    KinematicTensor(const Eigen::MatrixBase<OtherDerived> &other, int input){
        assignTensorType(input);
        assert(other.rows() == DIM && other.cols() == DIM);
        for(int i=0;i<DIM;i++){
            for(int j=0;j<DIM;j++){
                data_ptr[TENSOR_MAX_DIM * i + j] = other(i,j);
            }
        }
        for (int i=DIM;i<3;i++) {
            for (int j = DIM; j < 3; j++) {
                data_ptr[TENSOR_MAX_DIM * i + j] = 0;
            }
        }
    }

    /*------------------------------------------------------------------------*/
    //define operators
    operator EIGEN_MAP_OF_KINEMATIC_TENSOR () const {return EIGEN_MAP_OF_KINEMATIC_TENSOR(data_ptr,DIM,DIM);}
    KinematicTensor operator-();
    inline KinematicTensor& operator+= (const KinematicTensor &rhs);
    inline KinematicTensor& operator-= (const KinematicTensor &rhs);
    KinematicTensor& operator*= (const int &rhs);
    KinematicTensor& operator*= (const double &rhs);
    KinematicTensor& operator/= (const int &rhs);
    KinematicTensor& operator/= (const double &rhs);

    /*------------------------------------------------------------------------*/
    //set to standard tensors
    void setZero();
    void setIdentity();
    static KinematicTensor Identity(int input);

    /*------------------------------------------------------------------------*/
    //tensor contractions
    double dot(const MaterialTensor &rhs) const;
    double dot(const KinematicTensor &rhs) const;

    //tensor trace
    double trace() const;

    //tensor frobenius norm
    double norm() const;

    //tensor determinant
    double det() const;

    /*------------------------------------------------------------------------*/
    //inverse
    KinematicTensor inverse() const;

    //deviator
    MaterialTensor deviator() const;

    //transpose
    KinematicTensor transpose() const;

    //symetric and skew decomposition
    KinematicTensor sym() const;
    KinematicTensor skw() const;

    //tensor exponent
    KinematicTensor exp(int order = 5) const;
};

/*----------------------------------------------------------------------------*/
class KinematicTensor::Map : public KinematicTensor{
public:
    //point to other data (not safe)
    Map(const double* dataIN, int input){
        assignTensorType(input);
        data_ptr = (double*)(dataIN);
    }

    //assignment operator
    Map& operator= (const KinematicTensor &other);
};

/*----------------------------------------------------------------------------*/
//operations which return a material tensor
inline MaterialTensor operator+ (const MaterialTensor&, const MaterialTensor&);
inline MaterialTensor operator+ (const MaterialTensor&, const KinematicTensor&);
inline MaterialTensor operator+ (const KinematicTensor&, const MaterialTensor&);
inline MaterialTensor operator- (const MaterialTensor&, const MaterialTensor&);
inline MaterialTensor operator- (const MaterialTensor&, const KinematicTensor&);
inline MaterialTensor operator- (const KinematicTensor&, const MaterialTensor&);
inline MaterialTensor operator* (const double&, const MaterialTensor&);
inline MaterialTensor operator* (const int&, const MaterialTensor&);
inline MaterialTensor operator* (const MaterialTensor&, const double&);
inline MaterialTensor operator* (const MaterialTensor&, const int&);
inline MaterialTensor operator* (const MaterialTensor&, const MaterialTensor&);
inline MaterialTensor operator* (const MaterialTensor&, const KinematicTensor&);
inline MaterialTensor operator* (const KinematicTensor&, const MaterialTensor&);
inline MaterialTensor operator/ (const MaterialTensor&, const double&);
inline MaterialTensor operator/ (const MaterialTensor&, const int&);

/*----------------------------------------------------------------------------*/
//operations which return a material tensor
inline KinematicTensor operator+ (const KinematicTensor&, const KinematicTensor&);
inline KinematicTensor operator- (const KinematicTensor&, const KinematicTensor&);
inline KinematicTensor operator* (const double&, const KinematicTensor&);
inline KinematicTensor operator* (const int&, const KinematicTensor&);
inline KinematicTensor operator* (const KinematicTensor&, const double&);
inline KinematicTensor operator* (const KinematicTensor&, const int&);
inline KinematicTensor operator* (const KinematicTensor&, const KinematicTensor&);
inline KinematicTensor operator/ (const KinematicTensor&, const double&);
inline KinematicTensor operator/ (const KinematicTensor&, const int&);




/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/



/*------------------------------------------------------------------------*/
//construct MaterialTensor with pointer to data (not safe, but who cares right?)
inline MaterialTensor::MaterialTensor(const double* otherdata){
    for (int i=0; i<TENSOR_MAX_LENGTH; i++){
        data_ptr[i] = otherdata[i];
    }
}


//construct MaterialTensor from a KinematicTensor
inline MaterialTensor::MaterialTensor(const KinematicTensor& other){
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
inline MaterialTensor MaterialTensor::operator-(){
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp[i] = -data_ptr[i];
    }
    return MaterialTensor(tmp.data());
}

//define A += B
inline MaterialTensor& MaterialTensor::operator+= (const MaterialTensor &rhs){
    for (int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = data_ptr[i]+rhs(i);
    }
    return *this;
}

inline MaterialTensor& MaterialTensor::operator+= (const KinematicTensor &rhs){
    if (rhs.DIM == 1){
        data_ptr[XX] += rhs(XX);
    } else if (rhs.DIM ==2){
        data_ptr[XX] += rhs(XX); data_ptr[XY] += rhs(XY);
        data_ptr[YX] += rhs(YX); data_ptr[YY] += rhs(YY);
    } else {
        for (int i = 0; i < TENSOR_MAX_LENGTH; i++) {
            data_ptr[i] = data_ptr[i] + rhs(i);
        }
    }
    return *this;
}

//define A -= B
inline MaterialTensor& MaterialTensor::operator-= (const MaterialTensor &rhs){
    for (int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = data_ptr[i]-rhs(i);
    }
    return *this;
}

inline MaterialTensor& MaterialTensor::operator-= (const KinematicTensor &rhs){
    if (rhs.DIM == 1){
        data_ptr[XX] -= rhs(XX);
    } else if (rhs.DIM ==2){
        data_ptr[XX] -= rhs(XX); data_ptr[XY] -= rhs(XY);
        data_ptr[YX] -= rhs(YX); data_ptr[YY] -= rhs(YY);
    } else {
        for (int i = 0; i < TENSOR_MAX_LENGTH; i++) {
            data_ptr[i] = data_ptr[i] - rhs(i);
        }
    }
    return *this;
}

//define A *= s for integer values
inline MaterialTensor& MaterialTensor::operator*= (const int &rhs){
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = data_ptr[i]*rhs;
    }
    return *this;
}

//define A *= s for double values
inline MaterialTensor& MaterialTensor::operator*= (const double &rhs){
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = data_ptr[i]*rhs;
    }
    return *this;
}

//define A /= s for integer values
inline MaterialTensor& MaterialTensor::operator/= (const int &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = data_ptr[i]*tmp;
    }
    return *this;
}

//define A /= s for double values
inline MaterialTensor& MaterialTensor::operator/= (const double &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = data_ptr[i]*tmp;
    }
    return *this;
}

/*------------------------------------------------------------------------*/
//set A to zero tensor
inline void MaterialTensor::setZero(){
    for (int i=0;i<TENSOR_MAX_LENGTH;i++){
        data_ptr[i] = 0;
    }
    return;
}

//set A to identity tensor
inline void MaterialTensor::setIdentity(){
    data_ptr[XX] = 1; data_ptr[XY] = 0; data_ptr[XZ] = 0;
    data_ptr[YX] = 0; data_ptr[YY] = 1; data_ptr[YZ] = 0;
    data_ptr[ZX] = 0; data_ptr[ZY] = 0; data_ptr[ZZ] = 1;
    return;
}

//return identity tensor
inline MaterialTensor MaterialTensor::Identity() {
    MaterialTensor tmp;
    tmp[XX] = 1; tmp[XY] = 0; tmp[XZ] = 0;
    tmp[YX] = 0; tmp[YY] = 1; tmp[YZ] = 0;
    tmp[ZX] = 0; tmp[ZY] = 0; tmp[ZZ] = 1;
    return tmp;
}

/*------------------------------------------------------------------------*/
//tensor inner products
inline double MaterialTensor::dot(const MaterialTensor &rhs) const{
    double tmp=0;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp += data_ptr[i]*rhs[i];
    }
    return tmp;
}

inline double MaterialTensor::dot(const KinematicTensor &rhs) const{
    double tmp=0;
    for(int i=0;i<rhs.DIM;i++){
        for(int j=0;j<rhs.DIM;j++) {
            tmp += data_ptr[TENSOR_MAX_DIM*i+j] * rhs(i,j);
        }
    }
    return tmp;
}

//trace of MaterialTensor
inline double MaterialTensor::trace() const{
    return data_ptr[XX]+data_ptr[YY]+data_ptr[ZZ];
}

//tensor norm given by sqrt(A:A)
inline double MaterialTensor::norm() const{
    double tmp = 0;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp += data_ptr[i]*data_ptr[i];
    }
    return std::sqrt(tmp);
}

//determinant of MaterialTensor
inline double MaterialTensor::det() const{
    double a = data_ptr[XX]*(data_ptr[YY]*data_ptr[ZZ] - data_ptr[ZY]*data_ptr[YZ]);
    double b = data_ptr[XY]*(data_ptr[YX]*data_ptr[ZZ] - data_ptr[ZX]*data_ptr[YZ]);
    double c = data_ptr[XZ]*(data_ptr[YX]*data_ptr[ZY] - data_ptr[ZX]*data_ptr[YY]);
    return a-b+c;
}

/*------------------------------------------------------------------------*/
//tensor inverse of MaterialTensor (not safe)
inline MaterialTensor MaterialTensor::inverse() const{
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    double detA = det();
    tmp[XX] = 1/detA * (data_ptr[YY]*data_ptr[ZZ] - data_ptr[ZY]*data_ptr[YZ]);
    tmp[XY] = 1/detA * (data_ptr[XZ]*data_ptr[ZY] - data_ptr[ZZ]*data_ptr[XY]);
    tmp[XZ] = 1/detA * (data_ptr[XY]*data_ptr[YZ] - data_ptr[YY]*data_ptr[XZ]);
    tmp[YX] = 1/detA * (data_ptr[YZ]*data_ptr[ZX] - data_ptr[ZZ]*data_ptr[YX]);
    tmp[YY] = 1/detA * (data_ptr[XX]*data_ptr[ZZ] - data_ptr[ZX]*data_ptr[XZ]);
    tmp[YZ] = 1/detA * (data_ptr[XZ]*data_ptr[YX] - data_ptr[YZ]*data_ptr[XX]);
    tmp[ZX] = 1/detA * (data_ptr[YX]*data_ptr[ZY] - data_ptr[ZX]*data_ptr[YY]);
    tmp[ZY] = 1/detA * (data_ptr[XY]*data_ptr[ZX] - data_ptr[ZY]*data_ptr[XX]);
    tmp[ZZ] = 1/detA * (data_ptr[XX]*data_ptr[YY] - data_ptr[YX]*data_ptr[XY]);
    if (!std::isfinite(1/detA)){
        std::cerr << "WARNING: MaterialTensor::inverse() has resulted in non-finite value. So that's bad." << std::endl;
    }
    return MaterialTensor(tmp.data());
}

//tensor deviator (A_0 = A - 1/3 tr(A) I)
inline MaterialTensor MaterialTensor::deviator() const{
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    double trA = trace();
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp[i] = data_ptr[i];
    }
    tmp[XX] -= trA/3.0; tmp[YY] -= trA/3.0; tmp[ZZ] -= trA/3.0;
    return MaterialTensor(tmp.data());
}

//transpose (A^T)
inline MaterialTensor MaterialTensor::transpose() const{
    std::array<double,9> tmp;
    tmp[XX] = data_ptr[XX]; tmp[XY] = data_ptr[YX]; tmp[XZ] = data_ptr[ZX];
    tmp[YX] = data_ptr[XY]; tmp[YY] = data_ptr[YY]; tmp[YZ] = data_ptr[ZY];
    tmp[ZX] = data_ptr[XZ]; tmp[ZY] = data_ptr[YZ]; tmp[ZZ] = data_ptr[ZZ];
    return MaterialTensor(tmp.data());
}

//symmetric and skew decompositions of A
inline MaterialTensor MaterialTensor::sym() const{
    MaterialTensor tmp(this->data());
    return 0.5*(tmp + tmp.transpose());
}
inline MaterialTensor MaterialTensor::skw() const{
    MaterialTensor tmp(this->data());
    return 0.5*(tmp - tmp.transpose());
}

//material tensor exponent
inline MaterialTensor MaterialTensor::exp(int order) const{
    MaterialTensor tmp(this->data());
    MaterialTensor output = MaterialTensor::Identity();

    //exp(X) = I + \sum(1/k! * X^k)
    for (int k=order; k>=1; k--){
        output = MaterialTensor::Identity() + 1.0/k * tmp * output;
    }

    return output;
}


/*------------------------------------------------------------------------*/
//construct KinematicTensor from data pointer and input tensor_type
inline KinematicTensor::KinematicTensor(const double* otherdata, int input){
    assignTensorType(input);
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++){
            data_ptr[TENSOR_MAX_DIM * i + j] = otherdata[TENSOR_MAX_DIM * i + j];
        }
    }
    for (int i=DIM;i<3;i++) {
        for (int j = DIM; j < TENSOR_MAX_DIM; j++) {
            if (otherdata[3*i+j] != 0){
                std::cerr << "WARNING: ignoring non-zero entries in input data to KinematicTensor. (" << otherdata[3*i+j] << ")" << std::endl;
            }
            data_ptr[TENSOR_MAX_DIM * i + j] = 0;
        }
    }
}

//copy KinematicTensor from MaterialTensor
inline KinematicTensor::KinematicTensor(const MaterialTensor& other, int input){
    assignTensorType(input);
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++){
            data_ptr[TENSOR_MAX_DIM * i + j] = other[TENSOR_MAX_DIM * i + j];
        }
    }
    for (int i=DIM;i<3;i++) {
        for (int j = DIM; j < TENSOR_MAX_DIM; j++) {
            if (other[TENSOR_MAX_DIM * i + j] != 0){
                std::cerr << "WARNING: ignoring non-zero entries in input data to KinematicTensor. (" << other[TENSOR_MAX_DIM * i + j] << ")" << std::endl;
            }
            data_ptr[TENSOR_MAX_DIM * i + j] = 0;
        }
    }
}

/*------------------------------------------------------------------------*/
//define -A
inline KinematicTensor KinematicTensor::operator-(){
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp[i] = -data_ptr[i];
    }
    return KinematicTensor(tmp.data(),TENSOR_TYPE);
}

//define A += B
inline KinematicTensor& KinematicTensor::operator+= (const KinematicTensor &rhs){
    assert(TENSOR_TYPE == rhs.TENSOR_TYPE && "Addition failed.");
    if (DIM == 1){
        data_ptr[XX] += rhs(XX);
    } else if (DIM ==2){
        data_ptr[XX] += rhs(XX); data_ptr[XY] += rhs(XY);
        data_ptr[YX] += rhs(YX); data_ptr[YY] += rhs(YY);
    } else {
        for (int i = 0; i < TENSOR_MAX_LENGTH; i++) {
            data_ptr[i] = data_ptr[i] + rhs(i);
        }
    }
    return *this;
}

//define A -= B
inline KinematicTensor& KinematicTensor::operator-= (const KinematicTensor &rhs){
    if (DIM == 1){
        data_ptr[XX] -= rhs(XX);
    } else if (DIM ==2){
        data_ptr[XX] -= rhs(XX); data_ptr[XY] -= rhs(XY);
        data_ptr[YX] -= rhs(YX); data_ptr[YY] -= rhs(YY);
    } else {
        for (int i = 0; i < TENSOR_MAX_LENGTH; i++) {
            data_ptr[i] = data_ptr[i] - rhs(i);
        }
    }
    return *this;
}

//define A *= s for integer values
inline KinematicTensor& KinematicTensor::operator*= (const int &rhs){
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
inline KinematicTensor& KinematicTensor::operator*= (const double &rhs){
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
inline KinematicTensor& KinematicTensor::operator/= (const int &rhs){
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
inline KinematicTensor& KinematicTensor::operator/= (const double &rhs){
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
inline void KinematicTensor::setZero(){
    for (int i = 0; i < TENSOR_MAX_LENGTH; i++) {
        data_ptr[i] = 0;
    }
    return;
}

inline void KinematicTensor::setIdentity(){
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

//return identity tensor
inline KinematicTensor KinematicTensor::Identity(int input) {
    KinematicTensor tmp(input);
    if (tmp.DIM == 1){
        tmp[XX] = 1;
    } else if (tmp.DIM == 2){
        tmp[XX] = 1; tmp[XY] = 0;
        tmp[YX] = 0; tmp[YY] = 1;
    } else {
        tmp[XX] = 1; tmp[XY] = 0; tmp[XZ] = 0;
        tmp[YX] = 0; tmp[YY] = 1; tmp[YZ] = 0;
        tmp[ZX] = 0; tmp[ZY] = 0; tmp[ZZ] = 1;
    }
    return tmp;
}

/*------------------------------------------------------------------------*/
//tensor inner products
inline double KinematicTensor::dot(const MaterialTensor &rhs) const{
    double tmp=0;
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++) {
            tmp += data_ptr[TENSOR_MAX_DIM*i+j] * rhs(i,j);
        }
    }
    return tmp;
}

inline double KinematicTensor::dot(const KinematicTensor &rhs) const{
    double tmp=0;
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++) {
            tmp += data_ptr[TENSOR_MAX_DIM*i+j] * rhs(i,j);
        }
    }
    return tmp;
}

//trace of KinematicTensor
inline double KinematicTensor::trace() const{
    return data_ptr[XX]+data_ptr[YY]+data_ptr[ZZ];
}

//tensor norm given by sqrt(A:A)
inline double KinematicTensor::norm() const{
    double tmp = 0;
    for(int i=0;i<DIM;i++){
        for(int j=0;j<DIM;j++) {
            tmp += data_ptr[TENSOR_MAX_DIM * i + j] * data_ptr[TENSOR_MAX_DIM * i + j];
        }
    }
    return std::sqrt(tmp);
}

//determinant of KinematicTensor
inline double KinematicTensor::det() const{
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
inline KinematicTensor KinematicTensor::inverse() const{
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
        tmp[ZY] = 1 / detA * (data_ptr[XY] * data_ptr[ZX] - data_ptr[ZY] * data_ptr[XX]);
        tmp[ZZ] = 1 / detA * (data_ptr[XX] * data_ptr[YY] - data_ptr[YX] * data_ptr[XY]);
    }
    if (!std::isfinite(1/detA)){
        std::cerr << "WARNING: MaterialTensor::inverse() has resulted in non-finite value. So that's bad." << std::endl;
    }
    return KinematicTensor(tmp.data(),TENSOR_TYPE);
}

//tensor deviator (A_0 = A - 1/3 tr(A) I)
inline MaterialTensor KinematicTensor::deviator() const{
    std::array<double,TENSOR_MAX_LENGTH> tmp;
    double trA = trace();
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){
        tmp[i] = data_ptr[i];
    }
    tmp[XX] -= trA/3.0; tmp[YY] -= trA/3.0; tmp[ZZ] -= trA/3.0;
    return MaterialTensor(tmp.data());
}

//transpose (A^T)
inline KinematicTensor KinematicTensor::transpose() const{
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
inline KinematicTensor KinematicTensor::sym() const{
    KinematicTensor tmp(this->data(),this->TENSOR_TYPE);
    return 0.5*(tmp + tmp.transpose());
}
inline KinematicTensor KinematicTensor::skw() const{
    KinematicTensor tmp(this->data(),this->TENSOR_TYPE);
    return 0.5*(tmp - tmp.transpose());
}

//kinematic tensor exponent
inline KinematicTensor KinematicTensor::exp(int order) const{
    KinematicTensor tmp(this->data(), this->TENSOR_TYPE);
    KinematicTensor output = KinematicTensor::Identity(this->TENSOR_TYPE);

    //exp(X) = I + \sum(1/k! * X^k)
    for (int k=order; k>=1; k--){
        output = KinematicTensor::Identity(this->TENSOR_TYPE) + 1.0/k * tmp * output;
    }

    return output;
}

/*----------------------------------------------------------------------------*/
//Tensor Addition
inline MaterialTensor operator+ (const MaterialTensor& lhs, const MaterialTensor& rhs){
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    return MaterialTensor(tmp.data());
}

inline MaterialTensor operator+ (const MaterialTensor& lhs, const KinematicTensor& rhs){
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    return MaterialTensor(tmp.data());
}

inline MaterialTensor operator+ (const KinematicTensor& lhs, const MaterialTensor& rhs){
    return rhs+lhs;
}

inline KinematicTensor operator+ (const KinematicTensor& lhs, const KinematicTensor& rhs){
    assert(lhs.TENSOR_TYPE == rhs.TENSOR_TYPE && "Addition failed.");
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    return KinematicTensor(tmp.data(),lhs.TENSOR_TYPE);
}

/*----------------------------------------------------------------------------*/
//Tensor Subtraction
inline MaterialTensor operator- (const MaterialTensor& lhs, const MaterialTensor& rhs){
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialTensor(tmp.data());
}

inline MaterialTensor operator- (const MaterialTensor& lhs, const KinematicTensor& rhs){
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialTensor(tmp.data());
}

inline MaterialTensor operator- (const KinematicTensor& lhs, const MaterialTensor& rhs){
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialTensor(tmp.data());
}

inline KinematicTensor operator- (const KinematicTensor& lhs, const KinematicTensor& rhs){
    assert(lhs.TENSOR_TYPE == rhs.TENSOR_TYPE && "Addition failed.");
    std::array<double, MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return KinematicTensor(tmp.data(), lhs.TENSOR_TYPE);
}

/*----------------------------------------------------------------------------*/
//Tensor Multiplication (returning MaterialTensor)
inline MaterialTensor operator* (const double& lhs, const MaterialTensor& rhs){
    MaterialTensor tmp = MaterialTensor(rhs);
    return tmp*=lhs;
}

inline MaterialTensor operator* (const int& lhs, const MaterialTensor& rhs){
    MaterialTensor tmp = MaterialTensor(rhs);
    return tmp*=lhs;
}

inline MaterialTensor operator* (const MaterialTensor& lhs, const double& rhs){
    MaterialTensor tmp = MaterialTensor(lhs);
    return tmp*=rhs;
}

inline MaterialTensor operator* (const MaterialTensor& lhs, const int& rhs){
    MaterialTensor tmp = MaterialTensor(lhs);
    return tmp*=rhs;
}

inline MaterialTensor operator/ (const MaterialTensor& lhs, const double& rhs){
    MaterialTensor tmp = MaterialTensor(lhs);
    return tmp/=rhs;
}

inline MaterialTensor operator/ (const MaterialTensor& lhs, const int& rhs){
    MaterialTensor tmp = MaterialTensor(lhs);
    return tmp/=rhs;
}

inline MaterialTensor operator* (const MaterialTensor& lhs, const MaterialTensor& rhs){
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

inline MaterialTensor operator* (const MaterialTensor& lhs, const KinematicTensor& rhs){
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

inline MaterialTensor operator* (const KinematicTensor& lhs, const MaterialTensor& rhs){
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
inline KinematicTensor operator* (const double& lhs, const KinematicTensor& rhs){
    KinematicTensor tmp = KinematicTensor(rhs);
    return tmp*=lhs;
}

inline KinematicTensor operator* (const int& lhs, const KinematicTensor& rhs){
    KinematicTensor tmp = KinematicTensor(rhs);
    return tmp*=lhs;
}

inline KinematicTensor operator* (const KinematicTensor& lhs, const double& rhs){
    return rhs*lhs;
}

inline KinematicTensor operator* (const KinematicTensor& lhs, const int& rhs){
    return rhs*lhs;
}

inline KinematicTensor operator/ (const KinematicTensor& lhs, const double& rhs){
    KinematicTensor tmp = KinematicTensor(lhs);
    return tmp/=rhs;
}

inline KinematicTensor operator/ (const KinematicTensor& lhs, const int& rhs){
    KinematicTensor tmp = KinematicTensor(lhs);
    return tmp/=rhs;
}

inline KinematicTensor operator* (const KinematicTensor& lhs, const KinematicTensor& rhs){
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
inline MaterialTensor::Map& MaterialTensor::Map::operator= (const MaterialTensor &other){
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){data_ptr[i]=other(i);};
    return *this;
}

inline MaterialTensor::Map& MaterialTensor::Map::operator= (const KinematicTensor &other){
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){data_ptr[i]=other(i);};
    return *this;
}

inline KinematicTensor::Map& KinematicTensor::Map::operator= (const KinematicTensor &other){
    for(int i=0;i<TENSOR_MAX_LENGTH;i++){data_ptr[i]=other(i);};
    return *this;
}

#endif //MPM_V2_TENSOR_HPP

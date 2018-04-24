//
// Created by aaron on 4/22/18.
// algebra.hpp
//

#ifndef MPM_V2_TENSOR_HPP
#define MPM_V2_TENSOR_HPP

#include <stdlib.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <iostream>

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

    MPMTensor& operator= (const MPMTensor& other){
        for (int i=0;i<TENSOR_MAX_LENGTH;i++){
            buffer[i] = other.data_ptr[i];
        }
        data_ptr = buffer.data(); //make sure not to copy stale pointer
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
    double* data(){
        return data_ptr;
    }

    /*------------------------------------------------------------------------*/
    //return size of tensor
    double size(){
        return TENSOR_MAX_LENGTH;
    }

    /*------------------------------------------------------------------------*/
    //return rows and columns
    double rows(){
        return TENSOR_MAX_DIM;
    }

    double cols(){
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
    MaterialTensor(){};

    //copy constructor
    MaterialTensor(const MaterialTensor& other) : MPMTensor(other){}

    //construct tensor from data pointer
    MaterialTensor(double* otherdata);

    //copy from KinematicTensor
    MaterialTensor(KinematicTensor&);

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
    operator EIGEN_MAP_OF_MATERIAL_TENSOR () {return EIGEN_MAP_OF_MATERIAL_TENSOR(data_ptr);}
    MaterialTensor operator-();
    MaterialTensor& operator*= (const int &rhs);
    MaterialTensor& operator*= (const double &rhs);
    MaterialTensor& operator/= (const int &rhs);
    MaterialTensor& operator/= (const double &rhs);

    /*------------------------------------------------------------------------*/
    //set to standard tensors
    void setZero();
    void setIdentity();

    /*------------------------------------------------------------------------*/
    //tensor contractions
    double dot(const MaterialTensor &rhs);
    double dot(const KinematicTensor &rhs);

    //tensor trace
    double trace();

    //tensor frobenius norm
    double norm();

    //tensor determinant
    double det();

    /*------------------------------------------------------------------------*/
    //inverse
    MaterialTensor inverse();

    //deviator
    MaterialTensor deviator();

    //transpose
    MaterialTensor transpose();

    //symetric and skew decomposition
    MaterialTensor sym();
    MaterialTensor skw();
};

/*----------------------------------------------------------------------------*/
class MaterialTensor::Map : public MaterialTensor{
public:
    //point to other data (not safe)
    Map(double* dataIN){ data_ptr = dataIN; }

    //assignment operator
    Map& operator= (const MaterialTensor &other);
    Map& operator= (const KinematicTensor &other);
};

/*----------------------------------------------------------------------------*/
class KinematicTensor : public MPMTensor{
public:
    //KinematicTensor types
    static const int TENSOR_1D = 1, TENSOR_2D = 2, TENSOR_3D = 3, TENSOR_AXISYM = 4;

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
        } else {
            std::cerr << "KinematicTensor doesn't have defined type for input " << TENSOR_TYPE << "." << std::endl;
        }
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
    KinematicTensor(double* otherdata, int input);

    //copy from KinematicTensor
    //KinematicTensor(KinematicTensor&);
    //KinematicTensor(const KinematicTensor&, int input);

    //copy from MaterialTensor
    KinematicTensor(MaterialTensor&, int input);

    //forward declare Map class
    class Map;

    //construct from Eigen::Matrix
    template <typename OtherDerived>
    KinematicTensor(Eigen::MatrixBase<OtherDerived> &other, int input){
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
    operator EIGEN_MAP_OF_KINEMATIC_TENSOR () {return EIGEN_MAP_OF_KINEMATIC_TENSOR(data_ptr,DIM,DIM);}
    KinematicTensor operator-();
    KinematicTensor& operator*= (const int &rhs);
    KinematicTensor& operator*= (const double &rhs);
    KinematicTensor& operator/= (const int &rhs);
    KinematicTensor& operator/= (const double &rhs);

    /*------------------------------------------------------------------------*/
    //set to standard tensors
    void setZero();
    void setIdentity();

    /*------------------------------------------------------------------------*/
    //tensor contractions
    double dot(const MaterialTensor &rhs);
    double dot(const KinematicTensor &rhs);

    //tensor trace
    double trace();

    //tensor frobenius norm
    double norm();

    //tensor determinant
    double det();

    /*------------------------------------------------------------------------*/
    //inverse
    KinematicTensor inverse();

    //deviator
    MaterialTensor deviator();

    //transpose
    KinematicTensor transpose();

    //symetric and skew decomposition
    KinematicTensor sym();
    KinematicTensor skw();
};

/*----------------------------------------------------------------------------*/
class KinematicTensor::Map : public KinematicTensor{
public:
    //point to other data (not safe)
    Map(double* dataIN, int input){
        assignTensorType(input);
        data_ptr = dataIN;
    }

    //assignment operator
    Map& operator= (const KinematicTensor &other);
};

/*----------------------------------------------------------------------------*/
//operations which return a material tensor
MaterialTensor operator+ (const MaterialTensor&, const MaterialTensor&);
MaterialTensor operator+ (const MaterialTensor&, const KinematicTensor&);
MaterialTensor operator+ (const KinematicTensor&, const MaterialTensor&);
MaterialTensor operator- (const MaterialTensor&, const MaterialTensor&);
MaterialTensor operator- (const MaterialTensor&, const KinematicTensor&);
MaterialTensor operator- (const KinematicTensor&, const MaterialTensor&);
MaterialTensor operator* (const double&, const MaterialTensor&);
MaterialTensor operator* (const int&, const MaterialTensor&);
MaterialTensor operator* (const MaterialTensor&, const double&);
MaterialTensor operator* (const MaterialTensor&, const int&);
MaterialTensor operator* (const MaterialTensor&, const MaterialTensor&);
MaterialTensor operator* (const MaterialTensor&, const KinematicTensor&);
MaterialTensor operator* (const KinematicTensor&, const MaterialTensor&);
MaterialTensor operator/ (const MaterialTensor&, const double&);
MaterialTensor operator/ (const MaterialTensor&, const int&);

/*----------------------------------------------------------------------------*/
//operations which return a material tensor
KinematicTensor operator+ (const KinematicTensor&, const KinematicTensor&);
KinematicTensor operator- (const KinematicTensor&, const KinematicTensor&);
KinematicTensor operator* (const double&, const KinematicTensor&);
KinematicTensor operator* (const int&, const KinematicTensor&);
KinematicTensor operator* (const KinematicTensor&, const double&);
KinematicTensor operator* (const KinematicTensor&, const int&);
KinematicTensor operator* (const KinematicTensor&, const KinematicTensor&);
KinematicTensor operator/ (const KinematicTensor&, const double&);
KinematicTensor operator/ (const KinematicTensor&, const int&);

#endif //MPM_V2_TENSOR_HPP

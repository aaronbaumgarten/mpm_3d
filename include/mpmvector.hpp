//
// Created by aaron on 4/23/18.
// mpmvector.hpp
//

#ifndef MPM_V3_MPMVECTOR_HPP
#define MPM_V3_MPMVECTOR_HPP

#include "mpmtensor.hpp"

/*----------------------------------------------------------------------------*/
//Eigen::Map object which interfaces with KinematicVector object.

typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> EIGEN_MAP_OF_KINEMATIC_VECTOR;

/*----------------------------------------------------------------------------*/
//Eigen::Map object which interfaces with MaterialVector object.

typedef Eigen::Map<Eigen::Matrix<double,3,1>> EIGEN_MAP_OF_MATERIAL_VECTOR;

/*----------------------------------------------------------------------------*/
//base class for KinematicVector and MaterialVector.
class MPMVector {
public:
    MPMVector(){
        buffer = { {0} };
        data_ptr = buffer.data();
    }

    MPMVector(const MPMVector& other){
        for (int i=0;i<VECTOR_MAX_DIM;i++){
            buffer[i] = other.data_ptr[i];
        }
        data_ptr = buffer.data(); //make sure not to copy stale pointer
    }

    MPMVector& operator= (const MPMVector& other){
        for (int i=0;i<VECTOR_MAX_DIM;i++){
            buffer[i] = other.data_ptr[i];
        }
        data_ptr = buffer.data(); //make sure not to copy stale pointer
    }

    //set standard vector length
    static const int VECTOR_MAX_DIM = 3;

    //set standard vector access variables
    static const int X = 0,  Y = 1,  Z = 2;

    /*------------------------------------------------------------------------*/
    //functions for accessing vector data
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
    //return size of vector
    double size(){
        return VECTOR_MAX_DIM;
    }

    /*------------------------------------------------------------------------*/
    //return rows and columns
    double rows(){
        return VECTOR_MAX_DIM;
    }

    double cols(){
        return 1;
    }

    /*------------------------------------------------------------------------*/
private:
    //standard buffer for data
    std::array<double, VECTOR_MAX_DIM> buffer;

    /*------------------------------------------------------------------------*/
protected:
    //pointer to data (by default points to standard buffer)
    double* data_ptr;

};

/*----------------------------------------------------------------------------*/
class KinematicVector;

/*----------------------------------------------------------------------------*/
class MaterialVector : public MPMVector{
public:
    //default constructor
    MaterialVector(){};

    //copy constructor
    MaterialVector(const MaterialVector& other) : MPMVector(other){}

    //construct vector from data pointer
    MaterialVector(double* otherdata);

    //copy from KinematicVector
    MaterialVector(KinematicVector&);

    //forward declare Map class
    class Map;

    //construct from Eigen::Matrix
    template <typename OtherDerived>
    MaterialVector(Eigen::MatrixBase<OtherDerived> &other){
        assert(other.rows() == VECTOR_MAX_DIM && other.cols() == 1);
        for(int i=0;i<VECTOR_MAX_DIM;i++){
            data_ptr[i] = other(i);
        }
    }

    /*------------------------------------------------------------------------*/
    //define operators
    operator EIGEN_MAP_OF_MATERIAL_VECTOR () {return EIGEN_MAP_OF_MATERIAL_VECTOR(data_ptr);}
    MaterialVector operator-();
    MaterialVector& operator*= (const int &rhs);
    MaterialVector& operator*= (const double &rhs);
    MaterialVector& operator/= (const int &rhs);
    MaterialVector& operator/= (const double &rhs);

    /*------------------------------------------------------------------------*/
    //set to standard vector
    void setZero();
    void setOnes();

    /*------------------------------------------------------------------------*/
    //vector contractions
    double dot(const MaterialVector &rhs);
    double dot(const KinematicVector &rhs);

    //vector cross product
    MaterialVector cross(const MaterialVector &rhs);
    MaterialVector cross(const KinematicVector &rhs);

    //vector tensor product
    MaterialTensor tensor(const MaterialVector& rhs);
    MaterialTensor tensor(const KinematicVector& rhs);

    //vector 2-norm
    double norm();
};


/*----------------------------------------------------------------------------*/
class MaterialVector::Map : public MaterialVector{
public:
    //point to other data (not safe)
    Map(double* dataIN){ data_ptr = dataIN; }

    //assignment operator
    Map& operator= (const MaterialVector &other);
    Map& operator= (const KinematicVector &other);
};


/*----------------------------------------------------------------------------*/
class KinematicVector : public MPMVector{
public:
    //KinematicVector types derived from KinematicMatrix types
    static const int VECTOR_1D = KinematicTensor::TENSOR_1D;
    static const int VECTOR_2D = KinematicTensor::TENSOR_2D;
    static const int VECTOR_3D = KinematicTensor::TENSOR_3D;
    static const int VECTOR_AXISYM = KinematicTensor::TENSOR_AXISYM;

    //default dimension and type
    int DIM = VECTOR_MAX_DIM;
    int VECTOR_TYPE = VECTOR_3D;

    void assignVectorType(int input){
        VECTOR_TYPE = input;
        if (VECTOR_TYPE == VECTOR_1D){
            DIM = 1;
        } else if (VECTOR_TYPE == VECTOR_2D){
            DIM = 2;
        } else if (VECTOR_TYPE == VECTOR_3D){
            DIM = 3;
        } else {
            std::cerr << "KinematicVector doesn't have defined type for input " << VECTOR_TYPE << "." << std::endl;
        }
    }

    //default constructor
    KinematicVector(){};
    KinematicVector(int input){
        assignVectorType(input);
    }

    KinematicVector(const KinematicVector& other) : MPMVector(other){
        DIM = other.DIM;
        VECTOR_TYPE = other.VECTOR_TYPE;
    }

    //construct vector from data pointer
    KinematicVector(double* otherdata, int input);

    //copy from MaterialVector
    KinematicVector(MaterialVector&, int input);

    //forward declare Map class
    class Map;

    //construct from Eigen::Matrix
    template <typename OtherDerived>
    KinematicVector(Eigen::MatrixBase<OtherDerived> &other, int input){
        assignVectorType(input);
        assert(other.rows() == DIM && other.cols() == 1);
        for(int i=0;i<DIM;i++){
            data_ptr[i] = other(i);
        }
        for (int i=DIM;i<3;i++) {
            data_ptr[i] = 0;
        }
    }

    /*------------------------------------------------------------------------*/
    //define operators
    operator EIGEN_MAP_OF_KINEMATIC_VECTOR () {return EIGEN_MAP_OF_KINEMATIC_VECTOR(data_ptr,DIM);}
    KinematicVector operator-();
    KinematicVector& operator*= (const int &rhs);
    KinematicVector& operator*= (const double &rhs);
    KinematicVector& operator/= (const int &rhs);
    KinematicVector& operator/= (const double &rhs);

    /*------------------------------------------------------------------------*/
    //set to standard vectors
    void setZero();
    void setOnes();

    /*------------------------------------------------------------------------*/
    //vector contractions
    double dot(const MaterialVector &rhs);
    double dot(const KinematicVector &rhs);

    //vector cross product
    MaterialVector cross(const MaterialVector &rhs);
    MaterialVector cross(const KinematicVector &rhs);

    //vector tensor product
    MaterialTensor tensor(const MaterialVector& rhs);
    KinematicTensor tensor(const KinematicVector& rhs);

    //vector 2-norm
    double norm();
};


/*----------------------------------------------------------------------------*/
class KinematicVector::Map : public KinematicVector{
public:
    //point to other data (not safe)
    Map(double* dataIN, int input){
        assignVectorType(input);
        data_ptr = dataIN;
    }

    //assignment operator
    Map& operator= (const KinematicVector &other);
};

/*----------------------------------------------------------------------------*/
//operations which return a material vector
MaterialVector operator+ (const MaterialVector&, const MaterialVector&);
MaterialVector operator+ (const MaterialVector&, const KinematicVector&);
MaterialVector operator+ (const KinematicVector&, const MaterialVector&);
MaterialVector operator- (const MaterialVector&, const MaterialVector&);
MaterialVector operator- (const MaterialVector&, const KinematicVector&);
MaterialVector operator- (const KinematicVector&, const MaterialVector&);
MaterialVector operator* (const double&, const MaterialVector&);
MaterialVector operator* (const int&, const MaterialVector&);
MaterialVector operator* (const MaterialVector&, const double&);
MaterialVector operator* (const MaterialVector&, const int&);
MaterialVector operator* (const MaterialTensor&, const MaterialVector&);
MaterialVector operator* (const MaterialTensor&, const KinematicVector&);
MaterialVector operator* (const KinematicTensor&, const MaterialVector&);
MaterialVector operator/ (const MaterialVector&, const double&);
MaterialVector operator/ (const MaterialVector&, const int&);

/*----------------------------------------------------------------------------*/
//operations which return a kinematic vector
KinematicVector operator+ (const KinematicVector&, const KinematicVector&);
KinematicVector operator- (const KinematicVector&, const KinematicVector&);
KinematicVector operator* (const double&, const KinematicVector&);
KinematicVector operator* (const int&, const KinematicVector&);
KinematicVector operator* (const KinematicVector&, const double&);
KinematicVector operator* (const KinematicVector&, const int&);
KinematicVector operator* (const KinematicTensor&, const KinematicVector&);
KinematicVector operator/ (const KinematicVector&, const double&);
KinematicVector operator/ (const KinematicVector&, const int&);

#endif //MPM_V3_MPMVECTOR_HPP

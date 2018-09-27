//
// Created by aaron on 4/23/18.
// mpmvector.hpp
//

#ifndef MPM_V3_MPMVECTOR_HPP
#define MPM_V3_MPMVECTOR_HPP

#include "mpm_tensor.hpp"

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

    MPMVector(const MPMVector&& other){
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

    MPMVector& operator= (const MPMVector&& other){
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
    double* data() const {
        return data_ptr;
    }

    /*------------------------------------------------------------------------*/
    //return size of vector
    virtual double size() const {
        return VECTOR_MAX_DIM;
    }

    /*------------------------------------------------------------------------*/
    //return rows and columns
    virtual double rows() const {
        return VECTOR_MAX_DIM;
    }

    double cols() const {
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
    MaterialVector(const KinematicVector&);

    //forward declare Map class
    class Map;

    //construct from Eigen::Matrix
    template <typename OtherDerived>
    MaterialVector(const Eigen::MatrixBase<OtherDerived> &other){
        assert(other.rows() == VECTOR_MAX_DIM && other.cols() == 1);
        for(int i=0;i<VECTOR_MAX_DIM;i++){
            data_ptr[i] = other(i);
        }
    }

    /*------------------------------------------------------------------------*/
    //define operators
    operator EIGEN_MAP_OF_MATERIAL_VECTOR () const {return EIGEN_MAP_OF_MATERIAL_VECTOR(data_ptr);}
    MaterialVector operator-();
    inline MaterialVector& operator+= (const MaterialVector &rhs);
    inline MaterialVector& operator-= (const MaterialVector &rhs);
    inline MaterialVector& operator+= (const KinematicVector &rhs);
    inline MaterialVector& operator-= (const KinematicVector &rhs);
    inline MaterialVector& operator*= (const int &rhs);
    inline MaterialVector& operator*= (const double &rhs);
    inline MaterialVector& operator/= (const int &rhs);
    inline MaterialVector& operator/= (const double &rhs);

    /*------------------------------------------------------------------------*/
    //set to standard vector
    void setZero();
    void setOnes();

    /*------------------------------------------------------------------------*/
    //vector contractions
    double dot(const MaterialVector &rhs) const;
    double dot(const KinematicVector &rhs) const;

    //vector cross product
    MaterialVector cross(const MaterialVector &rhs) const;
    MaterialVector cross(const KinematicVector &rhs) const;

    //vector tensor product
    MaterialTensor tensor(const MaterialVector& rhs) const;
    MaterialTensor tensor(const KinematicVector& rhs) const;

    //vector 2-norm
    double norm() const;
};


/*----------------------------------------------------------------------------*/
class MaterialVector::Map : public MaterialVector{
public:
    //point to other data (not safe)
    Map(const double* dataIN){ data_ptr = (double*)(dataIN); }

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
    static const int VECTOR_2D_OOP = KinematicTensor::TENSOR_2D_OOP;

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
        } else if (VECTOR_TYPE == VECTOR_2D_OOP){
            DIM = 3;
        } else {
            std::cerr << "KinematicVector doesn't have defined type for input " << VECTOR_TYPE << "." << std::endl;
        }
    }

    /*------------------------------------------------------------------------*/
    //return size of vector
    double size() const {
        return DIM;
    }

    /*------------------------------------------------------------------------*/
    //return rows and columns
    double rows() const {
        return DIM;
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
    KinematicVector(const MaterialVector&, int input);

    //forward declare Map class
    class Map;

    //construct from Eigen::Matrix
    template <typename OtherDerived>
    KinematicVector(const Eigen::MatrixBase<OtherDerived> &other, int input){
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
    operator EIGEN_MAP_OF_KINEMATIC_VECTOR () const {return EIGEN_MAP_OF_KINEMATIC_VECTOR(data_ptr,DIM);}
    inline KinematicVector operator-();
    inline KinematicVector& operator+= (const KinematicVector &rhs);
    inline KinematicVector& operator-= (const KinematicVector &rhs);
    inline KinematicVector& operator*= (const int &rhs);
    inline KinematicVector& operator*= (const double &rhs);
    inline KinematicVector& operator/= (const int &rhs);
    inline KinematicVector& operator/= (const double &rhs);

    /*------------------------------------------------------------------------*/
    //set to standard vectors
    void setZero();
    void setOnes();

    /*------------------------------------------------------------------------*/
    //vector contractions
    double dot(const MaterialVector &rhs) const;
    double dot(const KinematicVector &rhs) const;

    //vector cross product
    MaterialVector cross(const MaterialVector &rhs) const;
    MaterialVector cross(const KinematicVector &rhs) const;

    //vector tensor product
    MaterialTensor tensor(const MaterialVector& rhs) const;
    KinematicTensor tensor(const KinematicVector& rhs) const;

    //vector 2-norm
    double norm() const;
};


/*----------------------------------------------------------------------------*/
class KinematicVector::Map : public KinematicVector{
public:
    //point to other data (not safe)
    Map(const double* dataIN, const int input){
        assignVectorType(input);
        data_ptr = (double*)dataIN;
    }

    //assignment operator
    Map& operator= (const KinematicVector &other);
};

/*----------------------------------------------------------------------------*/
//operations which return a material vector
inline MaterialVector operator+ (const MaterialVector&, const MaterialVector&);
inline MaterialVector operator+ (const MaterialVector&, const KinematicVector&);
inline MaterialVector operator+ (const KinematicVector&, const MaterialVector&);
inline MaterialVector operator- (const MaterialVector&, const MaterialVector&);
inline MaterialVector operator- (const MaterialVector&, const KinematicVector&);
inline MaterialVector operator- (const KinematicVector&, const MaterialVector&);
inline MaterialVector operator* (const double&, const MaterialVector&);
inline MaterialVector operator* (const int&, const MaterialVector&);
inline MaterialVector operator* (const MaterialVector&, const double&);
inline MaterialVector operator* (const MaterialVector&, const int&);
inline MaterialVector operator* (const MaterialTensor&, const MaterialVector&);
inline MaterialVector operator* (const MaterialTensor&, const KinematicVector&);
inline MaterialVector operator* (const KinematicTensor&, const MaterialVector&);
inline MaterialVector operator/ (const MaterialVector&, const double&);
inline MaterialVector operator/ (const MaterialVector&, const int&);

/*----------------------------------------------------------------------------*/
//operations which return a kinematic vector
inline KinematicVector operator+ (const KinematicVector&, const KinematicVector&);
inline KinematicVector operator- (const KinematicVector&, const KinematicVector&);
inline KinematicVector operator* (const double&, const KinematicVector&);
inline KinematicVector operator* (const int&, const KinematicVector&);
inline KinematicVector operator* (const KinematicVector&, const double&);
inline KinematicVector operator* (const KinematicVector&, const int&);
inline KinematicVector operator* (const KinematicTensor&, const KinematicVector&);
inline KinematicVector operator/ (const KinematicVector&, const double&);
inline KinematicVector operator/ (const KinematicVector&, const int&);





/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/




/*----------------------------------------------------------------------------*/
//construct MaterialVector with pointer to data (not safe, but who cares right?)
inline MaterialVector::MaterialVector(double* otherdata){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] = otherdata[i];
    }
}

//construct MaterialVector from KinematicVector
inline MaterialVector::MaterialVector(const KinematicVector& other){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] = other[i];
    }
}

/*----------------------------------------------------------------------------*/
//define -v
inline MaterialVector MaterialVector::operator-(){
    std::array<double, VECTOR_MAX_DIM> tmp;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        tmp[i] = -data_ptr[i];
    }
    return MaterialVector(tmp.data());
}

//define v += s and v -= s
MaterialVector& MaterialVector::operator+= (const MaterialVector &rhs){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] += rhs(i);
    }
    return *this;
}

MaterialVector& MaterialVector::operator-= (const MaterialVector &rhs){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] -= rhs(i);
    }
    return *this;
}

MaterialVector& MaterialVector::operator+= (const KinematicVector &rhs){
    for(int i=0;i<rhs.DIM;i++){
        data_ptr[i] += rhs(i);
    }
    return *this;
}

MaterialVector& MaterialVector::operator-= (const KinematicVector &rhs){
    for(int i=0;i<rhs.DIM;i++){
        data_ptr[i] += rhs(i);
    }
    return *this;
}

//define v *= s for integer scalar value
inline MaterialVector& MaterialVector::operator*= (const int &rhs){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] *= rhs;
    }
    return *this;
}

//define v *= s for double scalar
inline MaterialVector& MaterialVector::operator*= (const double &rhs){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] *= rhs;
    }
    return *this;
}

//define v /= s for integer scalar value
inline MaterialVector& MaterialVector::operator/= (const int &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] *= tmp;
    }
    return *this;
}

//define v /= s for double scalar
inline MaterialVector& MaterialVector::operator/= (const double &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] *= tmp;
    }
    return *this;
}

/*----------------------------------------------------------------------------*/
//set MaterialVector to standard vectors
inline void MaterialVector::setZero(){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] = 0;
    }
    return;
}
inline void MaterialVector::setOnes(){
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] = 1;
    }
    return;
}

/*----------------------------------------------------------------------------*/
//vector inner product
inline double MaterialVector::dot(const MaterialVector &rhs) const{
    double tmp = 0;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        tmp += data_ptr[i]*rhs[i];
    }
    return tmp;
}

inline double MaterialVector::dot(const KinematicVector &rhs) const{
    double tmp = 0;
    for(int i=0;i<rhs.DIM;i++){
        tmp += data_ptr[i]*rhs[i];
    }
    return tmp;
}

//vector cross product
inline MaterialVector MaterialVector::cross(const MaterialVector &rhs) const{
    std::array<double, VECTOR_MAX_DIM> tmp;
    tmp[X] = data_ptr[Y]*rhs[Z] - data_ptr[Z]*rhs[Y];
    tmp[Y] = data_ptr[Z]*rhs[X] - data_ptr[X]*rhs[Z];
    tmp[Z] = data_ptr[X]*rhs[Y] - data_ptr[Y]*rhs[X];
    return MaterialVector(tmp.data());
}

inline MaterialVector MaterialVector::cross(const KinematicVector &rhs) const{
    std::array<double, VECTOR_MAX_DIM> tmp;
    tmp[X] = data_ptr[Y]*rhs[Z] - data_ptr[Z]*rhs[Y];
    tmp[Y] = data_ptr[Z]*rhs[X] - data_ptr[X]*rhs[Z];
    tmp[Z] = data_ptr[X]*rhs[Y] - data_ptr[Y]*rhs[X];
    return MaterialVector(tmp.data());
}

//vector tensor product
inline MaterialTensor MaterialVector::tensor(const MaterialVector& rhs) const{
    std::array<double,MPMTensor::TENSOR_MAX_LENGTH> tmp;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        for(int j=0;j<VECTOR_MAX_DIM;j++){
            tmp[VECTOR_MAX_DIM * i + j] = data_ptr[i]*rhs[j];
        }
    }
    return MaterialTensor(tmp.data());
}
inline MaterialTensor MaterialVector::tensor(const KinematicVector& rhs) const{
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
inline double MaterialVector::norm() const{
    double tmp = 0;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        tmp += data_ptr[i]*data_ptr[i];
    }
    return std::sqrt(tmp);
}


/*----------------------------------------------------------------------------*/
//construct KinematicVector from pointer to data and input type
inline KinematicVector::KinematicVector(double* otherdata, int input){
    assignVectorType(input);
    for(int i=0;i<DIM;i++){
        data_ptr[i] = otherdata[i];
    }
    for(int i=DIM;i<VECTOR_MAX_DIM;i++){
        data_ptr[i] = 0;
    }
}

//construct KinematicVector from MaterialVector
inline KinematicVector::KinematicVector(const MaterialVector& other, int input){
    assignVectorType(input);
    for(int i=0;i<DIM;i++){
        data_ptr[i] = other[i];
    }
    for(int i=DIM;i<VECTOR_MAX_DIM;i++){
        if (other[i] != 0){
            std::cerr << "WARNING: ignoring non-zero entries in input data to KinematicVector. (" << other[i] << ")" << std::endl;
        }
        data_ptr[i] = 0;
    }
}

/*----------------------------------------------------------------------------*/
//define -v
inline KinematicVector KinematicVector::operator-(){
    std::array<double, VECTOR_MAX_DIM> tmp;
    for(int i=0;i<VECTOR_MAX_DIM;i++){
        tmp[i] = -data_ptr[i];
    }
    return KinematicVector(tmp.data(),VECTOR_TYPE);
}

//define v += s and v -= s
KinematicVector& KinematicVector::operator+= (const KinematicVector &rhs){
    assert(VECTOR_TYPE == rhs.VECTOR_TYPE && "Addition failed.");
    for(int i=0;i<DIM;i++){
        data_ptr[i] += rhs(i);
    }
    return *this;
}

KinematicVector& KinematicVector::operator-= (const KinematicVector &rhs){
    for(int i=0;i<DIM;i++){
        data_ptr[i] -= rhs(i);
    }
    return *this;
}

//define v *= s for integer scalar value
inline KinematicVector& KinematicVector::operator*= (const int &rhs){
    for(int i=0;i<DIM;i++){
        data_ptr[i] *= rhs;
    }
    return *this;
}

//define v *= s for double scalar
inline KinematicVector& KinematicVector::operator*= (const double &rhs){
    for(int i=0;i<DIM;i++){
        data_ptr[i] *= rhs;
    }
    return *this;
}

//define v /= s for integer scalar value
inline KinematicVector& KinematicVector::operator/= (const int &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<DIM;i++){
        data_ptr[i] *= tmp;
    }
    return *this;
}

//define v /= s for double scalar
inline KinematicVector& KinematicVector::operator/= (const double &rhs){
    double tmp = 1.0/rhs;
    for(int i=0;i<DIM;i++){
        data_ptr[i] *= tmp;
    }
    return *this;
}

/*----------------------------------------------------------------------------*/
//set KinematicVector to standard vectors
inline void KinematicVector::setZero(){
    for(int i=0;i<DIM;i++){
        data_ptr[i] = 0;
    }
    return;
}

inline void KinematicVector::setOnes(){
    for(int i=0;i<DIM;i++){
        data_ptr[i] = 1;
    }
    return;
}

/*----------------------------------------------------------------------------*/
//vector contraction
inline double KinematicVector::dot(const MaterialVector &rhs) const{
    double tmp = 0;
    for(int i=0;i<DIM;i++){
        tmp += data_ptr[i]*rhs[i];
    }
    return tmp;
}

inline double KinematicVector::dot(const KinematicVector &rhs) const{
    double tmp = 0;
    for(int i=0;i<DIM;i++){
        tmp += data_ptr[i]*rhs[i];
    }
    return tmp;
}

//vector cross product (always return MaterialVector)
inline MaterialVector KinematicVector::cross(const MaterialVector &rhs) const{
    std::array<double, VECTOR_MAX_DIM> tmp;
    tmp[X] = data_ptr[Y]*rhs[Z] - data_ptr[Z]*rhs[Y];
    tmp[Y] = data_ptr[Z]*rhs[X] - data_ptr[X]*rhs[Z];
    tmp[Z] = data_ptr[X]*rhs[Y] - data_ptr[Y]*rhs[X];
    return MaterialVector(tmp.data());
}

inline MaterialVector KinematicVector::cross(const KinematicVector &rhs) const{
    std::array<double, VECTOR_MAX_DIM> tmp;
    tmp[X] = data_ptr[Y]*rhs[Z] - data_ptr[Z]*rhs[Y];
    tmp[Y] = data_ptr[Z]*rhs[X] - data_ptr[X]*rhs[Z];
    tmp[Z] = data_ptr[X]*rhs[Y] - data_ptr[Y]*rhs[X];
    return MaterialVector(tmp.data());
}

//vector tensor product
inline MaterialTensor KinematicVector::tensor(const MaterialVector& rhs) const{
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

inline KinematicTensor KinematicVector::tensor(const KinematicVector& rhs) const{
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
inline double KinematicVector::norm() const{
    double tmp = 0;
    for(int i=0;i<DIM;i++){
        tmp += data_ptr[i]*data_ptr[i];
    }
    return std::sqrt(tmp);
}


/*----------------------------------------------------------------------------*/
//Vector Addition
inline MaterialVector operator+ (const MaterialVector& rhs, const MaterialVector& lhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    return MaterialVector(tmp.data());
}

inline MaterialVector operator+ (const MaterialVector& rhs, const KinematicVector& lhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = lhs[i]+rhs[i];
    }
    return MaterialVector(tmp.data());
}

inline MaterialVector operator+ (const KinematicVector& lhs, const MaterialVector& rhs){
    return (rhs+lhs);
}

inline KinematicVector operator+ (const KinematicVector& lhs, const KinematicVector& rhs){
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
inline MaterialVector operator- (const MaterialVector& lhs, const MaterialVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialVector(tmp.data());
}

inline MaterialVector operator- (const MaterialVector& lhs, const KinematicVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialVector(tmp.data());
}

inline MaterialVector operator- (const KinematicVector& lhs, const MaterialVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        tmp[i] = lhs[i]-rhs[i];
    }
    return MaterialVector(tmp.data());
}

inline KinematicVector operator- (const KinematicVector& lhs, const KinematicVector& rhs){
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
inline MaterialVector operator* (const double& lhs, const MaterialVector& rhs){
    MaterialVector tmp = MaterialVector(rhs);
    return tmp*=lhs;
}

inline MaterialVector operator* (const int& lhs, const MaterialVector& rhs){
    MaterialVector tmp = MaterialVector(rhs);
    return tmp*=lhs;
}

inline MaterialVector operator* (const MaterialVector& lhs, const double& rhs){
    MaterialVector tmp = MaterialVector(lhs);
    return tmp*=rhs;
}

inline MaterialVector operator* (const MaterialVector& lhs, const int& rhs){
    MaterialVector tmp = MaterialVector(lhs);
    return tmp*=rhs;
}

inline MaterialVector operator* (const MaterialTensor& lhs, const MaterialVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp;
    tmp[MPMVector::X] = lhs[MPMTensor::XX]*rhs[MPMVector::X] + lhs[MPMTensor::XY]*rhs[MPMVector::Y] + lhs[MPMTensor::XZ]*rhs[MPMVector::Z];
    tmp[MPMVector::Y] = lhs[MPMTensor::YX]*rhs[MPMVector::X] + lhs[MPMTensor::YY]*rhs[MPMVector::Y] + lhs[MPMTensor::YZ]*rhs[MPMVector::Z];
    tmp[MPMVector::Z] = lhs[MPMTensor::ZX]*rhs[MPMVector::X] + lhs[MPMTensor::ZY]*rhs[MPMVector::Y] + lhs[MPMTensor::ZZ]*rhs[MPMVector::Z];
    return MaterialVector(tmp.data());
}

inline MaterialVector operator* (const MaterialTensor& lhs, const KinematicVector& rhs){
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp = { {0} };
    for (int i=0;i<MPMVector::VECTOR_MAX_DIM;i++){
        for (int j=0;j<rhs.DIM;j++){
            tmp[i] += lhs(i,j) * rhs(j);
        }
    }
    return MaterialVector(tmp.data());
}

inline MaterialVector operator* (const KinematicTensor& lhs, const MaterialVector& rhs){
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
inline MaterialVector operator/ (const MaterialVector& lhs, const double& rhs){
    MaterialVector tmp = MaterialVector(lhs);
    return tmp/=rhs;
}

inline MaterialVector operator/ (const MaterialVector& lhs, const int& rhs){
    MaterialVector tmp = MaterialVector(lhs);
    return tmp/=rhs;
}

/*----------------------------------------------------------------------------*/
//Multiplication (returning KinematicTensor)
inline KinematicVector operator* (const double& lhs, const KinematicVector& rhs){
    KinematicVector tmp = KinematicVector(rhs);
    return tmp*=lhs;
}

inline KinematicVector operator* (const int& lhs, const KinematicVector& rhs){
    KinematicVector tmp = KinematicVector(rhs);
    return tmp*=lhs;
}

inline KinematicVector operator* (const KinematicVector& lhs, const double& rhs){
    KinematicVector tmp = KinematicVector(lhs);
    return tmp*=rhs;
}

inline KinematicVector operator* (const KinematicVector& lhs, const int& rhs){
    KinematicVector tmp = KinematicVector(lhs);
    return tmp*=rhs;
}

inline KinematicVector operator* (const KinematicTensor& lhs, const KinematicVector& rhs){
    assert(lhs.TENSOR_TYPE == rhs.VECTOR_TYPE && "Multiplication failed.");
    std::array<double,MPMVector::VECTOR_MAX_DIM> tmp = { {0} };
    for (int i=0;i<lhs.DIM;i++){
        for (int j=0;j<rhs.DIM;j++){
            tmp[i] += lhs(i,j) * rhs(j);
        }
    }
    return KinematicVector(tmp.data(), rhs.VECTOR_TYPE);
}

inline KinematicVector operator/ (const KinematicVector& lhs, const double& rhs){
    KinematicVector tmp = KinematicVector(lhs);
    return tmp/=rhs;
}
inline KinematicVector operator/ (const KinematicVector& lhs, const int& rhs){
    KinematicVector tmp = KinematicVector(lhs);
    return tmp/=rhs;
}


/*----------------------------------------------------------------------------*/
//Vector Maps
inline MaterialVector::Map& MaterialVector::Map::operator= (const MaterialVector &other){
    for(int i=0;i<VECTOR_MAX_DIM;i++){data_ptr[i]=other(i);};
    return *this;
}

inline MaterialVector::Map& MaterialVector::Map::operator= (const KinematicVector &other){
    for(int i=0;i<VECTOR_MAX_DIM;i++){data_ptr[i]=other(i);};
    return *this;
}

inline KinematicVector::Map& KinematicVector::Map::operator= (const KinematicVector &other){
    for(int i=0;i<VECTOR_MAX_DIM;i++){data_ptr[i]=other(i);};
    return *this;
}

#endif //MPM_V3_MPMVECTOR_HPP
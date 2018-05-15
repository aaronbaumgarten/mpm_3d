//
// Created by aaron on 4/23/18.
// mpmvectorarray.hpp
//

#ifndef MPM_V3_MPMVECTORARRAY_HPP
#define MPM_V3_MPMVECTORARRAY_HPP

#include <stdlib.h>
#include "mpm_tensor.hpp"
#include "mpm_vector.hpp"

/*----------------------------------------------------------------------------*/
//base class for KinematicVectorArray and MaterialVectorArray
class MPMVectorArray {
public:
    //default construction
    //MPMVectorArray(){}

    //construction with initial length (data stored in single vector)
    //MPMVectorArray(int i){ buffer.resize(MPMVector::VECTOR_MAX_DIM * i); }

    /*------------------------------------------------------------------------*/
    //access j(th) component of i(th) vector
    virtual double& operator() (int i, int j) = 0;

    //resize data to size of i tensors
    void resize(int i){ buffer.resize(MPMVector::VECTOR_MAX_DIM * i); return; }

    //clear data
    void clear(){ buffer.clear(); return; }

    //set all elements to zero
    void setZero(){
        for(int i=0;i<buffer.size();i++){
            buffer[i] = 0;
        }
        return;
    }

    //return buffer size/9 (number of tensors)
    int size() const { return (buffer.size()/MPMVector::VECTOR_MAX_DIM); }

    //return pointer to data
    double* data(){ return buffer.data(); }

    //protected data member
protected:
    std::vector<double> buffer;
};

class KinematicVectorArray; //forward declare class

/*----------------------------------------------------------------------------*/
//storage unit for MaterialVectors
class MaterialVectorArray : public MPMVectorArray{
public:
    //constructors
    MaterialVectorArray(){}

    //construct with initial size
    MaterialVectorArray(int i){ buffer.resize(MPMVector::VECTOR_MAX_DIM * i); }

    /*------------------------------------------------------------------------*/
    //access j(th) component of i(th) vector
    double& operator() (int i, int j);

    //MaterialTensor::Map to data stored in i(th) tensor
    MaterialVector::Map operator() (int i);

    //MaterialTensor::Map to data stored in i(th) tensor
    MaterialVector::Map operator[] (int i);

    //effectively add one MaterialVector
    void push_back(MaterialVector& vector);

    //effectively add one KinematicVector
    void push_back(KinematicVector& vector);

    friend class KinematicVectorArray;

    //addition, subtraction, and scaling
    MaterialVectorArray operator-();
    MaterialVectorArray& operator+= (const MaterialVectorArray &rhs);
    MaterialVectorArray& operator-= (const MaterialVectorArray &rhs);
    MaterialVectorArray& operator+= (const KinematicVectorArray &rhs);
    MaterialVectorArray& operator-= (const KinematicVectorArray &rhs);
    MaterialVectorArray& operator*= (const int &rhs);
    MaterialVectorArray& operator*= (const double &rhs);
    MaterialVectorArray& operator/= (const int &rhs);
    MaterialVectorArray& operator/= (const double &rhs);
};

/*----------------------------------------------------------------------------*/
//storage unit for KinematicVectors
class KinematicVectorArray : public MPMVectorArray{
public:
    //default dimensions for KinematicTensors
    int DIM = MPMVector::VECTOR_MAX_DIM;
    int VECTOR_TYPE = KinematicVector::VECTOR_3D;

    //assign KinematicTensor type by input (rely on dummy KinematicTensor)
    void assignVectorType(int input);

    //constructors
    KinematicVectorArray(){}

    //construct with initial size
    KinematicVectorArray(int i){ buffer.resize(MPMVector::VECTOR_MAX_DIM * i); }

    //construct with specified type
    KinematicVectorArray(int i, int input);

    /*------------------------------------------------------------------------*/
    //access j(th) component of i(th) vector
    double& operator() (int i, int j);

    //KinematicTensor::Map to data stored in i(th) vector
    KinematicVector::Map operator() (int i);

    //KinematicTensor::Map to data stored in i(th) vector
    KinematicVector::Map operator[] (int i);

    //effectively add on KinematicVector
    void push_back(KinematicVector& vector);

    friend class MaterialVectorArray;

    KinematicVectorArray operator-();
    KinematicVectorArray& operator+= (const KinematicVectorArray &rhs);
    KinematicVectorArray& operator-= (const KinematicVectorArray &rhs);
    KinematicVectorArray& operator*= (const int &rhs);
    KinematicVectorArray& operator*= (const double &rhs);
    KinematicVectorArray& operator/= (const int &rhs);
    KinematicVectorArray& operator/= (const double &rhs);
};

/*----------------------------------------------------------------------------*/
//math functions
inline MaterialVectorArray operator+ (const MaterialVectorArray&, const MaterialVectorArray&);
inline MaterialVectorArray operator+ (const MaterialVectorArray&, const KinematicVectorArray&);
inline MaterialVectorArray operator+ (const KinematicVectorArray&, const MaterialVectorArray&);
inline KinematicVectorArray operator+ (const KinematicVectorArray&, const KinematicVectorArray&);

inline MaterialVectorArray operator- (const MaterialVectorArray&, const MaterialVectorArray&);
inline MaterialVectorArray operator- (const MaterialVectorArray&, const KinematicVectorArray&);
inline MaterialVectorArray operator- (const KinematicVectorArray&, const MaterialVectorArray&);
inline KinematicVectorArray operator- (const KinematicVectorArray&, const KinematicVectorArray&);

inline MaterialVectorArray operator* (const MaterialVectorArray&, const int&);
inline MaterialVectorArray operator* (const MaterialVectorArray&, const double&);
inline MaterialVectorArray operator* (const int&, const MaterialVectorArray&);
inline MaterialVectorArray operator* (const double&, const MaterialVectorArray&);

inline KinematicVectorArray operator* (const KinematicVectorArray&, const int&);
inline KinematicVectorArray operator* (const KinematicVectorArray&, const double&);
inline KinematicVectorArray operator* (const int&, const KinematicVectorArray&);
inline KinematicVectorArray operator* (const double&, const KinematicVectorArray&);

inline MaterialVectorArray operator/ (const MaterialVectorArray&, const int&);
inline MaterialVectorArray operator/ (const MaterialVectorArray&, const double&);

inline KinematicVectorArray operator/ (const KinematicVectorArray&, const int&);
inline KinematicVectorArray operator/ (const KinematicVectorArray&, const double&);

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
//access j(th) component of i(th) vector
inline double& MaterialVectorArray::operator() (int i, int j){
    assert((MPMVector::VECTOR_MAX_DIM * i + j) < buffer.size());
    return buffer[MPMVector::VECTOR_MAX_DIM * i + j];
}

//MaterialTensor::Map to data stored in i(th) tensor
inline MaterialVector::Map MaterialVectorArray::operator() (int i){
    assert((MPMVector::VECTOR_MAX_DIM * (i + 1) <= buffer.size()));
    return MaterialVector::Map(buffer.data() + MPMVector::VECTOR_MAX_DIM * i);
}

//MaterialTensor::Map to data stored in i(th) tensor
inline MaterialVector::Map MaterialVectorArray::operator[] (int i){
    assert((MPMVector::VECTOR_MAX_DIM * (i + 1) <= buffer.size()));
    return MaterialVector::Map(buffer.data() + MPMVector::VECTOR_MAX_DIM * i);
}

/*----------------------------------------------------------------------------*/
//add one MaterialVector to the end of MaterialVectorArray
inline void MaterialVectorArray::push_back(MaterialVector& vector){
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++) {
        buffer.push_back(vector[i]);
    }
    return;
}

//add one KinematicVector to the end of MaterialVectorArray
inline void MaterialVectorArray::push_back(KinematicVector& vector){
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++) {
        buffer.push_back(vector[i]);
    }
    return;
}

/*----------------------------------------------------------------------------*/
//simple math operators for the MaterialVectorArray
inline MaterialVectorArray MaterialVectorArray::operator-(){
    MaterialVectorArray tmp(size());
    for(int i=0;i<buffer.size();i++){
        tmp.buffer[i] = -buffer[i];
    }
    return tmp;
}

inline MaterialVectorArray& MaterialVectorArray::operator+= (const MaterialVectorArray &rhs){
    assert(size() == rhs.size() && "MaterialVectorArray addition failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] + rhs.buffer[i];
    }
    return *this;
}

inline MaterialVectorArray& MaterialVectorArray::operator-= (const MaterialVectorArray &rhs){
    assert(size() == rhs.size() && "MaterialVectorArray subtraction failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] - rhs.buffer[i];
    }
    return *this;
}

inline MaterialVectorArray& MaterialVectorArray::operator+= (const KinematicVectorArray &rhs){
    assert(size() == rhs.size() && "MaterialVectorArray addition failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] + rhs.buffer[i];
    }
    return *this;
}

inline MaterialVectorArray& MaterialVectorArray::operator-= (const KinematicVectorArray &rhs){
    assert(size() == rhs.size() && "MaterialVectorArray subtraction failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] - rhs.buffer[i];
    }
    return *this;
}

inline MaterialVectorArray& MaterialVectorArray::operator*= (const int &rhs){
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] * rhs;
    }
    return *this;
}

inline MaterialVectorArray& MaterialVectorArray::operator*= (const double &rhs){
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] * rhs;
    }
    return *this;
}

inline MaterialVectorArray& MaterialVectorArray::operator/= (const int &rhs){
    double tmp = 1.0/rhs;
    for (int i=0;i<buffer.size();i++){
        //buffer[i] = buffer[i] / rhs;
        buffer[i] = buffer[i] * tmp;
    }
    return *this;
}

inline MaterialVectorArray& MaterialVectorArray::operator/= (const double &rhs){
    double tmp = 1.0/rhs;
    for (int i=0;i<buffer.size();i++){
        //buffer[i] = buffer[i] / rhs;
        buffer[i] = buffer[i] * tmp;
    }
    return *this;
}


/*----------------------------------------------------------------------------*/
//assign vector type for KinematicVectorArray by input (rely on dummy KinematicVector)
inline void KinematicVectorArray::assignVectorType(int input){
    KinematicVector tmp(input);
    VECTOR_TYPE = tmp.VECTOR_TYPE;
    DIM = tmp.DIM;
}

/*----------------------------------------------------------------------------*/
//construct KinematicVectorArray with specified tensor type
inline KinematicVectorArray::KinematicVectorArray(int i, int input){
    assignVectorType(input);
    buffer.resize(MPMVector::VECTOR_MAX_DIM * i);
}

/*------------------------------------------------------------------------*/
//access j(th) component of i(th) KinematicVector
inline double& KinematicVectorArray::operator() (int i, int j){
    assert((MPMVector::VECTOR_MAX_DIM * i + j) < buffer.size());
    return buffer[MPMVector::VECTOR_MAX_DIM * i + j];
}

//KinematicTensor::Map to data stored in i(th) vector
inline KinematicVector::Map KinematicVectorArray::operator() (int i){
    assert(MPMVector::VECTOR_MAX_DIM * (i+1) <= buffer.size());
    return KinematicVector::Map(buffer.data() + MPMVector::VECTOR_MAX_DIM * i, VECTOR_TYPE);
}

//KinematicTensor::Map to data stored in i(th) vector
inline KinematicVector::Map KinematicVectorArray::operator[] (int i){
    assert(MPMVector::VECTOR_MAX_DIM * (i+1) <= buffer.size());
    return KinematicVector::Map(buffer.data() + MPMVector::VECTOR_MAX_DIM * i, VECTOR_TYPE);
}

//effectively add on KinematicVector
inline void KinematicVectorArray::push_back(KinematicVector& vector){
    assert(VECTOR_TYPE == vector.VECTOR_TYPE && "Insert failed.");
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++) {
        buffer.push_back(vector[i]);
    }
    return;
}

/*----------------------------------------------------------------------------*/
//simple math operators for the MaterialVectorArray
inline KinematicVectorArray KinematicVectorArray::operator-(){
    KinematicVectorArray tmp(size(),VECTOR_TYPE);
    for(int i=0;i<buffer.size();i++){
        tmp.buffer[i] = -buffer[i];
    }
    return tmp;
}

inline KinematicVectorArray& KinematicVectorArray::operator+= (const KinematicVectorArray &rhs){
    assert(size() == rhs.size() && VECTOR_TYPE == rhs.VECTOR_TYPE && "KinematicVectorArray addition failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] + rhs.buffer[i];
    }
    return *this;
}

inline KinematicVectorArray& KinematicVectorArray::operator-= (const KinematicVectorArray &rhs){
    assert(size() == rhs.size() && VECTOR_TYPE == rhs.VECTOR_TYPE && "KinematicVectorArray addition failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] - rhs.buffer[i];
    }
    return *this;
}

inline KinematicVectorArray& KinematicVectorArray::operator*= (const int &rhs){
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] * rhs;
    }
    return *this;
}

inline KinematicVectorArray& KinematicVectorArray::operator*= (const double &rhs){
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] * rhs;
    }
    return *this;
}

inline KinematicVectorArray& KinematicVectorArray::operator/= (const int &rhs){
    double tmp = 1.0/rhs;
    for (int i=0;i<buffer.size();i++){
        //buffer[i] = buffer[i] / rhs;
        buffer[i] = buffer[i] * tmp;
    }
    return *this;
}

inline KinematicVectorArray& KinematicVectorArray::operator/= (const double &rhs){
    double tmp = 1.0/rhs;
    for (int i=0;i<buffer.size();i++){
        //buffer[i] = buffer[i] / rhs;
        buffer[i] = buffer[i] * tmp;
    }
    return *this;
}


/*----------------------------------------------------------------------------*/
//math functions
inline MaterialVectorArray operator+ (const MaterialVectorArray& lhs, const MaterialVectorArray& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(lhs);
    return tmp+=rhs;
}

inline MaterialVectorArray operator+ (const MaterialVectorArray& lhs, const KinematicVectorArray& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(lhs);
    return tmp+=rhs;
}

inline MaterialVectorArray operator+ (const KinematicVectorArray& lhs, const MaterialVectorArray& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(rhs);
    return tmp+=lhs;
}

inline KinematicVectorArray operator+ (const KinematicVectorArray& lhs, const KinematicVectorArray& rhs){
    KinematicVectorArray tmp = KinematicVectorArray(lhs);
    return tmp+=rhs;
}

inline MaterialVectorArray operator- (const MaterialVectorArray& lhs, const MaterialVectorArray& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(lhs);
    return tmp-=rhs;
}

inline MaterialVectorArray operator- (const MaterialVectorArray& lhs, const KinematicVectorArray& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(lhs);
    return tmp-=rhs;
}

inline MaterialVectorArray operator- (const KinematicVectorArray& lhs, const MaterialVectorArray& rhs){
    return -(rhs-lhs);
}

inline KinematicVectorArray operator- (const KinematicVectorArray& lhs, const KinematicVectorArray& rhs){
    KinematicVectorArray tmp = KinematicVectorArray(lhs);
    return tmp-=rhs;
}


inline MaterialVectorArray operator* (const MaterialVectorArray& lhs, const int& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(lhs);
    return tmp*=rhs;
}

inline MaterialVectorArray operator* (const MaterialVectorArray& lhs, const double& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(lhs);
    return tmp*=rhs;
}

inline MaterialVectorArray operator* (const int& lhs, const MaterialVectorArray& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(rhs);
    return tmp*=lhs;
}

inline MaterialVectorArray operator* (const double& lhs, const MaterialVectorArray& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(rhs);
    return tmp*=lhs;
}

inline KinematicVectorArray operator* (const KinematicVectorArray& lhs, const int& rhs){
    KinematicVectorArray tmp = KinematicVectorArray(lhs);
    return tmp*=rhs;
}

inline KinematicVectorArray operator* (const KinematicVectorArray& lhs, const double& rhs){
    KinematicVectorArray tmp = KinematicVectorArray(lhs);
    return tmp*=rhs;
}

inline KinematicVectorArray operator* (const int& lhs, const KinematicVectorArray& rhs){
    KinematicVectorArray tmp = KinematicVectorArray(rhs);
    return tmp*=lhs;
}

inline KinematicVectorArray operator* (const double& lhs, const KinematicVectorArray& rhs){
    KinematicVectorArray tmp = KinematicVectorArray(rhs);
    return tmp*=lhs;
}

inline MaterialVectorArray operator/ (const MaterialVectorArray& lhs, const int& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(lhs);
    return tmp/=rhs;
}

inline MaterialVectorArray operator/ (const MaterialVectorArray& lhs, const double& rhs){
    MaterialVectorArray tmp = MaterialVectorArray(lhs);
    return tmp/=rhs;
}

inline KinematicVectorArray operator/ (const KinematicVectorArray& lhs, const int& rhs){
    KinematicVectorArray tmp = KinematicVectorArray(lhs);
    return tmp/=rhs;
}

inline KinematicVectorArray operator/ (const KinematicVectorArray& lhs, const double& rhs){
    KinematicVectorArray tmp = KinematicVectorArray(lhs);
    return tmp/=rhs;
}

#endif //MPM_V3_MPMVECTORARRAY_HPP

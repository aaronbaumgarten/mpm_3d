//
// Created by aaron on 4/22/18.
// mpmtensorarray.hpp
//

#ifndef MPM_V3_MPMTENSORARRAY_HPP
#define MPM_V3_MPMTENSORARRAY_HPP

#include <stdlib.h>
#include "mpm_tensor.hpp"

/*----------------------------------------------------------------------------*/
//base class for KinematicTensorArray and MaterialTensorArray
class MPMTensorArray {
public:
    //default construction
    //MPMTensorArray(){}

    //construction with initial length (data stored in single vector)
    //MPMTensorArray(int i){ buffer.resize(MPMTensor::TENSOR_MAX_LENGTH*i); }

    /*------------------------------------------------------------------------*/
    //access j,k(th) component of i(th) tensor
    virtual double& operator() (int i, int j, int k) = 0;

    //access j(th) component of i(th) tensor
    virtual double& operator() (int i, int j) = 0;

    //resize data to size of i tensors
    void resize(int i){ buffer.resize(MPMTensor::TENSOR_MAX_LENGTH*i); return; }

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
    int size() const { return (buffer.size()/MPMTensor::TENSOR_MAX_LENGTH); }

    //return pointer to data
    double* data(){ return buffer.data(); }

    //protected data member
protected:
    std::vector<double> buffer;
};

class KinematicTensorArray; //forward declare class

/*----------------------------------------------------------------------------*/
//storage unit for MaterialTensors
class MaterialTensorArray : public MPMTensorArray{
public:
    //constructors
    MaterialTensorArray(){}

    //construct with initial size
    MaterialTensorArray(int i){ buffer.resize(MPMTensor::TENSOR_MAX_LENGTH*i); }

    /*------------------------------------------------------------------------*/
    //access j,k(th) component of i(th) tensor
    double& operator() (int i, int j, int k);

    //access j(th) component of i(th) tensor
    double& operator() (int i, int j);

    //MaterialTensor::Map to data stored in i(th) tensor
    MaterialTensor::Map operator() (int i) const;

    //MaterialTensor::Map to data stored in i(th) tensor
    MaterialTensor::Map operator[] (int i) const;

    //effectively add one MaterialTensor
    void push_back(const MaterialTensor& tensor);

    //effectively add one KinematicTensor
    void push_back(const KinematicTensor& tensor);

    friend class KinematicTensorArray;

    //addition, subtraction, and scaling
    MaterialTensorArray operator-();
    MaterialTensorArray& operator+= (const MaterialTensorArray &rhs);
    MaterialTensorArray& operator-= (const MaterialTensorArray &rhs);
    MaterialTensorArray& operator+= (const KinematicTensorArray &rhs);
    MaterialTensorArray& operator-= (const KinematicTensorArray &rhs);
    MaterialTensorArray& operator*= (const int &rhs);
    MaterialTensorArray& operator*= (const double &rhs);
    MaterialTensorArray& operator/= (const int &rhs);
    MaterialTensorArray& operator/= (const double &rhs);
};

/*----------------------------------------------------------------------------*/
//storage unit for KinematicTensors
class KinematicTensorArray : public MPMTensorArray{
public:
    //default dimensions for KinematicTensors
    int DIM = MPMTensor::TENSOR_MAX_DIM;
    int TENSOR_TYPE = KinematicTensor::TENSOR_3D;

    //assign KinematicTensor type by input (rely on dummy KinematicTensor)
    void assignTensorType(int input);

    //default constuctors (assume 3D)
    KinematicTensorArray(){}
    KinematicTensorArray(int i){ buffer.resize(MPMTensor::TENSOR_MAX_LENGTH*i); }

    //construct with specified type
    KinematicTensorArray(int i, int input);

    /*------------------------------------------------------------------------*/
    //access j,k(th) component of i(th) tensor
    double& operator() (int i, int j, int k);

    //access j(th) component of i(th) tensor
    double& operator() (int i, int j);

    //KinematicTensor::Map to data stored in i(th) tensor
    KinematicTensor::Map operator() (int i) const;

    //KinematicTensor::Map to data stored in i(th) tensor
    KinematicTensor::Map operator[] (int i) const;

    //effectively add on KinematicTensor
    void push_back(const KinematicTensor& tensor);

    friend class MaterialTensorArray;

    //addition subtraction and scaling
    KinematicTensorArray operator-();
    KinematicTensorArray& operator+= (const KinematicTensorArray &rhs);
    KinematicTensorArray& operator-= (const KinematicTensorArray &rhs);
    KinematicTensorArray& operator*= (const int &rhs);
    KinematicTensorArray& operator*= (const double &rhs);
    KinematicTensorArray& operator/= (const int &rhs);
    KinematicTensorArray& operator/= (const double &rhs);
};

/*----------------------------------------------------------------------------*/
//math functions
inline MaterialTensorArray operator+ (const MaterialTensorArray&, const MaterialTensorArray&);
inline MaterialTensorArray operator+ (const MaterialTensorArray&, const KinematicTensorArray&);
inline MaterialTensorArray operator+ (const KinematicTensorArray&, const MaterialTensorArray&);
inline KinematicTensorArray operator+ (const KinematicTensorArray&, const KinematicTensorArray&);

inline MaterialTensorArray operator- (const MaterialTensorArray&, const MaterialTensorArray&);
inline MaterialTensorArray operator- (const MaterialTensorArray&, const KinematicTensorArray&);
inline MaterialTensorArray operator- (const KinematicTensorArray&, const MaterialTensorArray&);
inline KinematicTensorArray operator- (const KinematicTensorArray&, const KinematicTensorArray&);

inline MaterialTensorArray operator* (const MaterialTensorArray&, const int&);
inline MaterialTensorArray operator* (const MaterialTensorArray&, const double&);
inline MaterialTensorArray operator* (const int&, const MaterialTensorArray&);
inline MaterialTensorArray operator* (const double&, const MaterialTensorArray&);

inline KinematicTensorArray operator* (const KinematicTensorArray&, const int&);
inline KinematicTensorArray operator* (const KinematicTensorArray&, const double&);
inline KinematicTensorArray operator* (const int&, const KinematicTensorArray&);
inline KinematicTensorArray operator* (const double&, const KinematicTensorArray&);

inline MaterialTensorArray operator/ (const MaterialTensorArray&, const int&);
inline MaterialTensorArray operator/ (const MaterialTensorArray&, const double&);

inline KinematicTensorArray operator/ (const KinematicTensorArray&, const int&);
inline KinematicTensorArray operator/ (const KinematicTensorArray&, const double&);

/*----------------------------------------------------------------------------*/
//access j,k(th) component of i(th) MaterialTensor in MaterialTensorArray
inline double& MaterialTensorArray::operator() (int i, int j, int k) {
    assert((MPMTensor::TENSOR_MAX_LENGTH * i + MPMTensor::TENSOR_MAX_DIM * j + k) < buffer.size());
    return buffer[MPMTensor::TENSOR_MAX_LENGTH * i + MPMTensor::TENSOR_MAX_DIM * j + k];
}

//access j(th) component of i(th) MaterialTensor in MaterialTensorArray (row major)
inline double& MaterialTensorArray::operator() (int i, int j) {
    assert((MPMTensor::TENSOR_MAX_LENGTH * i + j) < buffer.size());
    return buffer[MPMTensor::TENSOR_MAX_LENGTH * i + j];
}

//return MaterialTensor::Map to data stored in i(th) MaterialTensor
inline MaterialTensor::Map MaterialTensorArray::operator() (int i) const {
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return MaterialTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i);
}

//return MaterialTensor::Map to data stored in i(th) MaterialTensor
inline MaterialTensor::Map MaterialTensorArray::operator[] (int i) const {
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return MaterialTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i);
}

/*----------------------------------------------------------------------------*/
//add one MaterialTensor to the end of MaterialTensorArray
inline void MaterialTensorArray::push_back(const MaterialTensor& tensor){
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++) {
        buffer.push_back(tensor[i]);
    }
    return;
}

//add one KinematicTensor to the end of MaterialTensorArray
inline void MaterialTensorArray::push_back(const KinematicTensor& tensor){
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++) {
        buffer.push_back(tensor[i]);
    }
    return;
}

/*----------------------------------------------------------------------------*/
//simple math operators for the MaterialTensorArray
inline MaterialTensorArray MaterialTensorArray::operator-(){
    MaterialTensorArray tmp(size());
    for(int i=0;i<buffer.size();i++){
        tmp.buffer[i] = -buffer[i];
    }
    return tmp;
}

inline MaterialTensorArray& MaterialTensorArray::operator+= (const MaterialTensorArray &rhs){
    assert(size() == rhs.size() && "MaterialTensorArray addition failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] + rhs.buffer[i];
    }
    return *this;
}

inline MaterialTensorArray& MaterialTensorArray::operator-= (const MaterialTensorArray &rhs){
    assert(size() == rhs.size() && "MaterialTensorArray subtraction failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] - rhs.buffer[i];
    }
    return *this;
}

inline MaterialTensorArray& MaterialTensorArray::operator+= (const KinematicTensorArray &rhs){
    assert(size() == rhs.size() && "MaterialTensorArray addition failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] + rhs.buffer[i];
    }
    return *this;
}

inline MaterialTensorArray& MaterialTensorArray::operator-= (const KinematicTensorArray &rhs){
    assert(size() == rhs.size() && "MaterialTensorArray subtraction failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] - rhs.buffer[i];
    }
    return *this;
}

inline MaterialTensorArray& MaterialTensorArray::operator*= (const int &rhs){
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] * rhs;
    }
    return *this;
}

inline MaterialTensorArray& MaterialTensorArray::operator*= (const double &rhs){
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] * rhs;
    }
    return *this;
}

inline MaterialTensorArray& MaterialTensorArray::operator/= (const int &rhs){
    //double tmp = 1.0/rhs;
    for (int i=0;i<buffer.size();i++){
        //buffer[i] = buffer[i] / rhs;
        buffer[i] = buffer[i] / rhs;
    }
    return *this;
}

inline MaterialTensorArray& MaterialTensorArray::operator/= (const double &rhs){
    //double tmp = 1.0/rhs;
    for (int i=0;i<buffer.size();i++){
        //buffer[i] = buffer[i] / rhs;
        buffer[i] = buffer[i] / rhs;
    }
    return *this;
}


/*----------------------------------------------------------------------------*/
//assign tensor type for KinematicTensorArray by input (rely on dummy KinematicTensor)
inline void KinematicTensorArray::assignTensorType(int input){
    KinematicTensor tmp(input);
    TENSOR_TYPE = tmp.TENSOR_TYPE;
    DIM = tmp.DIM;
}

/*----------------------------------------------------------------------------*/
//construct KinematicTensorArray with specified tensor type
inline KinematicTensorArray::KinematicTensorArray(int i, int input){
    assignTensorType(input);
    buffer.resize(MPMTensor::TENSOR_MAX_LENGTH*i);
}

/*----------------------------------------------------------------------------*/
//access j,k(th) component of i(th) KinematicTensor in KinematicTensorArray
inline double& KinematicTensorArray::operator() (int i, int j, int k) {
    assert((MPMTensor::TENSOR_MAX_LENGTH * i + MPMTensor::TENSOR_MAX_DIM * j + k) < buffer.size());
    return buffer[MPMTensor::TENSOR_MAX_LENGTH * i + MPMTensor::TENSOR_MAX_DIM * j + k];
}

//access j(th) component of i(th) KinematicTensor in KinematicTensorArray (row major)
inline double& KinematicTensorArray::operator() (int i, int j) {
    assert((MPMTensor::TENSOR_MAX_LENGTH*i+j) < buffer.size());
    return buffer[MPMTensor::TENSOR_MAX_LENGTH*i + j];
}

//return KinematicTensor::Map to data stored in i(th) KinematicTensor
inline KinematicTensor::Map KinematicTensorArray::operator() (int i) const {
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return KinematicTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i, TENSOR_TYPE);
}

//same as previous, but with [] instead of ()
inline KinematicTensor::Map KinematicTensorArray::operator[] (int i) const {
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return KinematicTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i, TENSOR_TYPE);
}

//add one KinematicTensor to end of KinematicTensorArray
inline void KinematicTensorArray::push_back(const KinematicTensor& tensor){
    assert(TENSOR_TYPE == tensor.TENSOR_TYPE && "Insert failed.");
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++) {
        buffer.push_back(tensor[i]);
    }
    return;
}

/*----------------------------------------------------------------------------*/
//simple math operators for the MaterialTensorArray
inline KinematicTensorArray KinematicTensorArray::operator-(){
    KinematicTensorArray tmp(size(),TENSOR_TYPE);
    for(int i=0;i<buffer.size();i++){
        tmp.buffer[i] = -buffer[i];
    }
    return tmp;
}

inline KinematicTensorArray& KinematicTensorArray::operator+= (const KinematicTensorArray &rhs){
    assert(size() == rhs.size() && TENSOR_TYPE == rhs.TENSOR_TYPE && "KinematicTensorArray addition failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] + rhs.buffer[i];
    }
    return *this;
}

inline KinematicTensorArray& KinematicTensorArray::operator-= (const KinematicTensorArray &rhs){
    assert(size() == rhs.size() && TENSOR_TYPE == rhs.TENSOR_TYPE && "KinematicTensorArray addition failed.");
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] - rhs.buffer[i];
    }
    return *this;
}

inline KinematicTensorArray& KinematicTensorArray::operator*= (const int &rhs){
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] * rhs;
    }
    return *this;
}

inline KinematicTensorArray& KinematicTensorArray::operator*= (const double &rhs){
    for (int i=0;i<buffer.size();i++){
        buffer[i] = buffer[i] * rhs;
    }
    return *this;
}

inline KinematicTensorArray& KinematicTensorArray::operator/= (const int &rhs){
    //double tmp = 1.0/rhs;
    for (int i=0;i<buffer.size();i++){
        //buffer[i] = buffer[i] / rhs;
        buffer[i] = buffer[i] / rhs;
    }
    return *this;
}

inline KinematicTensorArray& KinematicTensorArray::operator/= (const double &rhs){
    //double tmp = 1.0/rhs;
    for (int i=0;i<buffer.size();i++){
        //buffer[i] = buffer[i] / rhs;
        buffer[i] = buffer[i] / rhs;
    }
    return *this;
}

/*----------------------------------------------------------------------------*/
//math functions
inline MaterialTensorArray operator+ (const MaterialTensorArray& lhs, const MaterialTensorArray& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(lhs);
    return tmp+=rhs;
}

inline MaterialTensorArray operator+ (const MaterialTensorArray& lhs, const KinematicTensorArray& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(lhs);
    return tmp+=rhs;
}

inline MaterialTensorArray operator+ (const KinematicTensorArray& lhs, const MaterialTensorArray& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(rhs);
    return tmp+=lhs;
}

inline KinematicTensorArray operator+ (const KinematicTensorArray& lhs, const KinematicTensorArray& rhs){
    KinematicTensorArray tmp = KinematicTensorArray(lhs);
    return tmp+=rhs;
}

inline MaterialTensorArray operator- (const MaterialTensorArray& lhs, const MaterialTensorArray& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(lhs);
    return tmp-=rhs;
}

inline MaterialTensorArray operator- (const MaterialTensorArray& lhs, const KinematicTensorArray& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(lhs);
    return tmp-=rhs;
}

inline MaterialTensorArray operator- (const KinematicTensorArray& lhs, const MaterialTensorArray& rhs){
    return -(rhs-lhs);
}

inline KinematicTensorArray operator- (const KinematicTensorArray& lhs, const KinematicTensorArray& rhs){
    KinematicTensorArray tmp = KinematicTensorArray(lhs);
    return tmp-=rhs;
}


inline MaterialTensorArray operator* (const MaterialTensorArray& lhs, const int& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(lhs);
    return tmp*=rhs;
}

inline MaterialTensorArray operator* (const MaterialTensorArray& lhs, const double& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(lhs);
    return tmp*=rhs;
}

inline MaterialTensorArray operator* (const int& lhs, const MaterialTensorArray& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(rhs);
    return tmp*=lhs;
}

inline MaterialTensorArray operator* (const double& lhs, const MaterialTensorArray& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(rhs);
    return tmp*=lhs;
}

inline KinematicTensorArray operator* (const KinematicTensorArray& lhs, const int& rhs){
    KinematicTensorArray tmp = KinematicTensorArray(lhs);
    return tmp*=rhs;
}

inline KinematicTensorArray operator* (const KinematicTensorArray& lhs, const double& rhs){
    KinematicTensorArray tmp = KinematicTensorArray(lhs);
    return tmp*=rhs;
}

inline KinematicTensorArray operator* (const int& lhs, const KinematicTensorArray& rhs){
    KinematicTensorArray tmp = KinematicTensorArray(rhs);
    return tmp*=lhs;
}

inline KinematicTensorArray operator* (const double& lhs, const KinematicTensorArray& rhs){
    KinematicTensorArray tmp = KinematicTensorArray(rhs);
    return tmp*=lhs;
}

inline MaterialTensorArray operator/ (const MaterialTensorArray& lhs, const int& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(lhs);
    return tmp/=rhs;
}

inline MaterialTensorArray operator/ (const MaterialTensorArray& lhs, const double& rhs){
    MaterialTensorArray tmp = MaterialTensorArray(lhs);
    return tmp/=rhs;
}

inline KinematicTensorArray operator/ (const KinematicTensorArray& lhs, const int& rhs){
    KinematicTensorArray tmp = KinematicTensorArray(lhs);
    return tmp/=rhs;
}

inline KinematicTensorArray operator/ (const KinematicTensorArray& lhs, const double& rhs){
    KinematicTensorArray tmp = KinematicTensorArray(lhs);
    return tmp/=rhs;
}


#endif //MPM_V3_MPMTENSORARRAY_HPP

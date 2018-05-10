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
    int size(){ return (buffer.size()/MPMTensor::TENSOR_MAX_LENGTH); }

    //return pointer to data
    double* data(){ return buffer.data(); }

    //protected data member
protected:
    std::vector<double> buffer;
};

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
    MaterialTensor::Map operator() (int i);

    //MaterialTensor::Map to data stored in i(th) tensor
    MaterialTensor::Map operator[] (int i);

    //effectively add one MaterialTensor
    void push_back(MaterialTensor& tensor);

    //effectively add one KinematicTensor
    void push_back(KinematicTensor& tensor);
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
    KinematicTensor::Map operator() (int i);

    //KinematicTensor::Map to data stored in i(th) tensor
    KinematicTensor::Map operator[] (int i);

    //effectively add on KinematicTensor
    void push_back(KinematicTensor& tensor);
};


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
inline MaterialTensor::Map MaterialTensorArray::operator() (int i){
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return MaterialTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i);
}

//return MaterialTensor::Map to data stored in i(th) MaterialTensor
inline MaterialTensor::Map MaterialTensorArray::operator[] (int i){
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return MaterialTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i);
}

/*----------------------------------------------------------------------------*/
//add one MaterialTensor to the end of MaterialTensorArray
inline void MaterialTensorArray::push_back(MaterialTensor& tensor){
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++) {
        buffer.push_back(tensor[i]);
    }
    return;
}

//add one KinematicTensor to the end of MaterialTensorArray
inline void MaterialTensorArray::push_back(KinematicTensor& tensor){
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++) {
        buffer.push_back(tensor[i]);
    }
    return;
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
inline KinematicTensor::Map KinematicTensorArray::operator() (int i){
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return KinematicTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i, TENSOR_TYPE);
}

//same as previous, but with [] instead of ()
inline KinematicTensor::Map KinematicTensorArray::operator[] (int i){
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return KinematicTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i, TENSOR_TYPE);
}

//add one KinematicTensor to end of KinematicTensorArray
inline void KinematicTensorArray::push_back(KinematicTensor& tensor){
    assert(TENSOR_TYPE == tensor.TENSOR_TYPE && "Insert failed.");
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++) {
        buffer.push_back(tensor[i]);
    }
    return;
}


#endif //MPM_V3_MPMTENSORARRAY_HPP

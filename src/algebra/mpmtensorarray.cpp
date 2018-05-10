//
// Created by aaron on 4/23/18.
// mpmtensorarray.cpp
//

#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

/*----------------------------------------------------------------------------*/
//access j,k(th) component of i(th) MaterialTensor in MaterialTensorArray
double& MaterialTensorArray::operator() (int i, int j, int k) {
    assert((MPMTensor::TENSOR_MAX_LENGTH * i + MPMTensor::TENSOR_MAX_DIM * j + k) < buffer.size());
    return buffer[MPMTensor::TENSOR_MAX_LENGTH * i + MPMTensor::TENSOR_MAX_DIM * j + k];
}

//access j(th) component of i(th) MaterialTensor in MaterialTensorArray (row major)
double& MaterialTensorArray::operator() (int i, int j) {
    assert((MPMTensor::TENSOR_MAX_LENGTH * i + j) < buffer.size());
    return buffer[MPMTensor::TENSOR_MAX_LENGTH * i + j];
}

//return MaterialTensor::Map to data stored in i(th) MaterialTensor
MaterialTensor::Map MaterialTensorArray::operator() (int i){
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return MaterialTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i);
}

//return MaterialTensor::Map to data stored in i(th) MaterialTensor
MaterialTensor::Map MaterialTensorArray::operator[] (int i){
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return MaterialTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i);
}

/*----------------------------------------------------------------------------*/
//add one MaterialTensor to the end of MaterialTensorArray
void MaterialTensorArray::push_back(MaterialTensor& tensor){
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++) {
        buffer.push_back(tensor[i]);
    }
    return;
}

//add one KinematicTensor to the end of MaterialTensorArray
void MaterialTensorArray::push_back(KinematicTensor& tensor){
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++) {
        buffer.push_back(tensor[i]);
    }
    return;
}


/*----------------------------------------------------------------------------*/
//assign tensor type for KinematicTensorArray by input (rely on dummy KinematicTensor)
void KinematicTensorArray::assignTensorType(int input){
    KinematicTensor tmp(input);
    TENSOR_TYPE = tmp.TENSOR_TYPE;
    DIM = tmp.DIM;
}

/*----------------------------------------------------------------------------*/
//construct KinematicTensorArray with specified tensor type
KinematicTensorArray::KinematicTensorArray(int i, int input){
    assignTensorType(input);
    buffer.resize(MPMTensor::TENSOR_MAX_LENGTH*i);
}

/*----------------------------------------------------------------------------*/
//access j,k(th) component of i(th) KinematicTensor in KinematicTensorArray
double& KinematicTensorArray::operator() (int i, int j, int k) {
    assert((MPMTensor::TENSOR_MAX_LENGTH * i + MPMTensor::TENSOR_MAX_DIM * j + k) < buffer.size());
    return buffer[MPMTensor::TENSOR_MAX_LENGTH * i + MPMTensor::TENSOR_MAX_DIM * j + k];
}

//access j(th) component of i(th) KinematicTensor in KinematicTensorArray (row major)
double& KinematicTensorArray::operator() (int i, int j) {
    assert((MPMTensor::TENSOR_MAX_LENGTH*i+j) < buffer.size());
    return buffer[MPMTensor::TENSOR_MAX_LENGTH*i + j];
}

//return KinematicTensor::Map to data stored in i(th) KinematicTensor
KinematicTensor::Map KinematicTensorArray::operator() (int i){
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return KinematicTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i, TENSOR_TYPE);
}

//same as previous, but with [] instead of ()
KinematicTensor::Map KinematicTensorArray::operator[] (int i){
    assert(MPMTensor::TENSOR_MAX_LENGTH * (i+1) <= buffer.size());
    return KinematicTensor::Map(buffer.data() + MPMTensor::TENSOR_MAX_LENGTH * i, TENSOR_TYPE);
}

//add one KinematicTensor to end of KinematicTensorArray
void KinematicTensorArray::push_back(KinematicTensor& tensor){
    assert(TENSOR_TYPE == tensor.TENSOR_TYPE && "Insert failed.");
    for(int i=0;i<MPMTensor::TENSOR_MAX_LENGTH;i++) {
        buffer.push_back(tensor[i]);
    }
    return;
}
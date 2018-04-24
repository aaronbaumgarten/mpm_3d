//
// Created by aaron on 4/23/18.
// mpmvectorarray.hpp
//

#ifndef MPM_V3_MPMVECTORARRAY_HPP
#define MPM_V3_MPMVECTORARRAY_HPP

#include <stdlib.h>
#include "mpmtensor.hpp"
#include "mpmvector.hpp"

/*----------------------------------------------------------------------------*/
//base class for KinematicTensorArray and MaterialTensorArray
class MPMVectorArray {
public:
    //default construction
    //MPMVectorArray(){}

    //construction with initial length (data stored in single vector)
    //MPMVectorArray(int i){ buffer.resize(MPMVector::VECTOR_MAX_DIM * i); }

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
    int size(){ return (buffer.size()/MPMVector::VECTOR_MAX_DIM); }

    //return pointer to data
    double* data(){ return buffer.data(); }

    //protected data member
protected:
    std::vector<double> buffer;
};

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
};

#endif //MPM_V3_MPMVECTORARRAY_HPP

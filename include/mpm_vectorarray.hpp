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



#endif //MPM_V3_MPMVECTORARRAY_HPP

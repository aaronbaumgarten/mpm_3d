//
// Created by aaron on 4/23/18.
// mpmvectorarray.cpp
//

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"

/*----------------------------------------------------------------------------*/
//access j(th) component of i(th) vector
double& MaterialVectorArray::operator() (int i, int j){
    assert((MPMVector::VECTOR_MAX_DIM * i + j) < buffer.size());
    return buffer[MPMVector::VECTOR_MAX_DIM * i + j];
}

//MaterialTensor::Map to data stored in i(th) tensor
MaterialVector::Map MaterialVectorArray::operator() (int i){
    assert((MPMVector::VECTOR_MAX_DIM * (i + 1) <= buffer.size()));
    return MaterialVector::Map(buffer.data() + MPMVector::VECTOR_MAX_DIM * i);
}

//MaterialTensor::Map to data stored in i(th) tensor
MaterialVector::Map MaterialVectorArray::operator[] (int i){
    assert((MPMVector::VECTOR_MAX_DIM * (i + 1) <= buffer.size()));
    return MaterialVector::Map(buffer.data() + MPMVector::VECTOR_MAX_DIM * i);
}

/*----------------------------------------------------------------------------*/
//add one MaterialVector to the end of MaterialVectorArray
void MaterialVectorArray::push_back(MaterialVector& vector){
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++) {
        buffer.push_back(vector[i]);
    }
    return;
}

//add one KinematicVector to the end of MaterialVectorArray
void MaterialVectorArray::push_back(KinematicVector& vector){
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++) {
        buffer.push_back(vector[i]);
    }
    return;
}


/*----------------------------------------------------------------------------*/
//assign vector type for KinematicVectorArray by input (rely on dummy KinematicVector)
void KinematicVectorArray::assignVectorType(int input){
    KinematicVector tmp(input);
    VECTOR_TYPE = tmp.VECTOR_TYPE;
    DIM = tmp.DIM;
}

/*----------------------------------------------------------------------------*/
//construct KinematicVectorArray with specified tensor type
KinematicVectorArray::KinematicVectorArray(int i, int input){
    assignVectorType(input);
    buffer.resize(MPMVector::VECTOR_MAX_DIM * i);
}

/*------------------------------------------------------------------------*/
//access j(th) component of i(th) KinematicVector
double& KinematicVectorArray::operator() (int i, int j){
    assert((MPMVector::VECTOR_MAX_DIM * i + j) < buffer.size());
    return buffer[MPMVector::VECTOR_MAX_DIM * i + j];
}

//KinematicTensor::Map to data stored in i(th) vector
KinematicVector::Map KinematicVectorArray::operator() (int i){
    assert(MPMVector::VECTOR_MAX_DIM * (i+1) <= buffer.size());
    return KinematicVector::Map(buffer.data() + MPMVector::VECTOR_MAX_DIM * i, VECTOR_TYPE);
}

//KinematicTensor::Map to data stored in i(th) vector
KinematicVector::Map KinematicVectorArray::operator[] (int i){
    assert(MPMVector::VECTOR_MAX_DIM * (i+1) <= buffer.size());
    return KinematicVector::Map(buffer.data() + MPMVector::VECTOR_MAX_DIM * i, VECTOR_TYPE);
}

//effectively add on KinematicVector
void KinematicVectorArray::push_back(KinematicVector& vector){
    assert(VECTOR_TYPE == vector.VECTOR_TYPE && "Insert failed.");
    for(int i=0;i<MPMVector::VECTOR_MAX_DIM;i++) {
        buffer.push_back(vector[i]);
    }
    return;
}


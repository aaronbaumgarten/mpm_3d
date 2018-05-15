//
// Created by aaron on 5/12/18.
// mpm_sparse.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

#include "mpm_sparse.hpp"

/*
 * This file defines the 'sparse matrix' operators.
 * The 'sparse matrix' classes defined in mpm_sparse.hpp hold 'tuples' determining the index and
 * value of a given position in the 'matrix'. The operators defined in this file determine how
 * a 'sparse matrix' with M rows and N cols operates on an 'array' with N rows to produce another
 * 'array' with M rows.
 *
 * Intuitively, A 'operator' v is determined by (Av)_i = A_ij 'operator' v_j.
 */

/*----------------------------------------------------------------------------*/
//sparse matrix multiplication operators
Eigen::VectorXd operator* (MPMScalarSparseMatrix& lhs, Eigen::VectorXd& rhs){
    assert(lhs.cols() == rhs.rows() && "Scalar sparse matrix multiplication failed.");
    Eigen::VectorXd tmp = Eigen::VectorXd::Zero(lhs.rows());
    int i_tmp, j_tmp;
    for (int k=0;k<lhs.size();k++){
        i_tmp = lhs.i_vec[k];
        j_tmp = lhs.j_vec[k];
        tmp[i_tmp] += lhs.buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialVectorArray operator* (MaterialVectorSparseMatrix& lhs, Eigen::VectorXd& rhs){
    assert(lhs.cols() == rhs.rows() && "Vector sparse matrix multiplication failed.");
    MaterialVectorArray tmp = MaterialVectorArray(lhs.rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<lhs.size();k++){
        i_tmp = lhs.i_vec[k];
        j_tmp = lhs.j_vec[k];
        tmp[i_tmp] += lhs.buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

KinematicVectorArray operator* (KinematicVectorSparseMatrix& lhs, Eigen::VectorXd& rhs){
    assert(lhs.cols() == rhs.rows() && "Vector sparse matrix multiplication failed.");
    KinematicVectorArray tmp = KinematicVectorArray(lhs.rows(),lhs.VECTOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<lhs.size();k++){
        i_tmp = lhs.i_vec[k];
        j_tmp = lhs.j_vec[k];
        tmp[i_tmp] += lhs.buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialVectorArray operator* (MaterialTensorSparseMatrix& lhs, MaterialVectorArray& rhs){
    assert(lhs.cols() == rhs.size() && "Tensor sparse matrix multiplication failed.");
    MaterialVectorArray tmp = MaterialVectorArray(lhs.rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<lhs.size();k++){
        i_tmp = lhs.i_vec[k];
        j_tmp = lhs.j_vec[k];
        tmp[i_tmp] += lhs.buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialVectorArray operator* (MaterialTensorSparseMatrix& lhs, KinematicVectorArray& rhs){
    assert(lhs.cols() == rhs.size() && "Tensor sparse matrix multiplication failed.");
    MaterialVectorArray tmp = MaterialVectorArray(lhs.rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<lhs.size();k++){
        i_tmp = lhs.i_vec[k];
        j_tmp = lhs.j_vec[k];
        tmp[i_tmp] += lhs.buffer[k] * rhs[j_tmp];
    }
    return tmp;
}


MaterialVectorArray operator* (KinematicTensorSparseMatrix& lhs, MaterialVectorArray& rhs){
    assert(lhs.cols() == rhs.size() && "Tensor sparse matrix multiplication failed.");
    MaterialVectorArray tmp = MaterialVectorArray(lhs.rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<lhs.size();k++){
        i_tmp = lhs.i_vec[k];
        j_tmp = lhs.j_vec[k];
        tmp[i_tmp] += lhs.buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

KinematicVectorArray operator* (KinematicTensorSparseMatrix& lhs, KinematicVectorArray& rhs){
    assert(lhs.cols() == rhs.size() && lhs.TENSOR_TYPE == rhs.VECTOR_TYPE && "Tensor sparse matrix multiplication failed.");
    KinematicVectorArray tmp = KinematicVectorArray(lhs.rows(),rhs.VECTOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<lhs.size();k++){
        i_tmp = lhs.i_vec[k];
        j_tmp = lhs.j_vec[k];
        tmp[i_tmp] += lhs.buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialTensorArray operator* (MaterialTensorSparseMatrix& lhs, Eigen::VectorXd& rhs){
    assert(lhs.cols() == rhs.rows() && "Vector sparse matrix multiplication failed.");
    MaterialTensorArray tmp = MaterialTensorArray(lhs.rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<lhs.size();k++){
        i_tmp = lhs.i_vec[k];
        j_tmp = lhs.j_vec[k];
        tmp[i_tmp] += lhs.buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

KinematicTensorArray operator* (KinematicTensorSparseMatrix& lhs, Eigen::VectorXd& rhs){
    assert(lhs.cols() == rhs.rows() && "Vector sparse matrix multiplication failed.");
    KinematicTensorArray tmp = KinematicTensorArray(lhs.rows(),lhs.TENSOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<lhs.size();k++){
        i_tmp = lhs.i_vec[k];
        j_tmp = lhs.j_vec[k];
        tmp[i_tmp] += lhs.buffer[k] * rhs[j_tmp];
    }
    return tmp;
}


/*----------------------------------------------------------------------------*/

//diagonalize
Eigen::VectorXd MPMScalarSparseMatrix::rowsum(){
    Eigen::VectorXd tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k];
    }
    return tmp;
}

Eigen::VectorXd MPMScalarSparseMatrix::colsum(){
    Eigen::VectorXd tmp(cols());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[j_tmp] += buffer[k];
    }
    return tmp;
}


/*----------------------------------------------------------------------------*/

//diagonalize
MaterialVectorArray MaterialVectorSparseMatrix::rowsum(){
    MaterialVectorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k];
    }
    return tmp;
}

MaterialVectorArray MaterialVectorSparseMatrix::colsum(){
    MaterialVectorArray tmp(cols());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[j_tmp] += buffer[k];
    }
    return tmp;
}

//contraction operator
Eigen::VectorXd MaterialVectorSparseMatrix::dot(MaterialVectorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k].dot(other[j_tmp]);
    }
    return tmp;
}

Eigen::VectorXd MaterialVectorSparseMatrix::dot(KinematicVectorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k].dot(other[j_tmp]);
    }
    return tmp;
}

//tensor product operator
MaterialTensorArray MaterialVectorSparseMatrix::tensor(MaterialVectorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor product) failed.");
    MaterialTensorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k].tensor(other[j_tmp]);
    }
    return tmp;
}

MaterialTensorArray MaterialVectorSparseMatrix::tensor(KinematicVectorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor product) failed.");
    MaterialTensorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k].tensor(other[j_tmp]);
    }
    return tmp;
}

//transpose of tensor product
MaterialTensorArray MaterialVectorSparseMatrix::tensor_transpose(MaterialVectorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor product) failed.");
    MaterialTensorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].tensor(buffer[k]);
    }
    return tmp;
}

MaterialTensorArray MaterialVectorSparseMatrix::tensor_transpose(KinematicVectorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor product) failed.");
    MaterialTensorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].tensor(buffer[k]);
    }
    return tmp;
}

//left multiply by tensor
MaterialVectorArray MaterialVectorSparseMatrix::left_multiply(MaterialTensorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor multiplication) failed.");
    MaterialVectorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp]*buffer[k];
    }
    return tmp;
}

MaterialVectorArray MaterialVectorSparseMatrix::left_multiply(KinematicTensorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor multiplication) failed.");
    MaterialVectorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp]*buffer[k];
    }
    return tmp;
}

//left multiply by tensor transpose
MaterialVectorArray MaterialVectorSparseMatrix::left_multiply_by_transpose(MaterialTensorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor transpose multiplication) failed.");
    MaterialVectorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].transpose()*buffer[k];
    }
    return tmp;
}

MaterialVectorArray MaterialVectorSparseMatrix::left_multiply_by_transpose(KinematicTensorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor transpose multiplication) failed.");
    MaterialVectorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].transpose()*buffer[k];
    }
    return tmp;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

//diagonalize
KinematicVectorArray KinematicVectorSparseMatrix::rowsum(){
    KinematicVectorArray tmp(rows(),VECTOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k];
    }
    return tmp;
}

KinematicVectorArray KinematicVectorSparseMatrix::colsum(){
    KinematicVectorArray tmp(cols(),VECTOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[j_tmp] += buffer[k];
    }
    return tmp;
}

//contraction operator
Eigen::VectorXd KinematicVectorSparseMatrix::dot(MaterialVectorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k].dot(other[j_tmp]);
    }
    return tmp;
}

Eigen::VectorXd KinematicVectorSparseMatrix::dot(KinematicVectorArray& other){
    assert(other.size() == cols() && other.VECTOR_TYPE == VECTOR_TYPE && "Vector sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k].dot(other[j_tmp]);
    }
    return tmp;
}

//tensor product operator
MaterialTensorArray KinematicVectorSparseMatrix::tensor(MaterialVectorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor product) failed.");
    MaterialTensorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k].tensor(other[j_tmp]);
    }
    return tmp;
}

KinematicTensorArray KinematicVectorSparseMatrix::tensor(KinematicVectorArray& other){
    assert(other.size() == cols() && other.VECTOR_TYPE == VECTOR_TYPE && "Vector sparse matrix multiplication (tensor product) failed.");
    KinematicTensorArray tmp(rows(),VECTOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k].tensor(other[j_tmp]);
    }
    return tmp;
}

//transpose of tensor product
MaterialTensorArray KinematicVectorSparseMatrix::tensor_transpose(MaterialVectorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor transpose product) failed.");
    MaterialTensorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].tensor(buffer[k]);
    }
    return tmp;
}

KinematicTensorArray KinematicVectorSparseMatrix::tensor_transpose(KinematicVectorArray& other){
    assert(other.size() == cols() && other.VECTOR_TYPE == VECTOR_TYPE && "Vector sparse matrix multiplication (tensor transpose product) failed.");
    KinematicTensorArray tmp(rows(),VECTOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].tensor(buffer[k]);
    }
    return tmp;
}

//left multiply by tensor
MaterialVectorArray KinematicVectorSparseMatrix::left_multiply(MaterialTensorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor multiplication) failed.");
    MaterialVectorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp]*buffer[k];
    }
    return tmp;
}

KinematicVectorArray KinematicVectorSparseMatrix::left_multiply(KinematicTensorArray& other){
    assert(other.size() == cols() && other.TENSOR_TYPE == VECTOR_TYPE && "Vector sparse matrix multiplication (tensor multiplication) failed.");
    KinematicVectorArray tmp(rows(),VECTOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp]*buffer[k];
    }
    return tmp;
}

//left multiply by tensor transpose
MaterialVectorArray KinematicVectorSparseMatrix::left_multiply_by_transpose(MaterialTensorArray& other){
    assert(other.size() == cols() && "Vector sparse matrix multiplication (tensor transpose multiplication) failed.");
    MaterialVectorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].transpose()*buffer[k];
    }
    return tmp;
}

KinematicVectorArray KinematicVectorSparseMatrix::left_multiply_by_transpose(KinematicTensorArray& other){
    assert(other.size() == cols() && other.TENSOR_TYPE == VECTOR_TYPE && "Vector sparse matrix multiplication (tensor transpose multiplication) failed.");
    KinematicVectorArray tmp(rows(),VECTOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].transpose()*buffer[k];
    }
    return tmp;
}

/*----------------------------------------------------------------------------*/

//diagonalize
MaterialTensorArray MaterialTensorSparseMatrix::rowsum(){
    MaterialTensorArray tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k];
    }
    return tmp;
}

MaterialTensorArray MaterialTensorSparseMatrix::colsum(){
    MaterialTensorArray tmp(cols());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[j_tmp] += buffer[k];
    }
    return tmp;
}

//contraction operator
Eigen::VectorXd MaterialTensorSparseMatrix::dot(MaterialTensorArray& other){
    assert(other.size() == cols() && "Tensor sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].dot(buffer[k]);
    }
    return tmp;
}

Eigen::VectorXd MaterialTensorSparseMatrix::dot(KinematicTensorArray& other){
    assert(other.size() == cols() && "Tensor sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].dot(buffer[k]);
    }
    return tmp;
}

/*----------------------------------------------------------------------------*/

//diagonalize
KinematicTensorArray KinematicTensorSparseMatrix::rowsum(){
    KinematicTensorArray tmp(rows(),TENSOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += buffer[k];
    }
    return tmp;
}

KinematicTensorArray KinematicTensorSparseMatrix::colsum(){
    KinematicTensorArray tmp(cols(),TENSOR_TYPE);
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[j_tmp] += buffer[k];
    }
    return tmp;
}

//contraction operator
Eigen::VectorXd KinematicTensorSparseMatrix::dot(MaterialTensorArray& other){
    assert(other.size() == cols() && "Tensor sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].dot(buffer[k]);
    }
    return tmp;
}

Eigen::VectorXd KinematicTensorSparseMatrix::dot(KinematicTensorArray& other){
    assert(other.size() == cols() && other.TENSOR_TYPE == TENSOR_TYPE && "Tensor sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows());
    tmp.setZero();
    int i_tmp, j_tmp;
    for (int k=0;k<size();k++){
        i_tmp = i_vec[k];
        j_tmp = j_vec[k];
        tmp[i_tmp] += other[j_tmp].dot(buffer[k]);
    }
    return tmp;
}

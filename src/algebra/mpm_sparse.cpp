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
Eigen::VectorXd MPMScalarSparseMatrix::operate(const Eigen::VectorXd& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.rows() && "Scalar sparse matrix multiplication failed.");
    Eigen::VectorXd tmp = Eigen::VectorXd::Zero(rows(SPEC));

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialVectorArray MPMScalarSparseMatrix::operate(const MaterialVectorArray& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.size() && "Scalar sparse matrix multiplication failed.");
    MaterialVectorArray tmp = MaterialVectorArray(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

//??


KinematicVectorArray MPMScalarSparseMatrix::operate(const KinematicVectorArray& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.size() && "Scalar sparse matrix multiplication failed.");
    KinematicVectorArray tmp = KinematicVectorArray(rows(SPEC),rhs.VECTOR_TYPE);
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialTensorArray MPMScalarSparseMatrix::operate(const MaterialTensorArray& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.size() && "Scalar sparse matrix multiplication failed.");
    MaterialTensorArray tmp = MaterialTensorArray(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

KinematicTensorArray MPMScalarSparseMatrix::operate(const KinematicTensorArray& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.size() && "Scalar sparse matrix multiplication failed.");
    KinematicTensorArray tmp = KinematicTensorArray(rows(SPEC), rhs.TENSOR_TYPE);
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

Eigen::VectorXd operator* (const MPMScalarSparseMatrix& lhs, const Eigen::VectorXd& rhs){
    return lhs.operate(rhs);
    /*assert(lhs.cols() == rhs.rows() && "Scalar sparse matrix multiplication failed.");
    Eigen::VectorXd tmp = Eigen::VectorXd::Zero(lhs.rows());
    int i_tmp, j_tmp;
    for (int k=0;k<lhs.size();k++){
        i_tmp = lhs.i_vec[k];
        j_tmp = lhs.j_vec[k];
        tmp[i_tmp] += lhs.buffer[k] * rhs[j_tmp];
    }
    return tmp;*/
}

MaterialVectorArray operator* (const MPMScalarSparseMatrix& lhs, const MaterialVectorArray& rhs){
    return lhs.operate(rhs);
}

KinematicVectorArray operator* (const MPMScalarSparseMatrix& lhs, const KinematicVectorArray& rhs){
    return lhs.operate(rhs);
}

MaterialTensorArray operator* (const MPMScalarSparseMatrix& lhs, const MaterialTensorArray& rhs){
    return lhs.operate(rhs);
}
KinematicTensorArray operator* (const MPMScalarSparseMatrix& lhs, const KinematicTensorArray& rhs){
    return lhs.operate(rhs);
}



/*----------------------------------------------------------------------------*/
//vector sparse matrix multiplication

MaterialVectorArray MaterialVectorSparseMatrix::operate(const Eigen::VectorXd& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.rows() && "Vector sparse matrix multiplication failed.");
    MaterialVectorArray tmp = MaterialVectorArray(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialVectorArray operator* (const MaterialVectorSparseMatrix& lhs, const Eigen::VectorXd& rhs){
    return lhs.operate(rhs);
}

KinematicVectorArray KinematicVectorSparseMatrix::operate(const Eigen::VectorXd& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.rows() && "Vector sparse matrix multiplication failed.");
    KinematicVectorArray tmp = KinematicVectorArray(rows(SPEC), VECTOR_TYPE);
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

KinematicVectorArray operator* (const KinematicVectorSparseMatrix& lhs, const Eigen::VectorXd& rhs){
    return lhs.operate(rhs);
}



/*----------------------------------------------------------------------------*/
//tensor sparse matrix multiplication

MaterialVectorArray MaterialTensorSparseMatrix::operate(const MaterialVectorArray& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.size() && "Tensor sparse matrix multiplication failed.");
    MaterialVectorArray tmp = MaterialVectorArray(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialVectorArray operator* (const MaterialTensorSparseMatrix& lhs, const MaterialVectorArray& rhs){
    return lhs.operate(rhs);
}

MaterialVectorArray MaterialTensorSparseMatrix::operate(const KinematicVectorArray& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.size() && "Tensor sparse matrix multiplication failed.");
    MaterialVectorArray tmp = MaterialVectorArray(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialVectorArray operator* (const MaterialTensorSparseMatrix& lhs, const KinematicVectorArray& rhs){
    return lhs.operate(rhs);
}


MaterialVectorArray KinematicTensorSparseMatrix::operate(const MaterialVectorArray& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.size() && "Tensor sparse matrix multiplication failed.");
    MaterialVectorArray tmp = MaterialVectorArray(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialVectorArray operator* (const KinematicTensorSparseMatrix& lhs, const MaterialVectorArray& rhs){
    return lhs.operate(rhs);
}

KinematicVectorArray KinematicTensorSparseMatrix::operate(const KinematicVectorArray& rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.size() && "Tensor sparse matrix multiplication failed.");
    KinematicVectorArray tmp = KinematicVectorArray(rows(SPEC), rhs.VECTOR_TYPE);
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

KinematicVectorArray operator* (const KinematicTensorSparseMatrix& lhs, const KinematicVectorArray& rhs){
    return lhs.operate(rhs);
}

MaterialTensorArray MaterialTensorSparseMatrix::operate(const Eigen::VectorXd &rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.size() && "Tensor sparse matrix multiplication failed.");
    MaterialTensorArray tmp = MaterialTensorArray(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

MaterialTensorArray operator* (const MaterialTensorSparseMatrix& lhs, const Eigen::VectorXd& rhs){
    return lhs.operate(rhs);
}

KinematicTensorArray KinematicTensorSparseMatrix::operate(const Eigen::VectorXd &rhs, int SPEC) const {
    assert(cols(SPEC) == rhs.size() && "Tensor sparse matrix multiplication failed.");
    KinematicTensorArray tmp = KinematicTensorArray(rows(SPEC), TENSOR_TYPE);
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k] * rhs[j_tmp];
    }
    return tmp;
}

KinematicTensorArray operator* (const KinematicTensorSparseMatrix& lhs, const Eigen::VectorXd& rhs){
    return lhs.operate(rhs);
}


/*----------------------------------------------------------------------------*/

//diagonalize
Eigen::VectorXd MPMScalarSparseMatrix::rowsum() const {
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

Eigen::VectorXd MPMScalarSparseMatrix::colsum() const {
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
MaterialVectorArray MaterialVectorSparseMatrix::rowsum() const {
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

MaterialVectorArray MaterialVectorSparseMatrix::colsum() const {
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
Eigen::VectorXd MaterialVectorSparseMatrix::dot(const MaterialVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k].dot(other[j_tmp]);
    }
    return tmp;
}

Eigen::VectorXd MaterialVectorSparseMatrix::dot(const KinematicVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k].dot(other[j_tmp]);
    }
    return tmp;
}

//tensor product operator
MaterialTensorArray MaterialVectorSparseMatrix::tensor(const MaterialVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor product) failed.");
    MaterialTensorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k].tensor(other[j_tmp]);
    }
    return tmp;
}

MaterialTensorArray MaterialVectorSparseMatrix::tensor(const KinematicVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor product) failed.");
    MaterialTensorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k].tensor(other[j_tmp]);
    }
    return tmp;
}

//transpose of tensor product
MaterialTensorArray MaterialVectorSparseMatrix::tensor_product_transpose(const MaterialVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor product) failed.");
    MaterialTensorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].tensor(buffer[k]);
    }
    return tmp;
}

MaterialTensorArray MaterialVectorSparseMatrix::tensor_product_transpose(const KinematicVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor product) failed.");
    MaterialTensorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].tensor(buffer[k]);
    }
    return tmp;
}

//left multiply by tensor
MaterialVectorArray MaterialVectorSparseMatrix::left_multiply_by_tensor(const MaterialTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor multiplication) failed.");
    MaterialVectorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp]*buffer[k];
    }
    return tmp;
}

MaterialVectorArray MaterialVectorSparseMatrix::left_multiply_by_tensor(const KinematicTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor multiplication) failed.");
    MaterialVectorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp]*buffer[k];
    }
    return tmp;
}

//left multiply by tensor transpose
MaterialVectorArray MaterialVectorSparseMatrix::left_multiply_by_tensor_transpose(const MaterialTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor transpose multiplication) failed.");
    MaterialVectorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].transpose()*buffer[k];
    }
    return tmp;
}

MaterialVectorArray MaterialVectorSparseMatrix::left_multiply_by_tensor_transpose(const KinematicTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor transpose multiplication) failed.");
    MaterialVectorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].transpose()*buffer[k];
    }
    return tmp;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/

//diagonalize
KinematicVectorArray KinematicVectorSparseMatrix::rowsum() const {
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

KinematicVectorArray KinematicVectorSparseMatrix::colsum() const {
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
Eigen::VectorXd KinematicVectorSparseMatrix::dot(const MaterialVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k].dot(other[j_tmp]);
    }
    return tmp;
}

Eigen::VectorXd KinematicVectorSparseMatrix::dot(const KinematicVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && other.VECTOR_TYPE == VECTOR_TYPE && "Vector sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k].dot(other[j_tmp]);
    }
    return tmp;
}

//tensor product operator
MaterialTensorArray KinematicVectorSparseMatrix::tensor(const MaterialVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor product) failed.");
    MaterialTensorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k].tensor(other[j_tmp]);
    }
    return tmp;
}

KinematicTensorArray KinematicVectorSparseMatrix::tensor(const KinematicVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && other.VECTOR_TYPE == VECTOR_TYPE && "Vector sparse matrix multiplication (tensor product) failed.");
    KinematicTensorArray tmp(rows(SPEC),VECTOR_TYPE);
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += buffer[k].tensor(other[j_tmp]);
    }
    return tmp;
}

//transpose of tensor product
MaterialTensorArray KinematicVectorSparseMatrix::tensor_product_transpose(const MaterialVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor transpose product) failed.");
    MaterialTensorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].tensor(buffer[k]);
    }
    return tmp;
}

KinematicTensorArray KinematicVectorSparseMatrix::tensor_product_transpose(const KinematicVectorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && other.VECTOR_TYPE == VECTOR_TYPE && "Vector sparse matrix multiplication (tensor transpose product) failed.");
    KinematicTensorArray tmp(rows(SPEC),VECTOR_TYPE);
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].tensor(buffer[k]);
    }
    return tmp;
}

//left multiply by tensor
MaterialVectorArray KinematicVectorSparseMatrix::left_multiply_by_tensor(const MaterialTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor multiplication) failed.");
    MaterialVectorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp]*buffer[k];
    }
    return tmp;
}

KinematicVectorArray KinematicVectorSparseMatrix::left_multiply_by_tensor(const KinematicTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && other.TENSOR_TYPE == VECTOR_TYPE && "Vector sparse matrix multiplication (tensor multiplication) failed.");
    KinematicVectorArray tmp(rows(SPEC),VECTOR_TYPE);
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp]*buffer[k];
    }
    return tmp;
}

//left multiply by tensor transpose
MaterialVectorArray KinematicVectorSparseMatrix::left_multiply_by_tensor_transpose(const MaterialTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Vector sparse matrix multiplication (tensor transpose multiplication) failed.");
    MaterialVectorArray tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].transpose()*buffer[k];
    }
    return tmp;
}

KinematicVectorArray KinematicVectorSparseMatrix::left_multiply_by_tensor_transpose(const KinematicTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && other.TENSOR_TYPE == VECTOR_TYPE && "Vector sparse matrix multiplication (tensor transpose multiplication) failed.");
    KinematicVectorArray tmp(rows(SPEC),VECTOR_TYPE);
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].transpose()*buffer[k];
    }
    return tmp;
}

/*----------------------------------------------------------------------------*/

//diagonalize
MaterialTensorArray MaterialTensorSparseMatrix::rowsum() const {
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

MaterialTensorArray MaterialTensorSparseMatrix::colsum() const {
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
Eigen::VectorXd MaterialTensorSparseMatrix::dot(const MaterialTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Tensor sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].dot(buffer[k]);
    }
    return tmp;
}

Eigen::VectorXd MaterialTensorSparseMatrix::dot(const KinematicTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Tensor sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].dot(buffer[k]);
    }
    return tmp;
}

/*----------------------------------------------------------------------------*/

//diagonalize
KinematicTensorArray KinematicTensorSparseMatrix::rowsum() const {
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

KinematicTensorArray KinematicTensorSparseMatrix::colsum() const {
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
Eigen::VectorXd KinematicTensorSparseMatrix::dot(const MaterialTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && "Tensor sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].dot(buffer[k]);
    }
    return tmp;
}

Eigen::VectorXd KinematicTensorSparseMatrix::dot(const KinematicTensorArray& other, int SPEC) const {
    assert(other.size() == cols(SPEC) && other.TENSOR_TYPE == TENSOR_TYPE && "Tensor sparse matrix multiplication (dot product) failed.");
    Eigen::VectorXd tmp(rows(SPEC));
    tmp.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = get_j_index_ref(SPEC);

    for (int k=0;k<size();k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        tmp[i_tmp] += other[j_tmp].dot(buffer[k]);
    }
    return tmp;
}

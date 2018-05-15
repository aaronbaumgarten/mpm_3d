//
// Created by aaron on 5/12/18.
// mpm_sparse.hpp
//

#ifndef MPM_V3_MPM_SPARSE_HPP
#define MPM_V3_MPM_SPARSE_HPP

#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"

/*----------------------------------------------------------------------------*/
//sparse matrix for operating with Eigen::VectorXd
//currently only supports creation and operation (no 'access')
class MPMScalarSparseMatrix{
public:
protected:
    //sparse matrix info
    int i_max, j_max;

public:
    //sparse matrix data storage
    std::vector<int> i_vec; //list of row indices
    std::vector<int> j_vec; //list of column indices
    std::vector<double> buffer; //data stored at i,j position

    /*------------------------------------------------------------------------*/
    //initialize with max rows, cols
    MPMScalarSparseMatrix(int i, int j){
        i_max = i;
        j_max = j;
    }

    /*------------------------------------------------------------------------*/
    //transpose
    MPMScalarSparseMatrix transpose(){
        MPMScalarSparseMatrix tmp(j_max, i_max);
        tmp.i_vec = j_vec;
        tmp.j_vec = i_vec;
        tmp.buffer = buffer;
        return tmp;
    }

    /*------------------------------------------------------------------------*/
    //return properties of sparse matrix
    int rows() const { return i_max; }
    int cols() const { return j_max; }
    int size() const { return i_vec.size(); }

    //clear sparse matrix
    void clear(){
        i_vec.clear();
        j_vec.clear();
        buffer.clear();
        return;
    }

    //push_back a tuple
    void push_back(int& i, int& j, double& val){
        i_vec.push_back(i);
        j_vec.push_back(j);
        buffer.push_back(val);
    }

    //diagonalize
    Eigen::VectorXd rowsum();
    Eigen::VectorXd colsum();
};


/*----------------------------------------------------------------------------*/
//sparse matrix for operating with MaterialVectorArray
//currently only supports creation and operation (no 'access')
class MaterialVectorSparseMatrix{
protected:
    //sparse matrix info
    int i_max, j_max;

public:
    //sparse matrix data storage
    std::vector<int> i_vec; //list of row indices
    std::vector<int> j_vec; //list of column indices
    MaterialVectorArray buffer; //data stored at i,j position

    /*------------------------------------------------------------------------*/
    //initialize with max rows, cols
    MaterialVectorSparseMatrix(int i, int j){
        i_max = i;
        j_max = j;
    }

    /*------------------------------------------------------------------------*/
    //transpose
    MaterialVectorSparseMatrix transpose(){
        MaterialVectorSparseMatrix tmp(j_max, i_max);
        tmp.i_vec = j_vec;
        tmp.j_vec = i_vec;
        tmp.buffer = buffer;
        return tmp;
    }

    /*------------------------------------------------------------------------*/
    //return properties of sparse matrix
    int rows() const { return i_max; }
    int cols() const { return j_max; }
    int size() const { return i_vec.size(); }

    //clear sparse matrix
    void clear(){
        i_vec.clear();
        j_vec.clear();
        buffer.clear();
        return;
    }

    //push_back a tuple
    void push_back(int& i, int& j, MaterialVector& val){
        i_vec.push_back(i);
        j_vec.push_back(j);
        buffer.push_back(val);
    }

    //diagonalize
    MaterialVectorArray rowsum();
    MaterialVectorArray colsum();

    //contraction operator
    Eigen::VectorXd dot(MaterialVectorArray& other);
    Eigen::VectorXd dot(KinematicVectorArray& other);

    //tensor product operator
    MaterialTensorArray tensor(MaterialVectorArray& other);
    MaterialTensorArray tensor(KinematicVectorArray& other);

    //transpose of tensor product
    MaterialTensorArray tensor_transpose(MaterialVectorArray& other);
    MaterialTensorArray tensor_transpose(KinematicVectorArray& other);

    //left multiply by tensor
    MaterialVectorArray left_multiply(MaterialTensorArray& other);
    MaterialVectorArray left_multiply(KinematicTensorArray& other);

    //left multiply by tensor transpose
    MaterialVectorArray left_multiply_by_transpose(MaterialTensorArray& other);
    MaterialVectorArray left_multiply_by_transpose(KinematicTensorArray& other);
};


/*----------------------------------------------------------------------------*/
//sparse matrix for operating with KinematicVectorArray
//currently only supports creation and operation (no 'access')
class KinematicVectorSparseMatrix{
protected:
    //sparse matrix info
    int i_max, j_max;

public:
    //sparse matrix data storage
    std::vector<int> i_vec; //list of row indices
    std::vector<int> j_vec; //list of column indices
    KinematicVectorArray buffer; //data stored at i,j position

    //default dimensions for KinematicVectorArray
    int DIM = MPMVector::VECTOR_MAX_DIM;
    int VECTOR_TYPE = KinematicVector::VECTOR_3D;

    //assign KinematicVector type by input (rely on dummy KinematicTensor)
    void assignVectorType(int input){
        buffer.assignVectorType(input);
        DIM = buffer.DIM;
        VECTOR_TYPE = buffer.VECTOR_TYPE;
    }

    /*------------------------------------------------------------------------*/
    //initialize with max rows, cols
    KinematicVectorSparseMatrix(int i, int j){
        i_max = i;
        j_max = j;
    }

    //initialize with max rows, cols, and type
    KinematicVectorSparseMatrix(int i, int j, int input_type){
        i_max = i;
        j_max = j;
        assignVectorType(input_type);
    }

    /*------------------------------------------------------------------------*/
    //transpose
    KinematicVectorSparseMatrix transpose(){
        KinematicVectorSparseMatrix tmp(j_max, i_max, VECTOR_TYPE);
        tmp.i_vec = j_vec;
        tmp.j_vec = i_vec;
        tmp.buffer = buffer;
        return tmp;
    }

    /*------------------------------------------------------------------------*/
    //return properties of sparse matrix
    int rows() const { return i_max; }
    int cols() const { return j_max; }
    int size() const { return i_vec.size(); }

    //clear sparse matrix
    void clear(){
        i_vec.clear();
        j_vec.clear();
        buffer.clear();
        return;
    }

    //push_back a tuple
    void push_back(int& i, int& j, KinematicVector& val){
        i_vec.push_back(i);
        j_vec.push_back(j);
        buffer.push_back(val);
    }

    //diagonalize
    KinematicVectorArray rowsum();
    KinematicVectorArray colsum();

    //contraction operator
    Eigen::VectorXd dot(MaterialVectorArray& other);
    Eigen::VectorXd dot(KinematicVectorArray& other);

    //tensor product operator
    MaterialTensorArray tensor(MaterialVectorArray& other);
    KinematicTensorArray tensor(KinematicVectorArray& other);

    //transpose of tensor product
    MaterialTensorArray tensor_transpose(MaterialVectorArray& other);
    KinematicTensorArray tensor_transpose(KinematicVectorArray& other);

    //left multiply by tensor
    MaterialVectorArray left_multiply(MaterialTensorArray& other);
    KinematicVectorArray left_multiply(KinematicTensorArray& other);

    //left multiply by tensor transpose
    MaterialVectorArray left_multiply_by_transpose(MaterialTensorArray& other);
    KinematicVectorArray left_multiply_by_transpose(KinematicTensorArray& other);
};


/*----------------------------------------------------------------------------*/
//sparse matrix for operating with MaterialTensorArray
//currently only supports creation and operation (no 'access')
class MaterialTensorSparseMatrix{
protected:
    //sparse matrix info
    int i_max, j_max;

public:
    //sparse matrix data storage
    std::vector<int> i_vec; //list of row indices
    std::vector<int> j_vec; //list of column indices
    MaterialTensorArray buffer; //data stored at i,j position

    /*------------------------------------------------------------------------*/
    //initialize with max rows, cols
    MaterialTensorSparseMatrix(int i, int j){
        i_max = i;
        j_max = j;
    }

    /*------------------------------------------------------------------------*/
    //transpose
    MaterialTensorSparseMatrix transpose(){
        MaterialTensorSparseMatrix tmp(j_max, i_max);
        tmp.i_vec = j_vec;
        tmp.j_vec = i_vec;
        tmp.buffer = buffer;
        return tmp;
    }

    /*------------------------------------------------------------------------*/
    //return properties of sparse matrix
    int rows() const { return i_max; }
    int cols() const { return j_max; }
    int size() const { return i_vec.size(); }

    //clear sparse matrix
    void clear(){
        i_vec.clear();
        j_vec.clear();
        buffer.clear();
        return;
    }

    //push_back a tuple
    void push_back(int& i, int& j, MaterialTensor& val){
        i_vec.push_back(i);
        j_vec.push_back(j);
        buffer.push_back(val);
    }

    //diagonalize
    MaterialTensorArray rowsum();
    MaterialTensorArray colsum();

    //contraction operator
    Eigen::VectorXd dot(MaterialTensorArray& other);
    Eigen::VectorXd dot(KinematicTensorArray& other);
};


/*----------------------------------------------------------------------------*/
//sparse matrix for operating on KinematicTensorArray
//currently only supports creation and operation (no 'access')
class KinematicTensorSparseMatrix{
protected:
    //sparse matrix info
    int i_max, j_max;

public:
    //sparse matrix data storage
    std::vector<int> i_vec; //list of row indices
    std::vector<int> j_vec; //list of column indices
    KinematicTensorArray buffer; //data stored at i,j position

    //default dimensions for KinematicTensorArray
    int DIM = MPMTensor::TENSOR_MAX_DIM;
    int TENSOR_TYPE = KinematicTensor::TENSOR_3D;

    //assign KinematicTensor type by input
    void assignTensorType(int input){
        buffer.assignTensorType(input);
        DIM = buffer.DIM;
        TENSOR_TYPE = buffer.TENSOR_TYPE;
    }

    /*------------------------------------------------------------------------*/
    //initialize with max rows, cols
    KinematicTensorSparseMatrix(int i, int j){
        i_max = i;
        j_max = j;
    }

    //initialize with max rows, cols, tensor_type
    KinematicTensorSparseMatrix(int i, int j, int input_type){
        i_max = i;
        j_max = j;
        assignTensorType(input_type);
    }

    /*------------------------------------------------------------------------*/
    //transpose
    KinematicTensorSparseMatrix transpose(){
        KinematicTensorSparseMatrix tmp(j_max, i_max, TENSOR_TYPE);
        tmp.i_vec = j_vec;
        tmp.j_vec = i_vec;
        tmp.buffer = buffer;
        return tmp;
    }

    /*------------------------------------------------------------------------*/
    //return properties of sparse matrix
    int rows() const { return i_max; }
    int cols() const { return j_max; }
    int size() const { return i_vec.size(); }

    //clear sparse matrix
    void clear(){
        i_vec.clear();
        j_vec.clear();
        buffer.clear();
        return;
    }

    //push_back a tuple
    void push_back(int& i, int& j, KinematicTensor& val){
        i_vec.push_back(i);
        j_vec.push_back(j);
        buffer.push_back(val);
    }

    //diagonalize
    KinematicTensorArray rowsum();
    KinematicTensorArray colsum();

    //contraction operator
    Eigen::VectorXd dot(MaterialTensorArray& other);
    Eigen::VectorXd dot(KinematicTensorArray& other);
};


/*----------------------------------------------------------------------------*/
//sparse matrix multiplication operators
Eigen::VectorXd operator* (MPMScalarSparseMatrix& lhs, Eigen::VectorXd& rhs);

MaterialVectorArray operator* (MaterialVectorSparseMatrix& lhs, Eigen::VectorXd& rhs);
KinematicVectorArray operator* (KinematicVectorSparseMatrix& lhs, Eigen::VectorXd& rhs);

MaterialVectorArray operator* (MaterialTensorSparseMatrix& lhs, MaterialVectorArray& rhs);
MaterialVectorArray operator* (MaterialTensorSparseMatrix& lhs, KinematicVectorArray& rhs);
MaterialVectorArray operator* (KinematicTensorSparseMatrix& lhs, MaterialVectorArray& rhs);
KinematicVectorArray operator* (KinematicTensorSparseMatrix& lhs, KinematicVectorArray& rhs);

MaterialTensorArray operator* (MaterialTensorSparseMatrix& lhs, Eigen::VectorXd& rhs);
KinematicTensorArray operator* (KinematicTensorSparseMatrix& lhs, Eigen::VectorXd& rhs);

#endif //MPM_V3_MPM_SPARSE_HPP

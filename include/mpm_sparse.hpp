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

class MPMSparseMatrixBase{
protected:
    //sparse matrix info
    int i_max, j_max;

public:
    static const int NORMAL = 0;
    static const int TRANSPOSED = 1;

    //sparse matrix data storage
    std::vector<int> i_vec; //list of row indices
    std::vector<int> j_vec; //list of column indices

    /*------------------------------------------------------------------------*/
    //initialize with max rows, cols
    MPMSparseMatrixBase(int i, int j) {
        i_max = i;
        j_max = j;
    }

    //initialize without max rows, cols
    MPMSparseMatrixBase() {
        i_max = 0;
        j_max = 0;
        std::cerr << "WARNING: Initializing sparse matrix with zero size!" << std::endl;
    }

    /*------------------------------------------------------------------------*/
    //return properties of sparse matrix
    inline int rows(int SPEC = NORMAL) const {
        if (SPEC == NORMAL) {
            return i_max;
        } else if (SPEC == TRANSPOSED) {
            return j_max;
        } else {
            std::cerr << "Unknown SPEC for MPMSparseMatrixBase: " << SPEC << " Exiting." << std::endl;
            exit(0);
            return i_max;
        }
    }

    inline int cols(int SPEC = NORMAL) const {
        if (SPEC == NORMAL) {
            return j_max;
        } else if (SPEC == TRANSPOSED) {
            return i_max;
        } else {
            std::cerr << "Unknown SPEC for MPMSparseMatrixBase: " << SPEC << " Exiting." << std::endl;
            exit(0);
            return j_max;
        }
    }

    int size() const { return i_vec.size(); }

    /*------------------------------------------------------------------------*/
    //return index pair
    inline const std::vector<int>& get_i_index_ref(int SPEC = NORMAL) const {
        if (SPEC == NORMAL) {
            return i_vec;
        } else if (SPEC == TRANSPOSED){
            return j_vec;
        } else {
            std::cerr << "Unknown SPEC for MPMSparseMatrixBase: " << SPEC << " Exiting." << std::endl;
            exit(0);
            return i_vec;
        }
    }

    inline const std::vector<int>& get_j_index_ref(int SPEC = NORMAL) const {
        if (SPEC == NORMAL) {
            return j_vec;
        } else if (SPEC == TRANSPOSED){
            return i_vec;
        } else {
            std::cerr << "Unknown SPEC for MPMSparseMatrixBase: " << SPEC << " Exiting." << std::endl;
            exit(0);
            return j_vec;
        }
    }
};

/*----------------------------------------------------------------------------*/
//sparse matrix for operating with Eigen::VectorXd
//currently only supports creation and operation (no 'access')
class MPMScalarSparseMatrix : public MPMSparseMatrixBase{
public:
    //sparse matrix data storage
    std::vector<double> buffer; //data stored at i,j position

    /*------------------------------------------------------------------------*/
    //initialize with max rows, cols
    MPMScalarSparseMatrix(int i, int j) : MPMSparseMatrixBase(i,j) {}

    //initialize without max rows, cols
    MPMScalarSparseMatrix() : MPMSparseMatrixBase() {}

    /*------------------------------------------------------------------------*/
    //return transpose
    MPMScalarSparseMatrix transpose() const {
        MPMScalarSparseMatrix tmp(j_max, i_max);
        tmp.i_vec = j_vec;
        tmp.j_vec = i_vec;
        tmp.buffer = buffer;
        return tmp;
    }

    /*------------------------------------------------------------------------*/
    //clear sparse matrix
    void clear() {
        i_vec.clear();
        j_vec.clear();
        buffer.clear();
        return;
    }

    //resize the sparse matrix
    void resize(int size){
        //I hope you know what you are doing.
        i_vec.resize(size);
        j_vec.resize(size);
        buffer.resize(size);
    }

    //push_back a tuple
    void push_back(const int &i, const int &j, const double &val) {
        i_vec.push_back(i);
        j_vec.push_back(j);
        buffer.push_back(val);
    }

    //insert a tuple
    void insert(const int &k, const int &i, const int &j, const double &val) {
        if (k < i_vec.size()) {
            i_vec[k] = i;
            j_vec[k] = j;
            buffer[k] = val;
        } else {
            std::cerr << "ERROR: MPMScalarSparseMatrix insert() out of range! " << k << " > " << i_vec.size()-1 << "!";
        }
    }

    //diagonalize
    Eigen::VectorXd rowsum() const;

    Eigen::VectorXd colsum() const;

    Eigen::VectorXd operate(const Eigen::VectorXd& rhs, int SPEC = NORMAL) const;
    MaterialVectorArray operate(const MaterialVectorArray& rhs, int SPEC = NORMAL) const;
    KinematicVectorArray operate(const KinematicVectorArray& rhs, int SPEC = NORMAL) const;
    MaterialTensorArray operate(const MaterialTensorArray& rhs, int SPEC = NORMAL) const;
    KinematicTensorArray operate(const KinematicTensorArray& rhs, int SPEC = NORMAL) const;
};


/*----------------------------------------------------------------------------*/
//sparse matrix for operating with MaterialVectorArray
//currently only supports creation and operation (no 'access')
class MaterialVectorSparseMatrix : public MPMSparseMatrixBase{
public:
    MaterialVectorArray buffer; //data stored at i,j position

    /*------------------------------------------------------------------------*/
    //initialize with max rows, cols
    MaterialVectorSparseMatrix(int i, int j) : MPMSparseMatrixBase(i,j) {}

    //initialize without max rows, cols
    MaterialVectorSparseMatrix() : MPMSparseMatrixBase() {}

    /*------------------------------------------------------------------------*/
    //transpose
    MaterialVectorSparseMatrix transpose() const {
        MaterialVectorSparseMatrix tmp(j_max, i_max);
        tmp.i_vec = j_vec;
        tmp.j_vec = i_vec;
        tmp.buffer = buffer;
        return tmp;
    }

    /*------------------------------------------------------------------------*/
    //clear sparse matrix
    void clear() {
        i_vec.clear();
        j_vec.clear();
        buffer.clear();
        return;
    }

    //resize the sparse matrix
    void resize(int size){
        //I hope you know what you are doing.
        i_vec.resize(size);
        j_vec.resize(size);
        buffer.resize(size);
    }

    //push_back a tuple
    void push_back(const int &i, const int &j, const MaterialVector &val) {
        i_vec.push_back(i);
        j_vec.push_back(j);
        buffer.push_back(val);
    }

    //insert a tuple
    void insert(const int &k, const int &i, const int &j, const MaterialVector &val) {
        if (k < i_vec.size()) {
            i_vec[k] = i;
            j_vec[k] = j;
            buffer[k] = val;
        } else {
            std::cerr << "ERROR: MaterialVectorSparseMatrix insert() out of range! " << k << " > " << i_vec.size()-1 << "!";
        }
    }

    //diagonalize
    MaterialVectorArray rowsum() const;

    MaterialVectorArray colsum() const;

    //multiplication base
    MaterialVectorArray operate(const Eigen::VectorXd& rhs, int SPEC = NORMAL) const;

    //contraction operator
    Eigen::VectorXd dot(const MaterialVectorArray &other, int SPEC = NORMAL) const;

    Eigen::VectorXd dot(const KinematicVectorArray &other, int SPEC = NORMAL) const;

    //tensor product operator
    MaterialTensorArray tensor(const MaterialVectorArray &other, int SPEC = NORMAL) const;

    MaterialTensorArray tensor(const KinematicVectorArray &other, int SPEC = NORMAL) const;

    //transpose of tensor product
    MaterialTensorArray tensor_product_transpose(const MaterialVectorArray &other, int SPEC = NORMAL) const;

    MaterialTensorArray tensor_product_transpose(const KinematicVectorArray &other, int SPEC = NORMAL) const;

    //left multiply by tensor
    MaterialVectorArray left_multiply_by_tensor(const MaterialTensorArray &other, int SPEC = NORMAL) const;

    MaterialVectorArray left_multiply_by_tensor(const KinematicTensorArray &other, int SPEC = NORMAL) const;

    //left multiply by tensor transpose
    MaterialVectorArray left_multiply_by_tensor_transpose(const MaterialTensorArray &other, int SPEC = NORMAL) const;

    MaterialVectorArray left_multiply_by_tensor_transpose(const KinematicTensorArray &other, int SPEC = NORMAL) const;
};


/*----------------------------------------------------------------------------*/
//sparse matrix for operating with KinematicVectorArray
//currently only supports creation and operation (no 'access')
class KinematicVectorSparseMatrix : public MPMSparseMatrixBase{
public:
    //sparse matrix data storage
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
    KinematicVectorSparseMatrix(int i, int j) : MPMSparseMatrixBase(i,j) {}

    //initialize with max rows, cols, and type
    KinematicVectorSparseMatrix(int i, int j, int input_type) : MPMSparseMatrixBase(i,j) {
        assignVectorType(input_type);
    }

    //initiialize without max rows, cols
    KinematicVectorSparseMatrix() : MPMSparseMatrixBase() {}

    /*------------------------------------------------------------------------*/
    //transpose
    KinematicVectorSparseMatrix transpose() const {
        KinematicVectorSparseMatrix tmp(j_max, i_max, VECTOR_TYPE);
        tmp.i_vec = j_vec;
        tmp.j_vec = i_vec;
        tmp.buffer = buffer;
        return tmp;
    }

    /*------------------------------------------------------------------------*/

    //clear sparse matrix
    void clear(){
        i_vec.clear();
        j_vec.clear();
        buffer.clear();
        return;
    }

    //resize the sparse matrix
    void resize(int size){
        //I hope you know what you are doing.
        //std::cout << size << std::endl;
        //exit(0);
        i_vec.resize(size);
        j_vec.resize(size);
        buffer.resize(size);
    }

    //push_back a tuple
    void push_back(const int& i, const int& j, const KinematicVector& val){
        i_vec.push_back(i);
        j_vec.push_back(j);
        buffer.push_back(val);
    }

    //insert a tuple
    void insert(const int &k, const int &i, const int &j, const KinematicVector& val) {
        if (k < i_vec.size()) {
            i_vec[k] = i;
            j_vec[k] = j;
            buffer[k] = val;
        } else {
            std::cerr << "ERROR: KinematicVectorSparseMatrix insert() out of range! " << k << " > " << i_vec.size()-1 << "!";
        }
    }

    //diagonalize
    KinematicVectorArray rowsum() const;
    KinematicVectorArray colsum() const;

    KinematicVectorArray operate(const Eigen::VectorXd& rhs, int SPEC = NORMAL) const;

    //contraction operator
    Eigen::VectorXd dot(const MaterialVectorArray& other, int SPEC = NORMAL) const;
    Eigen::VectorXd dot(const KinematicVectorArray& other, int SPEC = NORMAL) const;

    //tensor product operator
    MaterialTensorArray tensor(const MaterialVectorArray& other, int SPEC = NORMAL) const;
    KinematicTensorArray tensor(const KinematicVectorArray& other, int SPEC = NORMAL) const;

    //transpose of tensor product
    MaterialTensorArray tensor_product_transpose(const MaterialVectorArray& other, int SPEC = NORMAL) const;
    KinematicTensorArray tensor_product_transpose(const KinematicVectorArray& other, int SPEC = NORMAL) const;

    //left multiply by tensor
    MaterialVectorArray left_multiply_by_tensor(const MaterialTensorArray& other, int SPEC = NORMAL) const;
    KinematicVectorArray left_multiply_by_tensor(const KinematicTensorArray& other, int SPEC = NORMAL) const;

    //left multiply by tensor transpose
    MaterialVectorArray left_multiply_by_tensor_transpose(const MaterialTensorArray& other, int SPEC = NORMAL) const;
    KinematicVectorArray left_multiply_by_tensor_transpose(const KinematicTensorArray& other, int SPEC = NORMAL) const;
};


/*----------------------------------------------------------------------------*/
//sparse matrix for operating with MaterialTensorArray
//currently only supports creation and operation (no 'access')
class MaterialTensorSparseMatrix : public MPMSparseMatrixBase{
public:
    //sparse matrix data storage
    MaterialTensorArray buffer; //data stored at i,j position

    /*------------------------------------------------------------------------*/
    //initialize with max rows, cols
    MaterialTensorSparseMatrix(int i, int j) : MPMSparseMatrixBase(i,j) {}

    /*------------------------------------------------------------------------*/
    //transpose
    MaterialTensorSparseMatrix transpose() const {
        MaterialTensorSparseMatrix tmp(j_max, i_max);
        tmp.i_vec = j_vec;
        tmp.j_vec = i_vec;
        tmp.buffer = buffer;
        return tmp;
    }

    /*------------------------------------------------------------------------*/

    //clear sparse matrix
    void clear(){
        i_vec.clear();
        j_vec.clear();
        buffer.clear();
        return;
    }

    //push_back a tuple
    void push_back(const int& i, const int& j, const MaterialTensor& val){
        i_vec.push_back(i);
        j_vec.push_back(j);
        buffer.push_back(val);
    }

    //diagonalize
    MaterialTensorArray rowsum() const;
    MaterialTensorArray colsum() const;

    //contraction operator
    Eigen::VectorXd dot(const MaterialTensorArray& other, int SPEC = NORMAL) const;
    Eigen::VectorXd dot(const KinematicTensorArray& other, int SPEC = NORMAL) const;

    MaterialVectorArray operate(const MaterialVectorArray& rhs, int SPEC = NORMAL) const;
    MaterialVectorArray operate(const KinematicVectorArray& rhs, int SPEC = NORMAL) const;

    MaterialTensorArray operate(const Eigen::VectorXd& rhs, int SPEC = NORMAL) const;
};


/*----------------------------------------------------------------------------*/
//sparse matrix for operating on KinematicTensorArray
//currently only supports creation and operation (no 'access')
class KinematicTensorSparseMatrix : public MPMSparseMatrixBase{
public:
    //sparse matrix data storage
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
    KinematicTensorSparseMatrix(int i, int j) : MPMSparseMatrixBase(i,j) {}

    //initialize with max rows, cols, tensor_type
    KinematicTensorSparseMatrix(int i, int j, int input_type) : MPMSparseMatrixBase(i,j){
        assignTensorType(input_type);
    }

    /*------------------------------------------------------------------------*/
    //transpose
    KinematicTensorSparseMatrix transpose() const {
        KinematicTensorSparseMatrix tmp(j_max, i_max, TENSOR_TYPE);
        tmp.i_vec = j_vec;
        tmp.j_vec = i_vec;
        tmp.buffer = buffer;
        return tmp;
    }

    /*------------------------------------------------------------------------*/

    //clear sparse matrix
    void clear(){
        i_vec.clear();
        j_vec.clear();
        buffer.clear();
        return;
    }

    //push_back a tuple
    void push_back(const int& i, const int& j, const KinematicTensor& val){
        i_vec.push_back(i);
        j_vec.push_back(j);
        buffer.push_back(val);
    }

    //diagonalize
    KinematicTensorArray rowsum() const;
    KinematicTensorArray colsum() const;

    //contraction operator
    Eigen::VectorXd dot(const MaterialTensorArray& other, int SPEC = NORMAL) const;
    Eigen::VectorXd dot(const KinematicTensorArray& other, int SPEC = NORMAL) const;

    MaterialVectorArray operate(const MaterialVectorArray& rhs, int SPEC = NORMAL) const;
    KinematicVectorArray operate(const KinematicVectorArray& rhs, int SPEC = NORMAL) const;

    KinematicTensorArray operate(const Eigen::VectorXd& rhs, int SPEC = NORMAL) const;
};


/*----------------------------------------------------------------------------*/
//sparse matrix multiplication operators
Eigen::VectorXd operator* (const MPMScalarSparseMatrix& lhs, const Eigen::VectorXd& rhs);
MaterialVectorArray operator* (const MPMScalarSparseMatrix& lhs, const MaterialVectorArray& rhs);
KinematicVectorArray operator* (const MPMScalarSparseMatrix& lhs, const KinematicVectorArray& rhs);
MaterialTensorArray operator* (const MPMScalarSparseMatrix& lhs, const MaterialTensorArray& rhs);
KinematicTensorArray operator* (const MPMScalarSparseMatrix& lhs, const KinematicTensorArray& rhs);

MaterialVectorArray operator* (const MaterialVectorSparseMatrix& lhs, const Eigen::VectorXd& rhs);
KinematicVectorArray operator* (const KinematicVectorSparseMatrix& lhs, const Eigen::VectorXd& rhs);

MaterialVectorArray operator* (const MaterialTensorSparseMatrix& lhs, const MaterialVectorArray& rhs);
MaterialVectorArray operator* (const MaterialTensorSparseMatrix& lhs, const KinematicVectorArray& rhs);
MaterialVectorArray operator* (const KinematicTensorSparseMatrix& lhs, const MaterialVectorArray& rhs);
KinematicVectorArray operator* (const KinematicTensorSparseMatrix& lhs, const KinematicVectorArray& rhs);

MaterialTensorArray operator* (const MaterialTensorSparseMatrix& lhs, const Eigen::VectorXd& rhs);
KinematicTensorArray operator* (const KinematicTensorSparseMatrix& lhs, const Eigen::VectorXd& rhs);

#endif //MPM_V3_MPM_SPARSE_HPP

//
// Created by aaron on 5/31/19.
// parallel_explicit_usl.hpp
//

#ifndef MPM_V3_PARALLEL_EXPLICIT_USL_HPP
#define MPM_V3_PARALLEL_EXPLICIT_USL_HPP

#include "mpm_objects.hpp"
#include "mpm_vector.hpp"
#include "mpm_tensor.hpp"
#include "mpm_vectorarray.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_sparse.hpp"
#include "job.hpp"

#include "solvers.hpp"
#include "objects/bodies/bodies.hpp"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <thread>
#include <time.h>


//scalar sparse matrix multiply scalar vector
void ParallelExplicitUSL::ssmOperateStoS(const MPMScalarSparseMatrix &S,
                                         const Eigen::VectorXd &x,
                                         Eigen::VectorXd &lhs,
                                         int SPEC,
                                         int k_begin, int k_end){
    assert(S.cols(SPEC) == x.rows() && "Scalar sparse matrix multiplication failed.");

    //zero lhs output
    lhs.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = S.get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = S.get_j_index_ref(SPEC);

    for (int k=k_begin;k<=k_end;k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        lhs[i_tmp] += S.buffer[k] * x[j_tmp];
    }

    return;
}

//scalar sparse matrix multiply vector vector
void ParallelExplicitUSL::ssmOperateVtoV(const MPMScalarSparseMatrix &S,
                                         const KinematicVectorArray &x,
                                         KinematicVectorArray &lhs,
                                         int SPEC,
                                         int k_begin, int k_end){
    assert(S.cols(SPEC) == x.size() && "Scalar sparse matrix multiplication failed.");
    lhs.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = S.get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = S.get_j_index_ref(SPEC);

    for (int k=k_begin;k<=k_end;k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        lhs[i_tmp] += S.buffer[k] * x[j_tmp];
    }

    return;
}

//kinematic vector sparse matrix left multiply by tensor vector
void ParallelExplicitUSL::kvsmLeftMultiply(const KinematicVectorSparseMatrix &gradS,
                                           const MaterialTensorArray &T,
                                           MaterialVectorArray &lhs,
                                           int SPEC,
                                           int k_begin, int k_end){
    assert(T.size() == gradS.cols(SPEC) && "Vector sparse matrix multiplication (tensor multiplication) failed.");
    lhs.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = gradS.get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = gradS.get_j_index_ref(SPEC);

    for (int k=k_begin;k<=k_end;k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        lhs[i_tmp] += T[j_tmp]*gradS.buffer[k];
    }

    return;
}

//kinematic vector sparse matrix tensor product with vector vector
void ParallelExplicitUSL::kvsmTensorProductT(const KinematicVectorSparseMatrix &gradS,
                                             const KinematicVectorArray &x,
                                             KinematicTensorArray &L,
                                             int SPEC,
                                             int k_begin, int k_end){
    assert(x.size() == gradS.cols(SPEC) && x.VECTOR_TYPE == gradS.VECTOR_TYPE && "Vector sparse matrix multiplication (tensor transpose product) failed.");
    L.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = gradS.get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = gradS.get_j_index_ref(SPEC);

    for (int k=k_begin;k<=k_end;k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        L[i_tmp] += x[j_tmp].tensor(gradS.buffer[k]);
    }

    return;
}

//kinematic vecotr sparse matrix multiplying scalar vector
void ParallelExplicitUSL::kvsmOperateStoV(const KinematicVectorSparseMatrix &gradS,
                                             const Eigen::VectorXd &x,
                                             KinematicVectorArray &lhs,
                                             int SPEC,
                                             int k_begin, int k_end){
    assert(x.rows() == gradS.cols(SPEC) && lhs.VECTOR_TYPE == gradS.VECTOR_TYPE && "Vector sparse matrix multiplication (matrix multiplication) failed.");
    lhs.setZero();

    int i_tmp, j_tmp;
    const std::vector<int>& i_ref = gradS.get_i_index_ref(SPEC);
    const std::vector<int>& j_ref = gradS.get_j_index_ref(SPEC);

    for (int k=k_begin;k<=k_end;k++){
        i_tmp = i_ref[k];
        j_tmp = j_ref[k];
        lhs[i_tmp] += gradS.buffer[k] * x(j_tmp);
    }

    return;
}

/*----------------------------------------------------------------------------*/
//
//method for adding list of scalar vectors
void ParallelExplicitUSL::scalarAdd(const std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& list,
                                    Eigen::VectorXd &sum,
                                    int i_begin, int i_end,
                                    bool clear){
    //fill in i_begin to i_end components of sum
    assert(sum.size() == list[0].size() && "scalarAdd failed.");
    //set to zero
    if (clear) {
        for (int i = i_begin; i <= i_end; i++) {
            sum(i) = 0;
        }
    }
    //add components
    for (int l=0; l<list.size(); l++){
        for (int i=i_begin; i<=i_end; i++){
            sum(i) += list[l](i);
        }
    }
    return;
}

//method for adding list of kinematic vector vectors
void ParallelExplicitUSL::vectorAddK(const std::vector<KinematicVectorArray>& list,
                                     KinematicVectorArray &sum,
                                     int i_begin, int i_end,
                                     bool clear){
    //fill in i_begin to i_end components of sum
    assert(sum.size() == list[0].size() && "vectorAdd failed.");
    //set to zero
    if (clear) {
        for (int i = i_begin; i <= i_end; i++) {
            sum(i).setZero();
        }
    }
    //add components
    for (int l=0; l<list.size(); l++){
        for (int i=i_begin; i<=i_end; i++){
            sum(i) += list[l](i);
        }
    }
    return;
}

//method for adding list of material vector vectors
void ParallelExplicitUSL::vectorAddM(const std::vector<MaterialVectorArray>& list,
                                     MaterialVectorArray &sum,
                                     int i_begin, int i_end,
                                     bool clear){
    //fill in i_begin to i_end components of sum
    assert(sum.size() == list[0].size() && "vectorAdd failed.");
    //set to zero
    if (clear) {
        for (int i = i_begin; i <= i_end; i++) {
            sum(i).setZero();
        }
    }
    //add components
    for (int l=0; l<list.size(); l++){
        for (int i=i_begin; i<=i_end; i++){
            sum(i) += list[l](i);
        }
    }
    return;
}

//method for adding list of kinematic tensor vectors
void ParallelExplicitUSL::tensorAdd(const std::vector<KinematicTensorArray>& list,
                                    KinematicTensorArray &sum,
                                    int i_begin, int i_end,
                                    bool clear){
    //fill in i_begin to i_end components of sum
    assert(sum.size() == list[0].size() && "tensorAdd failed.");
    //set to zero
    if (clear) {
        for (int i = i_begin; i <= i_end; i++) {
            sum(i).setZero();
        }
    }
    //add components
    for (int l=0; l<list.size(); l++){
        for (int i=i_begin; i<=i_end; i++){
            sum(i) += list[l](i);
        }
    }
    return;
}

/*----------------------------------------------------------------------------*/
void ParallelExplicitUSL::parallelMultiply(const MPMScalarSparseMatrix &S,
                                           const Eigen::VectorXd &x,
                                           Eigen::VectorXd &lhs,
                                           int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = S.size() - 1;

    //get length of output vector
    int i_max = lhs.rows() - 1;

    //determine number of threads for matrix vector mult
    int thread_count = 0;
    if (S.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = S.size();
    }

    //intermediate storage vector
    std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>> lhs_vec(thread_count);

    //choose interval size
    int k_interval = (S.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        lhs_vec[t] = Eigen::VectorXd(lhs.rows());

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        threads[t] = std::thread(ssmOperateStoS, std::ref(S), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    //determine number of threads for addition
    if (lhs.rows() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.rows();
    }

    //choose interval size
    int i_interval = (lhs.rows()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(scalarAdd, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    //should be all done
    return;
}

void ParallelExplicitUSL::parallelMultiply(const MPMScalarSparseMatrix &S,
                                           const KinematicVectorArray &x,
                                           KinematicVectorArray &lhs,
                                           int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = S.size() - 1;

    //get length of output vector
    int i_max = lhs.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count = 0;
    if (S.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = S.size();
    }

    //intermediate storage vector
    std::vector<KinematicVectorArray> lhs_vec(thread_count);

    //choose interval size
    int k_interval = (S.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        lhs_vec[t] = KinematicVectorArray(lhs.size(),lhs.VECTOR_TYPE);

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        threads[t] = std::thread(ssmOperateVtoV, std::ref(S), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    //determine number of threads for addition
    if (lhs.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.size();
    }

    //choose interval size
    int i_interval = (lhs.size()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(vectorAddK, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    //should be all done
    return;
}

void ParallelExplicitUSL::parallelMultiply(const KinematicVectorSparseMatrix &gradS,
                                           const MaterialTensorArray &T,
                                           MaterialVectorArray &lhs,
                                           int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = gradS.size() - 1;

    //get length of output vector
    int i_max = lhs.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count = 0;
    if (gradS.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = gradS.size();
    }

    //intermediate storage vector
    std::vector<MaterialVectorArray> lhs_vec(thread_count);

    //choose interval size
    int k_interval = (gradS.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        lhs_vec[t] = MaterialVectorArray(lhs.size());

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        threads[t] = std::thread(kvsmLeftMultiply, std::ref(gradS), std::ref(T), std::ref(lhs_vec[t]), SPEC, k_begin, k_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    //determine number of threads for addition
    if (lhs.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.size();
    }

    //choose interval size
    int i_interval = (lhs.size()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(vectorAddM, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    //should be all done
    return;
}


void ParallelExplicitUSL::parallelMultiply(const KinematicVectorSparseMatrix &gradS,
                                           const KinematicVectorArray &x,
                                           KinematicTensorArray &L,
                                           int SPEC, bool clear, int memUnitID){
    //get length of sparse matrix storage
    int k_max = gradS.size() - 1;

    //get length of output vector
    int i_max = L.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count = 0;
    if (gradS.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = gradS.size();
    }

    //intermediate storage vector
    std::vector<KinematicTensorArray> lhs_vec(thread_count);

    //choose interval size
    int k_interval = (gradS.size()/thread_count) + 1;
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        lhs_vec[t] = KinematicTensorArray(L.size(),L.TENSOR_TYPE);

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        threads[t] = std::thread(kvsmTensorProductT, std::ref(gradS), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    //determine number of threads for addition
    if (L.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = L.size();
    }

    //choose interval size
    int i_interval = (L.size()/thread_count) + 1;
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(tensorAdd, std::ref(lhs_vec), std::ref(L), i_begin, i_end, clear);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    //should be all done
    return;
}


/*----------------------------------------------------------------------------*/
// methods for quick sums, divisions, and multiplications of vector objects
void ParallelExplicitUSL::multiplyKVbyS(KinematicVectorArray &kv,
                          Eigen::VectorXd &s,
                          double scale,
                          KinematicVectorArray &out,
                          int i_begin, int i_end){
    for (int i=i_begin; i<=i_end; i++){
        out(i) += scale * kv(i) * s(i);
    }
    return;
}

void ParallelExplicitUSL::divideKVbyS(KinematicVectorArray &kv,
                        Eigen::VectorXd &s,
                        double scale,
                        KinematicVectorArray &out,
                        int i_begin, int i_end){
    for (int i=i_begin; i<=i_end; i++){
        if (s(i) > 0){
            out(i) += scale * kv(i) / s(i);
        }
    }
    return;
}

void ParallelExplicitUSL::multiplySbyS(Eigen::VectorXd &a,
                         Eigen::VectorXd &b,
                         double scale,
                         Eigen::VectorXd &out,
                         int i_begin, int i_end){
    for (int i=i_begin; i<=i_end; i++){
        out(i) += scale * a(i) * b(i);
    }
    return;
}

void ParallelExplicitUSL::multiplyMTbyS(MaterialTensorArray &mt,
                          Eigen::VectorXd &s,
                          double scale,
                          MaterialTensorArray &out,
                          int i_begin, int i_end){
    for (int i=i_begin; i<=i_end; i++){
        out(i) += scale * mt(i) * s(i);
    }
    return;
}

void ParallelExplicitUSL::addStoS(Eigen::VectorXd &a,
                    Eigen::VectorXd &b,
                    double scale,
                    Eigen::VectorXd &out,
                    int i_begin, int i_end){
    for (int i=i_begin; i<=i_end; i++){
        out(i) = a(i) + scale*b(i);
    }
    return;
}

void ParallelExplicitUSL::subtractSfromS(Eigen::VectorXd &a,
                           Eigen::VectorXd &b,
                           double scale,
                           Eigen::VectorXd &out,
                           int i_begin, int i_end){
    for (int i=i_begin; i<=i_end; i++){
        out(i) += a(i) - scale*b(i);
    }
    return;
}

void ParallelExplicitUSL::zeroKV(KinematicVectorArray &kv, int i_begin, int i_end){
    for (int i=i_begin; i<=i_end; i++){
        kv(i).setZero();
    }
    return;
}

void ParallelExplicitUSL::zeroS(Eigen::VectorXd &s, int i_begin, int i_end){
    for (int i=i_begin; i<=i_end; i++){
        s(i) = 0.;
    }
    return;
}

void ParallelExplicitUSL::zeroMT(MaterialTensorArray &mt, int i_begin, int i_end){
    for (int i=i_begin; i<=i_end; i++){
        mt(i).setZero();
    }
    return;
}

/*----------------------------------------------------------------------------*/
//methods to manage threads during parallel operations

void ParallelExplicitUSL::parallelMultiply(KinematicVectorArray &kv,
                                                  Eigen::VectorXd &s,
                                                  double scale,
                                                  KinematicVectorArray &out,
                                                  bool clear){
    //determine maximum index
    int thread_count = 0;
    int i_max = kv.size() - 1;

    //determine number of threads for addition
    if (kv.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = kv.size();
    }

    //choose interval size
    int i_interval = (kv.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            threads[t] = std::thread(zeroKV, std::ref(out), i_begin, i_end);
        }

        //join threads
        for (int t=0; t<thread_count; t++){
            threads[t].join();
        }
    }

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(multiplyKVbyS, std::ref(kv), std::ref(s), scale, std::ref(out), i_begin, i_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    return;
}

void ParallelExplicitUSL::parallelDivide(KinematicVectorArray &kv,
                                                Eigen::VectorXd &s,
                                                double scale,
                                                KinematicVectorArray &out,
                                                bool clear){
    //determine maximum index
    int thread_count = 0;
    int i_max = kv.size() - 1;

    //determine number of threads for addition
    if (kv.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = kv.size();
    }

    //choose interval size
    int i_interval = (kv.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            threads[t] = std::thread(zeroKV, std::ref(out), i_begin, i_end);
        }

        //join threads
        for (int t=0; t<thread_count; t++){
            threads[t].join();
        }
    }

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(divideKVbyS, std::ref(kv), std::ref(s), scale, std::ref(out), i_begin, i_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    return;
}

void ParallelExplicitUSL::parallelMultiply(Eigen::VectorXd &a,
                                           Eigen::VectorXd &b,
                                           double scale,
                                           Eigen::VectorXd &out,
                                           bool clear){
    //determine maximum index
    int thread_count = 0;
    int i_max = a.size() - 1;

    //determine number of threads for addition
    if (a.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = a.size();
    }

    //choose interval size
    int i_interval = (a.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            threads[t] = std::thread(zeroS, std::ref(out), i_begin, i_end);
        }

        //join threads
        for (int t=0; t<thread_count; t++){
            threads[t].join();
        }
    }

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(multiplySbyS, std::ref(a), std::ref(b), scale, std::ref(out), i_begin, i_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    return;
}

void ParallelExplicitUSL::parallelMultiply(MaterialTensorArray &mt,
                                           Eigen::VectorXd &s,
                                           double scale,
                                           MaterialTensorArray &out,
                                           bool clear){
    //determine maximum index
    int thread_count = 0;
    int i_max = mt.size() - 1;

    //determine number of threads for addition
    if (mt.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = mt.size();
    }

    //choose interval size
    int i_interval = (mt.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            threads[t] = std::thread(zeroMT, std::ref(out), i_begin, i_end);
        }

        //join threads
        for (int t=0; t<thread_count; t++){
            threads[t].join();
        }
    }

    //multipy MaterialTensor by Scalar and add to output
    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(multiplyMTbyS, std::ref(mt), std::ref(s), scale, std::ref(out), i_begin, i_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    return;
}

void ParallelExplicitUSL::parallelAdd(Eigen::VectorXd &a,
                                      Eigen::VectorXd &b,
                                      double scale,
                                      Eigen::VectorXd &out,
                                      bool clear){

    //determine maximum index
    int thread_count = 0;
    int i_max = a.size() - 1;

    //determine number of threads for addition
    if (a.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = a.size();
    }

    //choose interval size
    int i_interval = (a.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            threads[t] = std::thread(zeroS, std::ref(out), i_begin, i_end);
        }

        //join threads
        for (int t=0; t<thread_count; t++){
            threads[t].join();
        }
    }

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(addStoS, std::ref(a), std::ref(b), scale, std::ref(out), i_begin, i_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    return;
}

void ParallelExplicitUSL::parallelSubtract(Eigen::VectorXd &a,
                                           Eigen::VectorXd &b,
                                           double scale,
                                           Eigen::VectorXd &out,
                                           bool clear){
    //determine maximum index
    int thread_count = 0;
    int i_max = a.size() - 1;

    //determine number of threads for addition
    if (a.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = a.size();
    }

    //choose interval size
    int i_interval = (a.size()/thread_count) + 1;
    int i_begin, i_end;

    //clear
    if (clear){
        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            threads[t] = std::thread(zeroS, std::ref(out), i_begin, i_end);
        }

        //join threads
        for (int t=0; t<thread_count; t++){
            threads[t].join();
        }
    }

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(subtractSfromS, std::ref(a), std::ref(b), scale, std::ref(out), i_begin, i_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    return;
}


#endif //MPM_V3_PARALLEL_EXPLICIT_USL_HPP

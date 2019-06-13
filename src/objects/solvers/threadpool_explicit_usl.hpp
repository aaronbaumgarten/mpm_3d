//
// Created by aaron on 6/10/19.
// threadpool_explicit_usl.hpp
//

#ifndef MPM_V3_THREADPOOL_EXPLICIT_USL_HPP
#define MPM_V3_THREADPOOL_EXPLICIT_USL_HPP

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

//scalar sparse matrix operations
void ThreadPoolExplicitUSL::ssmOperateStoSwithFlag(const MPMScalarSparseMatrix &S,
                                           const Eigen::VectorXd &x,
                                           Eigen::VectorXd &lhs,
                                           int SPEC, int k_begin, int k_end,
                                           bool& done){
    //do job
    ParallelExplicitUSL::ssmOperateStoS(S,x,lhs,SPEC,k_begin,k_end);
    //let everyone know it's done
    done = true;
    return;
}


void ThreadPoolExplicitUSL::ssmOperateVtoVwithFlag(const MPMScalarSparseMatrix &S,
                                           const KinematicVectorArray &x,
                                           KinematicVectorArray &lhs,
                                           int SPEC, int k_begin, int k_end,
                                           bool& done){
    //do job
    ParallelExplicitUSL::ssmOperateVtoV(S,x,lhs,SPEC,k_begin,k_end);
    //let everyone know it's done
    done = true;
    return;
}

//kinematic vector sparse matrix operations
void ThreadPoolExplicitUSL::kvsmLeftMultiplywithFlag(const KinematicVectorSparseMatrix &gradS,
                                             const MaterialTensorArray &T,
                                             MaterialVectorArray &lhs,
                                             int SPEC, int k_begin, int k_end,
                                             bool& done){
    //do job
    ParallelExplicitUSL::kvsmLeftMultiply(gradS, T, lhs, SPEC, k_begin, k_end);
    //say done
    done = true;
    return;
}

void ThreadPoolExplicitUSL::kvsmTensorProductTwithFlag(const KinematicVectorSparseMatrix &gradS,
                                               const KinematicVectorArray &x,
                                               KinematicTensorArray &L,
                                               int SPEC, int k_begin, int k_end,
                                               bool& done){
    //do job
    ParallelExplicitUSL::kvsmTensorProductT(gradS, x, L, SPEC, k_begin, k_end);
    //say done
    done = true;
    return;
}

//parallel scalar add
void ThreadPoolExplicitUSL::scalarAddwithFlag(const std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& list,
                                      Eigen::VectorXd &sum,
                                      int i_begin, int i_end, bool clear,
                                      bool& done){
    //do job
    ParallelExplicitUSL::scalarAdd(list, sum, i_begin, i_end, clear);
    //say done
    done = true;
    return;
}

void ThreadPoolExplicitUSL::vectorAddKwithFlag(const std::vector<KinematicVectorArray>& list,
                                       KinematicVectorArray &sum,
                                       int i_begin, int i_end, bool clear,
                                       bool& done){
    //do job
    ParallelExplicitUSL::vectorAddK(list, sum, i_begin, i_end, clear);
    //say done
    done = true;
    return;
}

void ThreadPoolExplicitUSL::vectorAddMwithFlag(const std::vector<MaterialVectorArray>& list,
                                       MaterialVectorArray &sum,
                                       int i_begin, int i_end, bool clear,
                                       bool& done){
    //do job
    ParallelExplicitUSL::vectorAddM(list, sum, i_begin, i_end, clear);
    //say done
    done = true;
    return;
}

void ThreadPoolExplicitUSL::tensorAddwithFlag(const std::vector<KinematicTensorArray>& list,
                                      KinematicTensorArray &sum,
                                      int i_begin, int i_end, bool clear,
                                      bool& done){
    //do job
    ParallelExplicitUSL::tensorAdd(list, sum, i_begin, i_end, clear);
    //say done
    done = true;
    return;
}

//vectorized operations
void ThreadPoolExplicitUSL::multiplyKVbySwithFlag(KinematicVectorArray &kv,
                                          Eigen::VectorXd &s,
                                          double scale,
                                          KinematicVectorArray &out,
                                          int i_begin, int i_end,
                                          bool& done){
    //do job
    ParallelExplicitUSL::multiplyKVbyS(kv, s, scale, out, i_begin, i_end);
    //say done
    done = true;
    return;
}

void ThreadPoolExplicitUSL::divideKVbySwithFlag(KinematicVectorArray &kv,
                                        Eigen::VectorXd &s,
                                        double scale,
                                        KinematicVectorArray &out,
                                        int i_begin, int i_end,
                                        bool& done){
    //do job
    ParallelExplicitUSL::divideKVbyS(kv, s, scale, out, i_begin, i_end);
    //say done
    done = true;
    return;
}

void ThreadPoolExplicitUSL::multiplySbySwithFlag(Eigen::VectorXd &a,
                                         Eigen::VectorXd &b,
                                         double scale,
                                         Eigen::VectorXd &out,
                                         int i_begin, int i_end,
                                         bool& done){
    //do job
    ParallelExplicitUSL::multiplySbyS(a, b, scale, out, i_begin, i_end);
    //say done
    done = true;
    return;
}

void ThreadPoolExplicitUSL::multiplyMTbySwithFlag(MaterialTensorArray &mt,
                                          Eigen::VectorXd &s,
                                          double scale,
                                          MaterialTensorArray &out,
                                          int i_begin, int i_end,
                                          bool& done){
    //do job
    ParallelExplicitUSL::multiplyMTbyS(mt, s, scale, out, i_begin, i_end);
    //say done
    done = true;
    return;
}

void ThreadPoolExplicitUSL::addStoSwithFlag(Eigen::VectorXd &a,
                                    Eigen::VectorXd &b,
                                    double scale,
                                    Eigen::VectorXd &out,
                                    int i_begin, int i_end,
                                    bool& done){
    //do job
    ParallelExplicitUSL::addStoS(a, b, scale, out, i_begin, i_end);
    //say done
    done = true;
    return;
}

void ThreadPoolExplicitUSL::subtractSfromSwithFlag(Eigen::VectorXd &a,
                                           Eigen::VectorXd &b,
                                           double scale,
                                           Eigen::VectorXd &out,
                                           int i_begin, int i_end,
                                           bool& done){
    //do job
    ParallelExplicitUSL::subtractSfromS(a, b, scale, out, i_begin, i_end);
    //say done
    done = true;
    return;
}

#endif //MPM_V3_THREADPOOL_EXPLICIT_USL_HPP

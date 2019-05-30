//
// Created by aaron on 5/27/19.
// parallel_explicit_usl.cpp
//

//
// Created by aaron on 5/16/18.
// explicit_usl.cpp
//

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

//scalar sparse matrix operations
void ParallelExplicitUSL::ssmOperate(const MPMScalarSparseMatrix &S,
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


void ParallelExplicitUSL::ssmOperate(const MPMScalarSparseMatrix &S,
                                     const KinematicVectorArray &x,
                                     KinematicVectorArray &lhs,
                                     int SPEC,
                                     int k_begin, int k_end){
    assert(S.cols(SPEC) == x.size() && "Scalar sparse matrix multiplication failed.");
    lhs = KinematicVectorArray(S.rows(SPEC),x.VECTOR_TYPE);
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

//kinematic vector sparse matrix operations
void ParallelExplicitUSL::kvsmLeftMultiply(const KinematicVectorSparseMatrix &gradS,
                                           const MaterialTensorArray &T,
                                           MaterialVectorArray &lhs,
                                           int SPEC,
                                           int k_begin, int k_end){
    assert(T.size() == gradS.cols(SPEC) && "Vector sparse matrix multiplication (tensor multiplication) failed.");
    lhs = MaterialVectorArray(gradS.rows(SPEC));
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

void ParallelExplicitUSL::kvsmTensorProductT(const KinematicVectorSparseMatrix &gradS,
                                             const KinematicVectorArray &x,
                                             KinematicTensorArray L,
                                             int SPEC,
                                             int k_begin, int k_end){
    assert(x.size() == gradS.cols(SPEC) && x.VECTOR_TYPE == gradS.VECTOR_TYPE && "Vector sparse matrix multiplication (tensor transpose product) failed.");
    L = KinematicTensorArray(gradS.rows(SPEC),gradS.VECTOR_TYPE);
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

void ParallelExplicitUSL::scalarAdd(const std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>>& list,
                                    Eigen::VectorXd &sum,
                                    int i_begin, int i_end){
    //fill in i_begin to i_end components of sum
    assert(sum.size() == list[0].size() && "scalarAdd failed.");
    //set to zero
    for (int i=i_begin; i<=i_end; i++){
        sum(i) = 0;
    }
    //add components
    for (int l=0; l<list.size(); l++){
        for (int i=i_begin; i<=i_end; i++){
            sum(i) += list[l](i);
        }
    }
    return;
}

void ParallelExplicitUSL::vectorAdd(const std::vector<KinematicVectorArray>& list,
                                    KinematicVectorArray &sum,
                                    int i_begin, int i_end){
    //fill in i_begin to i_end components of sum
    assert(sum.size() == list[0].size() && "vectorAdd failed.");
    //set to zero
    for (int i=i_begin; i<=i_end; i++){
        sum(i).setZero();
    }
    //add components
    for (int l=0; l<list.size(); l++){
        for (int i=i_begin; i<=i_end; i++){
            sum(i) += list[l](i);
        }
    }
    return;
}

void ParallelExplicitUSL::vectorAdd(const std::vector<MaterialVectorArray>& list,
                                    MaterialVectorArray &sum,
                                    int i_begin, int i_end){
    //fill in i_begin to i_end components of sum
    assert(sum.size() == list[0].size() && "vectorAdd failed.");
    //set to zero
    for (int i=i_begin; i<=i_end; i++){
        sum(i).setZero();
    }
    //add components
    for (int l=0; l<list.size(); l++){
        for (int i=i_begin; i<=i_end; i++){
            sum(i) += list[l](i);
        }
    }
    return;
}

void ParallelExplicitUSL::tensorAdd(const std::vector<KinematicTensorArray>& list,
                                    KinematicTensorArray &sum,
                                    int i_begin, int i_end){
    //fill in i_begin to i_end components of sum
    assert(sum.size() == list[0].size() && "tensorAdd failed.");
    //set to zero
    for (int i=i_begin; i<=i_end; i++){
        sum(i).setZero();
    }
    //add components
    for (int l=0; l<list.size(); l++){
        for (int i=i_begin; i<=i_end; i++){
            sum(i) += list[l](i);
        }
    }
    return;
}


void ParallelExplicitUSL::parallelMultiply(const MPMScalarSparseMatrix &S,
                                           const Eigen::VectorXd &x,
                                           Eigen::VectorXd &lhs,
                                           int SPEC){
    //get length of sparse matrix storage
    int k_max = S.size() - 1;

    //get length of output vector
    int i_max = lhs.rows() - 1;

    //determine number of threads for matrix vector mult
    int thread_count;
    if (S.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = S.size();
    }

    //intermediate storage vector
    std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd>> lhs_vec(thread_count);

    //choose interval size
    int k_interval = (int)ceil(S.size()/thread_count);
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
        threads[t] = std::thread(ssmOperate, S, x, lhs_vec[t], SPEC, k_begin, k_end);
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
    int i_interval = (int)ceil(lhs.rows()/thread_count);
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(scalarAdd,lhs_vec,lhs,i_begin,i_end);
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
                                           int SPEC){
    //get length of sparse matrix storage
    int k_max = S.size() - 1;

    //get length of output vector
    int i_max = lhs.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count;
    if (S.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = S.size();
    }

    //intermediate storage vector
    std::vector<KinematicVectorArray> lhs_vec(thread_count);

    //choose interval size
    int k_interval = (int)ceil(S.size()/thread_count);
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        lhs_vec[t] = Eigen::VectorXd(lhs.size());

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        threads[t] = std::thread(ssmOperate, S, x, lhs_vec[t], SPEC, k_begin, k_end);
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
    int i_interval = (int)ceil(lhs.size()/thread_count);
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(vectorAdd,lhs_vec,lhs,i_begin,i_end);
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
                                           int SPEC){
    //get length of sparse matrix storage
    int k_max = gradS.size() - 1;

    //get length of output vector
    int i_max = lhs.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count;
    if (gradS.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = gradS.size();
    }

    //intermediate storage vector
    std::vector<MaterialVectorArray> lhs_vec(thread_count);

    //choose interval size
    int k_interval = (int)ceil(gradS.size()/thread_count);
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        lhs_vec[t] = Eigen::VectorXd(lhs.size());

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        threads[t] = std::thread(kvsmLeftMultiply, gradS, T, lhs_vec[t], SPEC, k_begin, k_end);
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
    int i_interval = (int)ceil(lhs.size()/thread_count);
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(vectorAdd,lhs_vec,lhs,i_begin,i_end);
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
                                           KinematicTensorArray L,
                                           int SPEC){
    //get length of sparse matrix storage
    int k_max = gradS.size() - 1;

    //get length of output vector
    int i_max = L.size() - 1;

    //determine number of threads for matrix vector mult
    int thread_count;
    if (gradS.size() >= num_threads){
        thread_count = num_threads;
    } else {
        thread_count = gradS.size();
    }

    //intermediate storage vector
    std::vector<KinematicTensorArray> lhs_vec(thread_count);

    //choose interval size
    int k_interval = (int)ceil(gradS.size()/thread_count);
    int k_begin, k_end;

    for (int t=0; t<thread_count; t++) {
        //initialize lhs
        lhs_vec[t] = Eigen::VectorXd(L.size());

        //set interval
        k_begin = t * k_interval;
        k_end = k_begin + k_interval - 1;
        if (k_end > k_max){
            k_end = k_max;
        }
        //begin threads
        threads[t] = std::thread(kvsmTensorProductT, gradS, x, lhs_vec[t], SPEC, k_begin, k_end);
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
    int i_interval = (int)ceil(L.size()/thread_count);
    int i_begin, i_end;

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        threads[t] = std::thread(tensorAdd,lhs_vec,L,i_begin,i_end);
    }

    //join threads
    for (int t=0; t<thread_count; t++){
        threads[t].join();
    }

    //should be all done
    return;
}


/*----------------------------------------------------------------------------*/
//
void ParallelExplicitUSL::init(Job* job){
    //check that contact properties are set
    if (int_props.size() == 0){
        //do nothing
        cpdi_spec = DefaultBody::CPDI_ON;
        contact_spec = Contact::IMPLICIT;
        num_threads = 1;
    } if (int_props.size() == 1){
        //cpdi_spec given as argument
        cpdi_spec = int_props[0];
        contact_spec = Contact::IMPLICIT;
        num_threads = 1;
    } if (int_props.size() >= 2){
        //cpdi_spec and contact_spec given
        cpdi_spec = int_props[0];
        contact_spec = int_props[1];
        num_threads = 1;
    } if (int_props.size() >= 3){
        //cpdi_spec, contact_spec, and num_threads given
        cpdi_spec = int_props[0];
        contact_spec = int_props[1];
        num_threads = int_props[2];
    }

    //declare threads
    threads = std::vector<std::thread>(num_threads);

    printf("Solver properties (cpdi_spec = %i, contact_spec = %i, num_threads = %i).\n", cpdi_spec, contact_spec, num_threads);
    std::cout << "Solver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
std::string ParallelExplicitUSL::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

/*----------------------------------------------------------------------------*/
//
int ParallelExplicitUSL::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}


/*----------------------------------------------------------------------------*/
//

void ParallelExplicitUSL::mapPointsToNodes(Job* job){
    Body *body;
    Points *points;
    Nodes *nodes;
    Eigen::VectorXd pval;
    Eigen::VectorXd nval;
    KinematicVectorArray pvec;
    MaterialVectorArray nvec;
    MaterialTensorArray tmpMat;
    double tmpVAL;

    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = job->bodies[b].get();
        points = job->bodies[b]->points.get();
        nodes = job->bodies[b]->nodes.get();

        //map mass
        nodes->m = body->S * points->m; //m_i = S_ip * m_p

        //map momentum
        for (int i = 0; i < points->mx_t.size(); i++) {
            points->mx_t(i) = points->m(i) * points->x_t(i);
        }
        nodes->mx_t = body->S * points->mx_t;

        //calculate velocity
        for (int i = 0; i < nodes->x_t.size(); i++) {
            if (nodes->m(i) > 0) {
                nodes->x_t(i) = nodes->mx_t(i) / nodes->m(i);
            } else {
                nodes->x_t(i).setZero();
            }
        }

        //map body force
        pvec = KinematicVectorArray(points->b.size(), points->b.VECTOR_TYPE);
        for (int i = 0; i < points->b.size(); i++) {
            pvec(i) = points->m(i) * points->b(i);
        }
        nodes->f = body->S * pvec;

        //map divergence of stress
        tmpMat = points->T;
        for (int i = 0; i < tmpMat.size(); i++) {
            tmpMat[i] *= points->v[i];
        }
        nvec = body->gradS.left_multiply_by_tensor(tmpMat); //f_i = v_p T_p * gradS_ip
        for (int i = 0; i < nvec.size(); i++) {
            nodes->f[i] -= KinematicVector(nvec[i], nodes->f.VECTOR_TYPE);
        }

        if (job->JOB_TYPE == job->JOB_AXISYM){
            pval = Eigen::VectorXd(points->x.size());

            //scale s_tt by r and add contribution to f_r
            for (int i=0;i<pval.rows();i++){
                pval(i) = tmpMat(i,2,2) / points->x(i,0);
            }
            nval = body->S * pval;
            for (int i=0; i<nval.rows(); i++){
                nodes->f(i,0) -= nval(i);
            }

            //scale s_rt by r and add contribution to f_t
            for (int i=0; i<pval.rows(); i++){
                pval(i) = tmpMat(i,0,2) / points->x(i,0);
            }
            nval = body->S * pval;
            for (int i=0; i<nval.rows(); i++){
                nodes->f(i,2) += nval(i);
            }
        }
        //std::cout << points->m.sum() << " ?= " << nodes->m.sum() << std::endl;
    }
    return;
}

void ParallelExplicitUSL::moveGrid(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //update momentum
        job->bodies[b]->nodes->mx_t += job->dt * job->bodies[b]->nodes->f;

        //calculate velocity
        for (int i=0;i<job->bodies[b]->nodes->x_t.size();i++){
            if (job->bodies[b]->nodes->m(i) > 0) {
                job->bodies[b]->nodes->x_t(i) = job->bodies[b]->nodes->mx_t(i) / job->bodies[b]->nodes->m(i);
            } else {
                job->bodies[b]->nodes->x_t(i).setZero();
            }
        }

        //set displacement
        job->bodies[b]->nodes->u = job->dt * job->bodies[b]->nodes->x_t;

        //calculate difference in velocity
        for (int i=0;i<job->bodies[b]->nodes->diff_x_t.size();i++){
            if (job->bodies[b]->nodes->m(i) > 0) {
                job->bodies[b]->nodes->diff_x_t(i) = job->dt * job->bodies[b]->nodes->f(i) / job->bodies[b]->nodes->m(i);
            } else {
                job->bodies[b]->nodes->diff_x_t(i).setZero();
            }
        }
    }
    return;
}

void ParallelExplicitUSL::movePoints(Job* job){
    Body* body;
    Points* points;
    Nodes* nodes;
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = job->bodies[b].get();
        points = job->bodies[b]->points.get();
        nodes = job->bodies[b]->nodes.get();

        //map nodal displacement to point positions
        points->x += body->S.operate(nodes->u, MPMSparseMatrixBase::TRANSPOSED);
        points->u += body->S.operate(nodes->u, MPMSparseMatrixBase::TRANSPOSED);

        //fix position for out of plane dimension
        if (job->grid->GRID_DIM < job->DIM){
            for (int i=0; i<points->x.size(); i++){
                for (int pos=job->grid->GRID_DIM; pos<job->DIM; pos++){
                    points->x(i,pos) = 0;
                }
            }
        }

        //map nodal velocity diff to points
        points->x_t += body->S.operate(nodes->diff_x_t, MPMSparseMatrixBase::TRANSPOSED);

        //calculate momentum
        for (int i=0;i<points->mx_t.size();i++){
            points->mx_t(i) = points->m(i) * points->x_t(i);
        }
    }
    return;
}

void ParallelExplicitUSL::calculateStrainRate(Job* job){
    Body* body;
    Points* points;
    Nodes* nodes;
    KinematicVectorArray pvec;

    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = job->bodies[b].get();
        points = job->bodies[b]->points.get();
        nodes = job->bodies[b]->nodes.get();

        //calculate gradient of nodal velocity at points
        //body->bodyCalcPointGradient(job,points->L,nodes->x_t,Body::SET);
        points->L = body->gradS.tensor_product_transpose(nodes->x_t, MPMSparseMatrixBase::TRANSPOSED);

        //correct L if axisymmetric
        if (job->JOB_TYPE == job->JOB_AXISYM){
            pvec = body->S.operate(nodes->x_t, MPMSparseMatrixBase::TRANSPOSED);
            for (int i=0; i<points->L.size(); i++){
                points->L(i,0,2) = -pvec(i,2) / points->x(i,0);
                points->L(i,2,2) = pvec(i,0) / points->x(i,0);
            }
        }
    }
    return;
}

void ParallelExplicitUSL::updateDensity(Job* job){
    for (int b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (int i=0;i<job->bodies[b]->points->v.rows();i++) {
            job->bodies[b]->points->v(i) *= std::exp(job->dt * job->bodies[b]->points->L(i).trace());
        }

        //this is new, but maybe useful
        job->bodies[b]->points->updateIntegrators(job,job->bodies[b].get());
    }
    return;
}
//
// Created by aaron on 6/10/19.
// threadpool_explicit_usl.cpp
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
#include "threadpool_explicit_usl.hpp"

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <thread>
#include <time.h>
#include <functional>
#include <binders.h>

/*----------------------------------------------------------------------------*/
//
void ThreadPoolExplicitUSL::init(Job* job){
    //check that contact properties are set
    if (int_props.size() == 0){
        //do nothing
        cpdi_spec = DefaultBody::CPDI_ON;
        contact_spec = Contact::IMPLICIT;
        num_threads = job->thread_count;
    } else if (int_props.size() == 1){
        //cpdi_spec given as argument
        cpdi_spec = int_props[0];
        contact_spec = Contact::IMPLICIT;
        num_threads = job->thread_count;
    } else if (int_props.size() == 2){
        //cpdi_spec and contact_spec given
        cpdi_spec = int_props[0];
        contact_spec = int_props[1];
        num_threads = job->thread_count;
    } else if (int_props.size() >= 3){
        //cpdi_spec, contact_spec, num_threads and debug given
        cpdi_spec = int_props[0];
        contact_spec = int_props[1];
        num_threads = job->thread_count;
        if (int_props[2] == 1){
            debug = true;
        }
    }

    //get pointer to job->threadPool
    jobThreadPool = &job->threadPool;

    printf("Solver properties (cpdi_spec = %i, contact_spec = %i, num_threads = %i).\n", cpdi_spec, contact_spec, num_threads);
    std::cout << "Solver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
//
std::string ThreadPoolExplicitUSL::saveState(Job* job, Serializer* serializer, std::string filepath){
    return "err";
}

/*----------------------------------------------------------------------------*/
//
int ThreadPoolExplicitUSL::loadState(Job* job, Serializer* serializer, std::string fullpath){
    return 0;
}

/*----------------------------------------------------------------------------*/
void ThreadPoolExplicitUSL::parallelMultiply(const MPMScalarSparseMatrix &S,
                                           const Eigen::VectorXd &x,
                                           Eigen::VectorXd &lhs,
                                           int SPEC, bool clear){
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

    //boolean of completion status
    //std::vector<bool> taskComplete = std::vector<bool>(thread_count);
    /*for (int t=0; t<thread_count; t++){
        taskComplete[t] = false;
    }*/
    volatile bool firstTaskComplete[thread_count] = {false};

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
        //threads[t] = std::thread(ssmOperateStoS, std::ref(S), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end);
        //std::function<void (void)> tmpFunc = std::bind(ssmOperateStoSwithFlag, std::ref(S), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end, taskComplete[t]);
        jobThreadPool->doJob(std::bind(ssmOperateStoSwithFlag, std::ref(S), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end, std::ref(firstTaskComplete[t])));
        //std::cout << "Started task: " << t << "." << std::endl;
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            //std::cout << firstTaskComplete[t] << "?" << std::endl;
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }
    //std::cout << "Tasks completed." << std::endl;

    //determine number of threads for addition
    if (lhs.rows() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.rows();
    }

    //boolean for completion
    /*taskComplete = std::vector<bool>(thread_count);
    for (int t=0; t<thread_count; t++){
        taskComplete[t] = false;
    }*/
    volatile bool secondTaskComplete[thread_count] = {false};

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
        //threads[t] = std::thread(scalarAdd, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear);
        jobThreadPool->doJob(std::bind(scalarAddwithFlag, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear, std::ref(secondTaskComplete[t])));
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //should be all done
    return;
}

void ThreadPoolExplicitUSL::parallelMultiply(const MPMScalarSparseMatrix &S,
                                           const KinematicVectorArray &x,
                                           KinematicVectorArray &lhs,
                                           int SPEC, bool clear){
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

    //boolean of completion status
    volatile bool firstTaskComplete[thread_count] = {false};

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
        //threads[t] = std::thread(ssmOperateVtoV, std::ref(S), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end);
        jobThreadPool->doJob(std::bind(ssmOperateVtoVwithFlag, std::ref(S), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end, std::ref(firstTaskComplete[t])));
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //determine number of threads for addition
    if (lhs.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.size();
    }

    //boolean of completion status
    volatile bool secondTaskComplete[thread_count] = {false};

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
        //threads[t] = std::thread(vectorAddK, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear);
        jobThreadPool->doJob(std::bind(vectorAddKwithFlag, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear, std::ref(secondTaskComplete[t])));
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //should be all done
    return;
}

void ThreadPoolExplicitUSL::parallelMultiply(const KinematicVectorSparseMatrix &gradS,
                                           const MaterialTensorArray &T,
                                           MaterialVectorArray &lhs,
                                           int SPEC, bool clear){
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

    //boolean of completion status
    /*std::vector<bool> taskComplete = std::vector<bool>(thread_count);
    for (int t=0; t<thread_count; t++){
        taskComplete[t] = false;
    }*/
    volatile bool firstTaskComplete[thread_count] = {false};

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
        //threads[t] = std::thread(kvsmLeftMultiply, std::ref(gradS), std::ref(T), std::ref(lhs_vec[t]), SPEC, k_begin, k_end);
        jobThreadPool->doJob(std::bind(kvsmLeftMultiplywithFlag, std::ref(gradS), std::ref(T), std::ref(lhs_vec[t]), SPEC, k_begin, k_end, std::ref(firstTaskComplete[t])));
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //determine number of threads for addition
    if (lhs.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = lhs.size();
    }

    //boolean for task completion
    /*taskComplete = std::vector<bool>(thread_count);
    for (int t=0; t<thread_count; t++){
        taskComplete[t] = false;
    }*/
    volatile bool secondTaskComplete[thread_count] = {false};

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
        //threads[t] = std::thread(vectorAddM, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear);
        jobThreadPool->doJob(std::bind(vectorAddMwithFlag, std::ref(lhs_vec), std::ref(lhs), i_begin, i_end, clear, std::ref(secondTaskComplete[t])));
    }

    //join threads
    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //should be all done
    return;
}


void ThreadPoolExplicitUSL::parallelMultiply(const KinematicVectorSparseMatrix &gradS,
                                           const KinematicVectorArray &x,
                                           KinematicTensorArray &L,
                                           int SPEC, bool clear){
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

    //boolean array to track task completion
    volatile bool firstTaskComplete[thread_count] = {false};

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
        jobThreadPool->doJob(std::bind(kvsmTensorProductTwithFlag, std::ref(gradS), std::ref(x), std::ref(lhs_vec[t]), SPEC, k_begin, k_end, std::ref(firstTaskComplete[t])));
    }

    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!firstTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //determine number of threads for addition
    if (L.size() > num_threads){
        thread_count = num_threads;
    } else {
        thread_count = L.size();
    }

    //boolean array for second task
    volatile bool secondTaskComplete[thread_count] = {false};

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
        jobThreadPool->doJob(std::bind(tensorAddwithFlag, std::ref(lhs_vec), std::ref(L), i_begin, i_end, clear, std::ref(secondTaskComplete[t])));
    }

    //join threads
    //join threads
    /*for (int t=0; t<thread_count; t++){
        threads[t].join();
    }*/
    taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    //should be all done
    return;
}

/*----------------------------------------------------------------------------*/
//methods to manage threads during parallel operations

void ThreadPoolExplicitUSL::parallelMultiply(KinematicVectorArray &kv,
                                           Eigen::VectorXd &s,
                                           double scale,
                                           KinematicVectorArray &out,
                                           bool clear){
    //determine maximum index
    int thread_count;
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
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroKVwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        /*for (int t=0; t<thread_count; t++){
            threads[t].join();
        }*/
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //boolean array for second task
    volatile bool secondTaskComplete[thread_count] = {false};

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(multiplyKVbySwithFlag, std::ref(kv), std::ref(s), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    return;
}

void ThreadPoolExplicitUSL::parallelDivide(KinematicVectorArray &kv,
                                         Eigen::VectorXd &s,
                                         double scale,
                                         KinematicVectorArray &out,
                                         bool clear){
    //determine maximum index
    int thread_count;
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
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroKVwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //boolean array for second task completion
    volatile bool secondTaskComplete[thread_count] = {false};

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(divideKVbySwithFlag, std::ref(kv), std::ref(s), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    return;
}

void ThreadPoolExplicitUSL::parallelMultiply(Eigen::VectorXd &a,
                                           Eigen::VectorXd &b,
                                           double scale,
                                           Eigen::VectorXd &out,
                                           bool clear){
    //determine maximum index
    int thread_count;
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
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroSwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //boolean array for second task
    volatile bool secondTaskComplete[thread_count] = {false};

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(multiplySbySwithFlag, std::ref(a), std::ref(b), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    return;
}

void ThreadPoolExplicitUSL::parallelMultiply(MaterialTensorArray &mt,
                                           Eigen::VectorXd &s,
                                           double scale,
                                           MaterialTensorArray &out,
                                           bool clear){
    //determine maximum index
    int thread_count;
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
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroMTwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //second boolean array
    volatile bool secondTaskComplete[thread_count] = {false};

    //multipy MaterialTensor by Scalar and add to output
    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(multiplyMTbySwithFlag, std::ref(mt), std::ref(s), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    return;
}

void ThreadPoolExplicitUSL::parallelAdd(Eigen::VectorXd &a,
                                      Eigen::VectorXd &b,
                                      double scale,
                                      Eigen::VectorXd &out,
                                      bool clear){

    //determine maximum index
    int thread_count;
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
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroSwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //second task completion boolean array
    volatile bool secondTaskComplete[thread_count] = {false};

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(addStoSwithFlag, std::ref(a), std::ref(b), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    return;
}

void ThreadPoolExplicitUSL::parallelSubtract(Eigen::VectorXd &a,
                                           Eigen::VectorXd &b,
                                           double scale,
                                           Eigen::VectorXd &out,
                                           bool clear){
    //determine maximum index
    int thread_count;
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
        volatile bool firstTaskComplete[thread_count] = {false};

        for (int t=0; t<thread_count; t++){
            //set interval
            i_begin = t*i_interval;
            i_end = i_begin + i_interval-1;
            if (i_end > i_max){
                i_end = i_max;
            }
            //begin threads
            jobThreadPool->doJob(std::bind(zeroSwithFlag, std::ref(out), i_begin, i_end, std::ref(firstTaskComplete[t])));
        }

        //join threads
        bool taskDone = false;
        //wait for task to complete
        while (!taskDone){
            //set flag to true
            taskDone = true;
            for (int t=0; t<thread_count; t++){
                if (!firstTaskComplete[t]){
                    //if any task is not done, set flag to false
                    taskDone = false;
                    break;
                }
            }
        }
    }

    //second boolean array for task completion
    volatile bool secondTaskComplete[thread_count] = {false};

    for (int t=0; t<thread_count; t++){
        //set interval
        i_begin = t*i_interval;
        i_end = i_begin + i_interval-1;
        if (i_end > i_max){
            i_end = i_max;
        }
        //begin threads
        jobThreadPool->doJob(std::bind(subtractSfromSwithFlag, std::ref(a), std::ref(b), scale, std::ref(out), i_begin, i_end, std::ref(secondTaskComplete[t])));
    }

    //join threads
    bool taskDone = false;
    //wait for task to complete
    while (!taskDone){
        //set flag to true
        taskDone = true;
        for (int t=0; t<thread_count; t++){
            if (!secondTaskComplete[t]){
                //if any task is not done, set flag to false
                taskDone = false;
                break;
            }
        }
    }

    return;
}
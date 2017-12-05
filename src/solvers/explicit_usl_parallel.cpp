//
// Created by aaron on 8/18/17.
// explicit_usl_parallel.cpp
//

//
// Created by aaron on 5/25/17.
// explicit_usl.cpp
//

#include <stdlib.h>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include <pthread.h>
#include <atomic>
#include <driver.hpp>

#include "job.hpp"
#include "serializer.hpp"

#include "solver.hpp"

#include "grid.hpp"
#include "contact.hpp"

#include "body.hpp"
#include "nodes.hpp"
#include "points.hpp"

#include "material.hpp"
#include "boundary.hpp"

int use_cpdi;
int num_threads;
int phi_iterator;
int gradphi_iterator;
std::atomic<bool> use_gradphi_pushback, use_phi_pushback;
pthread_mutex_t phi_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t gradphi_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t phi_barrier, gradphi_barrier;

extern "C" void solverInit(Job* job); //initialize solver
extern "C" void solverStep(Job* job); //one forward mpm step

extern "C" std::string solverSaveState(Job* job, Serializer* serializer, std::string filepath); //save solver state to returned filename in serializer folder
extern "C" int solverLoadState(Job* job, Serializer* serializer, std::string fullpath); //load state from given full path

/*----------------------------------------------------------------------------*/

struct threadtask_t {
    Job *job;

} threadtask_s;

struct mappingtask_t {
    Job* job;
    Body* body;

    int p_blockstart;
    int p_blocksize;

}mappingtask_s;

/*----------------------------------------------------------------------------*/

void solverInit(Job* job){
    if (job->solver.int_props.size() < 2) {
        std::cout << job->solver.int_props.size() << "\n";
        fprintf(stderr,
                "%s:%s: Need at least 2 property defined (num_threads, use_cpdi).\n",
                __FILE__, __func__);
        exit(0);
    } else {
        //store stop_time
        num_threads = job->solver.int_props[0];
        use_cpdi = job->solver.int_props[1];
    }

    std::cout << "Solver Initialized." << std::endl;
    return;
}

/*----------------------------------------------------------------------------*/
void * parallel_map_method(void* _task){
    mappingtask_t *task = (mappingtask_t*)_task;

    //calculate phi and grad phi
    std::vector<int> nvec(0);
    std::vector<double> valvec(0);
    std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>> gradvec(0);
    Eigen::VectorXd tmpGrad = task->job->jobVector<double>();
    int ith_cpdi;
    int kth_update;

    //for (size_t i = 0; i < task->body->points.x.rows(); i++) {
    for (size_t i = task->p_blockstart; i < (task->p_blockstart + task->p_blocksize); i++) {
        if (use_cpdi == 1) {
            ith_cpdi = Body::CPDI_ON;
            //check whether point is active:
            if (task->body->points.active(i) == 0) {
                continue;
            } else if (!task->job->grid.gridInDomain(task->job, task->body->points.x.row(i).transpose())) {
                task->body->points.active(i) = 0;
                continue;
            } else if (ith_cpdi == Body::CPDI_ON) {
                //check that corners are in domain
                for (size_t c = 0; c < task->body->A.rows(); c++) {
                    if (!task->job->grid.gridInDomain(task->job, (task->body->points.x.row(i) +
                                                                  0.5 * task->body->points.extent(i) *
                                                                  task->body->A.row(c).cast<double>()).transpose())) {
                        //corner out of domain
                        ith_cpdi = Body::CPDI_OFF;
                        break;
                    }
                }
            }
        } else {
            ith_cpdi = Body::CPDI_OFF;
        }

        //calculate map for ith point
        if (ith_cpdi == Body::CPDI_OFF) {
            nvec.resize(0);
            valvec.resize(0);
            task->job->grid.gridEvaluateShapeFnValue(task->job, task->body->points.x.row(i).transpose(), nvec, valvec);

            for (size_t j=0; j< nvec.size(); j++) {
                //get lock on phi
                pthread_mutex_lock(&phi_mutex);
                //set k
                kth_update = phi_iterator;
                phi_iterator += 1;
                //unlock phi
                pthread_mutex_unlock(&phi_mutex);

                if (kth_update < task->body->phi.size()) {
                    //set kth value
                    task->body->pval[kth_update] = i;
                    task->body->nval[kth_update] = nvec[j];
                    task->body->phi[kth_update] = valvec[j];
                } else if (use_phi_pushback == false) {
                    //wait for all threads to catch up
                    //std::cout << "Waiting..." << std::endl;
                    pthread_barrier_wait(&phi_barrier);
                    //std::cout << "Released." << std::endl;
                    //set pushback to true
                    use_phi_pushback = true;
                    //get lock
                    pthread_mutex_lock(&phi_mutex);
                    //pushback values
                    task->body->pval.push_back((int) i);
                    task->body->nval.push_back(nvec[j]);
                    task->body->phi.push_back(valvec[j]);
                    //unlock
                    pthread_mutex_unlock(&phi_mutex);
                } else if (use_phi_pushback == true) {
                    //get lock
                    pthread_mutex_lock(&phi_mutex);
                    //pushback values
                    task->body->pval.push_back((int) i);
                    task->body->nval.push_back(nvec[j]);
                    task->body->phi.push_back(valvec[j]);
                    //unlock
                    pthread_mutex_unlock(&phi_mutex);
                }
            }
            /*
            //write phi
            for (size_t j = 0; j < nvec.size(); j++) {
                //write body map
                task->body->pval.push_back((int) i);
                task->body->nval.push_back(nvec[j]);
                task->body->phi.push_back(valvec[j]);
            }
            //release lock on phi
            pthread_mutex_unlock(&phi_mutex);
             */
        } else if (ith_cpdi == Body::CPDI_ON) {
            nvec.resize(0);
            valvec.resize(0);
            for (size_t c = 0; c < task->body->A.rows(); c++) {
                //spread influence
                task->job->grid.gridEvaluateShapeFnValue(task->job, (task->body->points.x.row(i) + 0.5 * task->body->points.extent(i) *
                                                                                                  task->body->A.row(c).cast<double>()).transpose(),
                                                   nvec, valvec);
            }

            for (size_t j=0; j< nvec.size(); j++) {
                //get lock on phi
                pthread_mutex_lock(&phi_mutex);
                //set k
                kth_update = phi_iterator;
                phi_iterator += 1;
                //unlock phi
                pthread_mutex_unlock(&phi_mutex);

                if (kth_update < task->body->phi.size()) {
                    //set kth value
                    task->body->pval[kth_update] = i;
                    task->body->nval[kth_update] = nvec[j];
                    task->body->phi[kth_update] = valvec[j]/task->body->A.rows();
                } else if (use_phi_pushback == false) {
                    //wait for all threads to catch up
                    //std::cout << "Waiting..." << std::endl;
                    pthread_barrier_wait(&phi_barrier);
                    //std::cout << "Released." << std::endl;
                    //set pushback to true
                    use_phi_pushback = true;
                    //get lock
                    pthread_mutex_lock(&phi_mutex);
                    //pushback values
                    task->body->pval.push_back((int) i);
                    task->body->nval.push_back(nvec[j]);
                    task->body->phi.push_back(valvec[j]/task->body->A.rows());
                    //unlock
                    pthread_mutex_unlock(&phi_mutex);
                } else if (use_phi_pushback == true) {
                    //get lock
                    pthread_mutex_lock(&phi_mutex);
                    //pushback values
                    task->body->pval.push_back((int) i);
                    task->body->nval.push_back(nvec[j]);
                    task->body->phi.push_back(valvec[j]/task->body->A.rows());
                    //unlock
                    pthread_mutex_unlock(&phi_mutex);
                }
            }
            /*
            //get lock on phi
            pthread_mutex_lock(&phi_mutex);
            //write phi
            for (size_t j = 0; j < nvec.size(); j++) {
                //write body map
                task->body->pval.push_back((int) i);
                task->body->nval.push_back(nvec[j]);
                task->body->phi.push_back(valvec[j]/task->body->A.rows());
            }
            //release lock on phi
            pthread_mutex_unlock(&phi_mutex);
             */
        }
    }

    if (use_phi_pushback == false) {
        //wait for all threads to catch up
        //std::cout << "Waiting..." << std::endl;
        pthread_barrier_wait(&phi_barrier);
        //std::cout << "Released." << std::endl;
    }

    for (size_t i = task->p_blockstart; i < (task->p_blockstart + task->p_blocksize); i++) {
        if (use_cpdi == 1) {
            ith_cpdi = Body::CPDI_ON;
            //check whether point is active:
            if (task->body->points.active(i) == 0) {
                continue;
            } else if (!task->job->grid.gridInDomain(task->job, task->body->points.x.row(i).transpose())) {
                task->body->points.active(i) = 0;
                continue;
            } else if (ith_cpdi == Body::CPDI_ON) {
                //check that corners are in domain
                for (size_t c = 0; c < task->body->A.rows(); c++) {
                    if (!task->job->grid.gridInDomain(task->job, (task->body->points.x.row(i) +
                                                                  0.5 * task->body->points.extent(i) *
                                                                  task->body->A.row(c).cast<double>()).transpose())) {
                        //corner out of domain
                        ith_cpdi = Body::CPDI_OFF;
                        break;
                    }
                }
            }
        } else {
            ith_cpdi = Body::CPDI_OFF;
        }

        //calculate map for ith point
        if (ith_cpdi == Body::CPDI_OFF) {
            nvec.resize(0);
            gradvec.resize(0);
            task->job->grid.gridEvaluateShapeFnGradient(task->job, task->body->points.x.row(i).transpose(), nvec, gradvec);

            for (size_t j=0; j< nvec.size(); j++) {
                //get lock on gradphi
                pthread_mutex_lock(&gradphi_mutex);
                //set k
                kth_update = gradphi_iterator;
                gradphi_iterator += 1;
                //unlock phi
                pthread_mutex_unlock(&gradphi_mutex);

                if (kth_update < task->body->gradphi.size()) {
                    //set kth value
                    task->body->pgrad[kth_update] = i;
                    task->body->ngrad[kth_update] = nvec[j];
                    task->body->gradphi[kth_update] = gradvec[j];
                } else if (use_gradphi_pushback == false) {
                    //wait for all threads to catch up
                    //std::cout << "Waiting..." << std::endl;
                    pthread_barrier_wait(&gradphi_barrier);
                    //std::cout << "Released." << std::endl;
                    //set pushback to true
                    use_gradphi_pushback = true;
                    //get lock
                    pthread_mutex_lock(&gradphi_mutex);
                    //pushback values
                    task->body->pgrad.push_back((int) i);
                    task->body->ngrad.push_back(nvec[j]);
                    task->body->gradphi.push_back(gradvec[j]);
                    //unlock
                    pthread_mutex_unlock(&gradphi_mutex);
                } else if (use_gradphi_pushback == true) {
                    //get lock
                    pthread_mutex_lock(&gradphi_mutex);
                    //pushback values
                    task->body->pgrad.push_back((int) i);
                    task->body->ngrad.push_back(nvec[j]);
                    task->body->gradphi.push_back(gradvec[j]);
                    //unlock
                    pthread_mutex_unlock(&gradphi_mutex);
                }
            }
            /*
            //get lock on gradphi
            pthread_mutex_lock(&gradphi_mutex);
            //write graphi
            for (size_t j = 0; j < nvec.size(); j++) {
                task->body->pgrad.push_back((int) i);
                task->body->ngrad.push_back(nvec[j]);
                task->body->gradphi.push_back(gradvec[j]);
            }
            //release lock on graphi
            pthread_mutex_unlock(&gradphi_mutex);
             */
        } else if (ith_cpdi == Body::CPDI_ON) {
            //average gradients along sides of extent
            nvec.resize(0);
            valvec.resize(0);
            gradvec.resize(0);
            for (size_t c = 0; c < task->body->A.rows(); c++) {
                //find shape function value at offset point location
                //add node ids to nodevec
                //add values to valvec
                valvec.resize(0);
                task->job->grid.gridEvaluateShapeFnValue(task->job, (task->body->points.x.row(i) + 0.5 * task->body->points.extent(i) *
                                                                                                   task->body->A.row(c).cast<double>()).transpose(),
                                                         nvec, valvec);
                for (size_t v = 0; v < valvec.size(); v++) {
                    //gradient contribution from corner
                    //G(x) = (S(x+a) - S(x-a))/(2a)
                    tmpGrad = valvec[v] / (task->body->A.rows() * 0.5 * task->body->points.extent(i)) * task->body->A.row(c).cast<double>().transpose();
                    gradvec.push_back(tmpGrad);
                }
            }
            for (size_t j=0; j< nvec.size(); j++) {
                //get lock on gradphi
                pthread_mutex_lock(&gradphi_mutex);
                //set k
                kth_update = gradphi_iterator;
                gradphi_iterator += 1;
                //unlock phi
                pthread_mutex_unlock(&gradphi_mutex);

                if (kth_update < task->body->gradphi.size()) {
                    //set kth value
                    task->body->pgrad[kth_update] = i;
                    task->body->ngrad[kth_update] = nvec[j];
                    task->body->gradphi[kth_update] = gradvec[j];
                } else if (use_gradphi_pushback == false) {
                    //wait for all threads to catch up
                    //std::cout << "Waiting..." << std::endl;
                    pthread_barrier_wait(&gradphi_barrier);
                    //std::cout << "Released." << std::endl;
                    //set pushback to true
                    use_gradphi_pushback = true;
                    //get lock
                    pthread_mutex_lock(&gradphi_mutex);
                    //pushback values
                    task->body->pgrad.push_back((int) i);
                    task->body->ngrad.push_back(nvec[j]);
                    task->body->gradphi.push_back(gradvec[j]);
                    //unlock
                    pthread_mutex_unlock(&gradphi_mutex);
                } else if (use_gradphi_pushback == true) {
                    //get lock
                    pthread_mutex_lock(&gradphi_mutex);
                    //pushback values
                    task->body->pgrad.push_back((int) i);
                    task->body->ngrad.push_back(nvec[j]);
                    task->body->gradphi.push_back(gradvec[j]);
                    //unlock
                    pthread_mutex_unlock(&gradphi_mutex);
                }
            }/*
            //get lock on gradphi
            pthread_mutex_lock(&gradphi_mutex);
            //write graphi
            for (size_t j = 0; j < nvec.size(); j++) {
                task->body->pgrad.push_back((int) i);
                task->body->ngrad.push_back(nvec[j]);
                task->body->gradphi.push_back(gradvec[j]);
            }
            //release lock on graphi
            pthread_mutex_unlock(&gradphi_mutex);
            */
        }
    }

    if (use_gradphi_pushback == false) {
        //wait for all threads to catch up
        //std::cout << "Waiting..." << std::endl;
        pthread_barrier_wait(&gradphi_barrier);
        //std::cout << "Released." << std::endl;
    }

    return NULL;
}


void createMappings(Job* job) {
    //setup threads
    std::vector<mappingtask_t> tasks(num_threads);
    std::vector<pthread_t> threads(num_threads-1);

    int max_p_blocksize;
    int min_p_blocksize;

    for (size_t b=0; b<job->bodies.size(); b++) {
        //clear prior map
        /*
        job->bodies[b].pval.clear();
        job->bodies[b].nval.clear();
        job->bodies[b].phi.clear();

        job->bodies[b].pgrad.clear();
        job->bodies[b].ngrad.clear();
        job->bodies[b].gradphi.clear();
         */

        //set size of map to zero
        phi_iterator = 0;
        use_phi_pushback = false;
        gradphi_iterator = 0;
        use_gradphi_pushback = false;

        //divide point domain
        max_p_blocksize = std::ceil(job->bodies[b].points.m.size() / num_threads);

        //setup thread structure
        for (size_t i = 0; i < num_threads; i++) {
            tasks[i].job = job;
            tasks[i].body = &(job->bodies[b]);

            tasks[i].p_blockstart = i*max_p_blocksize;
            min_p_blocksize = job->bodies[b].points.m.size() - i*max_p_blocksize;
            tasks[i].p_blocksize = std::min(max_p_blocksize,min_p_blocksize);
        }

        //initialize barrier
        pthread_barrier_init(&phi_barrier,NULL,num_threads);
        pthread_barrier_init(&gradphi_barrier,NULL,num_threads);

        //create threads
        for (size_t i = 0; i < num_threads - 1; i++) {
            pthread_create(&(threads[i]), NULL, &parallel_map_method, &(tasks[i]));
        }
        parallel_map_method(&(tasks[num_threads - 1]));

        for (size_t i = 0; i < num_threads - 1; i++) {
            //wait for thread to end
            //std::cout << "where is thread " << i << std::endl;
            pthread_join(threads[i], NULL);
            //std::cout << "found " << i << std::endl;
        }

        //resize vector if larger than iterator
        if (job->bodies[b].nval.size() > phi_iterator){
            job->bodies[b].pval.resize(phi_iterator);
            job->bodies[b].nval.resize(phi_iterator);
            job->bodies[b].phi.resize(phi_iterator);
        }
        if (job->bodies[b].ngrad.size() > gradphi_iterator){
            job->bodies[b].pgrad.resize(gradphi_iterator);
            job->bodies[b].ngrad.resize(gradphi_iterator);
            job->bodies[b].gradphi.resize(gradphi_iterator);
        }

        //destroy barriers
        pthread_barrier_destroy(&phi_barrier);
        pthread_barrier_destroy(&gradphi_barrier);

    }

    return;
}

void mapPointsToNodes(Job* job){
    Body *body;
    Points *points;
    Nodes *nodes;
    Eigen::MatrixXd pvec;
    Eigen::MatrixXd nvec;
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = &(job->bodies[b]);
        points = &(job->bodies[b].points);
        nodes = &(job->bodies[b].nodes);

        //map mass
        body->bodyCalcNodalValues(job,nodes->m,points->m,Body::SET);

        //map momentum
        for (size_t i=0;i<points->mx_t.cols();i++){
            points->mx_t.col(i) = points->m.array() * points->x_t.col(i).array();
        }
        body->bodyCalcNodalValues(job,nodes->mx_t,points->mx_t,Body::SET);

        //calculate velocity
        for (size_t i=0;i<nodes->x_t.rows();i++){
            if (nodes->m(i) > 0){
                nodes->x_t.row(i) = nodes->mx_t.row(i) / nodes->m(i);
            } else {
                nodes->x_t.row(i).setZero();
            }
        }

        //map body force
        pvec = job->jobVectorArray<double>(points->b.rows());
        for (size_t i=0;i<points->b.cols();i++){
            pvec.col(i) = points->m.array() * points->b.col(i).array();
        }
        body->bodyCalcNodalValues(job,nodes->f,pvec,Body::SET);

        //map divergence of stress
        body->bodyCalcNodalDivergence(job,nodes->f,points->T,Body::ADD);
        //nvec = job->jobVectorArray<double>(nodes->f.rows());
        //body->bodyCalcNodalDivergence(job,nvec,points->T);
        //nodes->f += nvec;
    }
    return;
}

void generateContacts(Job* job){
    for (size_t c=0;c<job->contacts.size();c++){
        if (job->activeContacts[c] == 0){
            continue;
        }
        job->contacts[c].contactGenerateRules(job);
    }
    return;
}

void addContacts(Job* job){
    for (size_t c=0;c<job->contacts.size();c++){
        if (job->activeContacts[c] == 0){
            continue;
        }
        job->contacts[c].contactApplyRules(job,Contact::EXPLICIT);
    }
    return;
}

void generateBoundaryConditions(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b].activeBoundary == 0){
            continue;
        }
        job->bodies[b].boundary.boundaryGenerateRules(job,&(job->bodies[b]));
    }
    return;
}

void addBoundaryConditions(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b].activeBoundary == 0){
            continue;
        }
        job->bodies[b].boundary.boundaryApplyRules(job,&(job->bodies[b]));
    }
    return;
}

void moveGrid(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }

        //update momentum
        job->bodies[b].nodes.mx_t += job->dt * job->bodies[b].nodes.f;

        //calculate velocity
        for (size_t i=0;i<job->bodies[b].nodes.x_t.rows();i++){
            if (job->bodies[b].nodes.m(i) > 0) {
                job->bodies[b].nodes.x_t.row(i) = job->bodies[b].nodes.mx_t.row(i) / job->bodies[b].nodes.m(i);
            } else {
                job->bodies[b].nodes.x_t.row(i).setZero();
            }
        }

        //set displacement
        job->bodies[b].nodes.u = job->dt * job->bodies[b].nodes.x_t;

        //calculate difference in velocity
        for (size_t i=0;i<job->bodies[b].nodes.diff_x_t.rows();i++){
            if (job->bodies[b].nodes.m(i) > 0) {
                job->bodies[b].nodes.diff_x_t.row(i) = job->dt * job->bodies[b].nodes.f.row(i) / job->bodies[b].nodes.m(i);
            } else {
                job->bodies[b].nodes.diff_x_t.row(i).setZero();
            }
        }
    }
    return;
}

void movePoints(Job* job){
    Body* body;
    Points* points;
    Nodes* nodes;
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = &(job->bodies[b]);
        points = &(job->bodies[b].points);
        nodes = &(job->bodies[b].nodes);

        //map nodal displacement to point positions
        body->bodyCalcPointValues(job,points->x,nodes->u,Body::ADD);
        body->bodyCalcPointValues(job,points->u,nodes->u,Body::ADD);

        //map nodal velocity diff to points
        body->bodyCalcPointValues(job,points->x_t,nodes->diff_x_t,Body::ADD);

        //calculate momentum
        for (size_t i=0;i<points->mx_t.cols();i++){
            points->mx_t.col(i) = points->m.array() * points->x_t.col(i).array();
        }
    }
    return;
}

void calculateStrainRate(Job* job){
    Body* body;
    Points* points;
    Nodes* nodes;
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        body = &(job->bodies[b]);
        points = &(job->bodies[b].points);
        nodes = &(job->bodies[b].nodes);

        //calculate gradient of nodal velocity at points
        body->bodyCalcPointGradient(job,points->L,nodes->x_t,Body::SET);
    }
    return;
}

void updateDensity(Job* job){
    Eigen::MatrixXd L;
    Eigen::VectorXd tmpVec;

    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0){
            continue;
        }
        for (size_t i=0;i<job->bodies[b].points.v.rows();i++) {
            tmpVec = job->bodies[b].points.L.row(i).transpose();
            L = job->jobTensor<double>(tmpVec.data());
            job->bodies[b].points.v(i) *= std::exp(job->dt * L.trace());
        }
    }
    return;
}

void updateStress(Job* job){
    for (size_t b=0;b<job->bodies.size();b++){
        if (job->activeBodies[b] == 0 || job->bodies[b].activeMaterial == 0){
            continue;
        }
        job->bodies[b].material.materialCalculateStress(job,&(job->bodies[b]),Material::UPDATE);
    }
    return;
}

/*----------------------------------------------------------------------------*/

void solverStep(Job* job){

    //setup threads
    /*tasks = std::vector<threadtask_t>(num_threads);
    threads = std::vector<pthread_t> (num_threads-1);

    for (size_t i=0;i<num_threads;i++){
        tasks[i].job = job;
    }

    for (size_t i=0;i<num_threads-1;i++){
        pthread_create(&(threads[i]),NULL,&test_pthread,&(tasks[i]));
    }
    test_pthread(&(tasks[num_threads-1]));
    for (size_t i=0;i<num_threads-1;i++){
        pthread_join(threads[i], NULL);
    }*/


    //this can be parallelized as its own functions
    //create map
    createMappings(job);

    //this can be split between bodies
    //map particles to grid
    mapPointsToNodes(job);

    //this must be run in series
    //add contact forces
    generateContacts(job);
    addContacts(job);

    //this can be split between bodies (in theory)
    //enforce boundary conditions
    generateBoundaryConditions(job);
    addBoundaryConditions(job);

    //this can be split between bodies
    //move grid
    moveGrid(job);

    //this can be split between bodies
    //move particles
    movePoints(job);

    //maybe parallelize if i do the math instead of the body?

    //this can be split between bodies
    //calculate strainrate
    calculateStrainRate(job);

    //this can be split between bodies
    //update density
    updateDensity(job);

    //add body forces
    job->driver.driverGenerateGravity(job);
    job->driver.driverApplyGravity(job);

    //update stress
    updateStress(job);

    return;
}

std::string solverSaveState(Job* job, Serializer* serializer, std::string filepath){
    // current date/time based on current system
    time_t now = time(0);

    // convert now to tm struct for UTC
    tm *gmtm = gmtime(&now);
    std::string filename = "ERR";

    //create filename/directory name
    std::stringstream s;
    s << "mpm_v2.solver." << gmtm->tm_mday << "." << gmtm->tm_mon << "." << gmtm->tm_year << ".";
    s << gmtm->tm_hour << "." << gmtm->tm_min << "." << gmtm->tm_sec << ".txt";

    //create filename
    filename = s.str();

    std::ofstream ffile((filepath+filename), std::ios::trunc);

    //write data
    if (ffile.is_open()) {
        ffile << "# mpm_v2 Solver\n";
        ffile << num_threads << "\n";
        ffile.close();
    } else {
        std::cout << "Unable to open \"" << filename << "\" !\n";
        return "ERR";
    }

    return filename;
}

int solverLoadState(Job* job, Serializer* serializer, std::string fullpath){
    std::string line; //read line

    std::ifstream fin(fullpath); //file to load from

    if (fin.is_open()) {
        //if open, read lines
        std::getline(fin,line); //header
        std::getline(fin,line); //num_threads
        num_threads = std::stoi(line);
        fin.close();
    } else {
        std::cout << "ERROR: Unable to open file: " << fullpath << std::endl;
        return 0;
    }

    std::cout << "Solver Loaded." << std::endl;
    return 1;
}

//
// Created by aaron on 4/23/18.
// test.hpp
//

#ifndef MPM_V3_TEST_HPP
#define MPM_V3_TEST_HPP

//test linear algebra functions
void algebra_test();
void map_test();

class Job;
void fvm_test(Job*);
void fvm_mpm_drag_test(Job*);
void fvm_mpm_buoyancy_test(Job*);
void fvm_mpm_porosity_test(Job*);

#endif //MPM_V3_TEST_HPP

//
// Created by aaron on 8/20/19.
// factorial.hpp
//

#ifndef MPM_V3_FACTORIAL_HPP
#define MPM_V3_FACTORIAL_HPP

inline int factorial(int n){
    if (n < 0){
        return -1;
    } else if (n == 0){
        return 1;
    } else {
        return n*factorial(n-1);
    }
}

#endif //MPM_V3_FACTORIAL_HPP

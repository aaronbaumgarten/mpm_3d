//
// Created by aaron on 12/3/18.
// map_test.cpp
//


#include "mpm_tensor.hpp"
#include "mpm_tensorarray.hpp"
#include "mpm_vector.hpp"
#include "mpm_vectorarray.hpp"
#include <stdlib.h>
#include "test.hpp"

void map_test(){
    int vector_type = KinematicVector::VECTOR_3D;

    KinematicVector a = KinematicVector(vector_type);
    KinematicVector b = KinematicVector(vector_type);
    KinematicVector c = KinematicVector(vector_type);
    KinematicVector d = KinematicVector(vector_type);

    a[0] = 1; a[1] = 2; a[2] = -1.414;
    b[0] = -1; b[1] = 2.5; b[2] = 7;
    c[0] = 0; c[1] = 1; c[2] = 18;
    d[0] = 14; d[1] = -2; d[2] = 0.001;

    MaterialVectorArray AB = MaterialVectorArray(0);
    KinematicVectorArray CD = KinematicVectorArray(0,vector_type);

    AB.push_back(a);
    AB.push_back(b);

    CD.push_back(c);
    CD.push_back(d);

    AB[1] = CD[1];
    std::cout << EIGEN_MAP_OF_MATERIAL_VECTOR(AB[1]).transpose() << std::endl;
    std::cout << EIGEN_MAP_OF_KINEMATIC_VECTOR(CD[1]).transpose() << std::endl;

    AB[1] = d;
    std::cout << EIGEN_MAP_OF_MATERIAL_VECTOR(AB[1]).transpose() << std::endl;

    return;
}
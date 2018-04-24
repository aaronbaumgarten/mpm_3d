//
// Created by aaron on 4/23/18.
// algebra_test.cpp
//

#include "mpmtensor.hpp"
#include "mpmtensorarray.hpp"
#include "mpmvector.hpp"
#include "mpmvectorarray.hpp"
#include <stdlib.h>
#include "test.hpp"

void algebra_test(){
    //an arbitrary set of tensors and vectors
    std::array<double,9> a = {{1, 2, -2, 4, 5, 7, 7, -8, 10}};
    std::array<double,9> b = {{0.6, 10.5, 13, 0.01, -13.2, 1, 1, 9.3, 7}};
    std::array<double,3> u = {{1, 2, 3}};
    std::array<double,3> v = {{10.1, 0.3, -4}};

    //choose dimension of kinematic tensors
    int tensor_type = KinematicTensor::TENSOR_2D;

    /*------------------------------------------------------------------------*/
    //test initialization of tensors
    std::cout << "Checking initialization of tensors." << std::endl;

    KinematicTensor A = KinematicTensor(a.data(),tensor_type);
    assert(A(0,1) == a[1] && "KinematicTensor initialization failed.");

    MaterialTensor S = MaterialTensor(b.data());
    assert(S(2,1) == b[7] && "MaterialTensor initialization failed.");


    /*------------------------------------------------------------------------*/
    //test copying of tensors
    std::cout << "Checking copying of tensors." << std::endl;

    KinematicTensor B;
    B = A;
    assert(A.data() != B.data() && A(1,1) == B(1,1) && "KinematicTensor copy failed.");

    MaterialTensor T;
    T = S;
    assert(S.data() != T.data() && T(0,2) == S(0,2) && "MaterialTensor copy failed.");


    /*------------------------------------------------------------------------*/
    //test negative
    std::cout << "Checking negative of tensors." << std::endl;

    assert((-B)(1,1) == -(B(1,1)) && "KinematicTensor negative failed.");

    assert((-S)(1,1) == -(S(1,1)) && "MaterialTensor negative failed.");


    /*------------------------------------------------------------------------*/
    //test scalar *=
    std::cout << "Checking *= of tensors." << std::endl;

    B *= 4; B*=2.2;
    assert(B(1,0) == 4*2.2*A(1,0) && "KinematicTensor scalar multiply failed.");

    T *= 2; T*=1.4;
    assert(T(2,2) == 2*1.4*S(2,2) && "MaterialTensor scalar multiply failed.");


    /*------------------------------------------------------------------------*/
    //test scalar division
    std::cout << "Checking scalar division of tensors." << std::endl;

    KinematicTensor C = B;
    C /= 2; C /= 0.25;
    assert(C(0,0) == B(0,0)/(2*0.25) && "KinematicTensor scalar divide failed.");

    MaterialTensor U = T;
    U /= 3; U /= 0.2;
    assert(U(0,2) == (T(0,2)/3)/0.2 && "MaterialTensor scalar divide failed.");


    /*------------------------------------------------------------------------*/
    //test deviator
    std::cout << "Checking deviator of tensors." << std::endl;

    MaterialTensor C0 = C.deviator();
    assert(C0(0,0) == C(0,0) - 1.0/3.0*(C(0,0) + C(1,1)) && "KinematicTensor deviator failed.");
    MaterialTensor U0 = U.deviator();
    assert(U0(0,0) == U(0,0) - 1.0/3.0*(U(0,0) + U(1,1) + U(2,2)) && "MaterialTensor deviator failed.");


    /*------------------------------------------------------------------------*/
    //test transpose
    std::cout << "Checking transpose of tensors." << std::endl;

    KinematicTensor AT = A.transpose();
    assert(AT(0,1) == A(1,0) && "KinematicTensor transpose failed.");
    MaterialTensor UT = U.transpose();
    assert(UT(2,1) == U(1,2) && "MaterialTensor transpose failed.");


    /*------------------------------------------------------------------------*/
    //test dot products
    std::cout << "Checking contraction of tensors." << std::endl;

    assert(A.dot(U) == U.dot(A) && A.dot(B) == B.dot(A) && T.dot(U) == U.dot(T) && "Tensor contraction failed.");


    /*------------------------------------------------------------------------*/
    //test trace
    std::cout << "Checking trace of tensors." << std::endl;

    assert(A.trace() == A(0,0)+A(1,1) && T.trace() == T(0,0)+T(1,1)+T(2,2) && "Tensor trace failed.");


    /*------------------------------------------------------------------------*/
    //test addition
    std::cout << "Checking addition of tensors." << std::endl;

    KinematicTensor D = AT + A + C + B;
    assert(D(0,0) == AT(0,0)+A(0,0)+C(0,0)+B(0,0) && "KinematicTensor addition failed.");
    MaterialTensor W = U + T + D + UT;
    assert(W(0,1) == U(0,1) + T(0,1) + D(0,1) + UT(0,1) && "MaterialTensor addition failed.");


    /*------------------------------------------------------------------------*/
    //test subtraction
    std::cout << "Checking subtraction of tensors." << std::endl;

    KinematicTensor E = -D - A;
    assert(E(0,1) == -(D(0,1)) - A(0,1) && "KinematicTensor subtraction failed.");
    MaterialTensor X = -W - T;
    assert(X(2,2) == -(W(2,2)) - T(2,2) && "MaterialTensor subtraction failed.");


    /*------------------------------------------------------------------------*/
    //test multiplication
    std::cout << "Checking multiplication of tensors." << std::endl;

    KinematicTensor F = D*A + B*2.0 + 3*C;
    assert(F(0,0) == D(0,1)*A(1,0) + D(0,0)*A(0,0) + B(0,0)*2.0 + 3*C(0,0) && "KinematicTensor multiplication fialed.");
    MaterialTensor Y = X*W + X*D + 3.0*U;
    assert(Y(0,1) == X(0,0)*W(0,1) + X(0,1)*W(1,1) + X(0,2)*W(2,1) + X(0,0)*D(0,1) + X(0,1)*D(1,1) + 3*U(0,1) && "MaterialTensor multiplication failed.");


    /*------------------------------------------------------------------------*/
    //test division
    std::cout << "Checking division of tensors." << std::endl;

    KinematicTensor G = D/0.2 + F/3;
    assert(G(1,1) == D(1,1)/0.2 + F(1,1)/3 && "KinematicTensor division failed.");
    MaterialTensor Z = U/1.4 + T/2;
    assert(Z(2,1) == (U(2,1)/1.4 + T(2,1)/2) && "MaterialTensor division failed.");


    /*------------------------------------------------------------------------*/
    //test mapping and array
    std::cout << "Checking mapping of tensors." << std::endl;

    KinematicTensorArray a_array = KinematicTensorArray(2,tensor_type);
    a_array[0] = A; a_array[1] = B; a_array.push_back(D);
    assert(a_array(0,1,1) == A(1,1) && a_array(1,0,1) == B(0,1) && a_array(2,7) == D(7) && "KinematicTensorMap and KinematicTensorArray failed.");

    MaterialTensorArray b_array = MaterialTensorArray(3);
    b_array[0] = T; b_array[1] = U; b_array.push_back(W);
    assert(b_array(0,2,1) == T(2,1) && b_array(1,0,0) == U(0,0) && b_array(3,7) == W(7) && "MaterialTensorMap and MaterialTensorArray failed.");


    /*------------------------------------------------------------------------*/
    //test eigen mapping
    std::cout << "Checking Eigen::Map of tensors." << std::endl;

    EIGEN_MAP_OF_KINEMATIC_TENSOR EB = EIGEN_MAP_OF_KINEMATIC_TENSOR(B);
    assert(EB(1,1) == B(1,1) && "KinematicTensor to Eigen::Map failed.");

    EIGEN_MAP_OF_MATERIAL_TENSOR ET = EIGEN_MAP_OF_MATERIAL_TENSOR(T);
    assert(ET(2,1) == T(2,1) && "KinematicTensor to Eigen::Map failed.");




    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    /*-------------------------Vector Algebra---------------------------------*/
    /*------------------------------------------------------------------------*/
    /*------------------------------------------------------------------------*/
    //choose dimension of kinematic vectors
    int vector_type = KinematicVector::VECTOR_2D;

    /*------------------------------------------------------------------------*/
    //test initialization of tensors
    std::cout << "Checking initializtion of vectors." << std::endl;

    KinematicVector av = KinematicVector(u.data(),tensor_type);
    assert(av(1) == u[1] && "KinematicVector initialization failed.");

    MaterialVector sv = MaterialVector(v.data());
    assert(sv(2) == v[2] && "MaterialTensor initialization failed.");


    /*------------------------------------------------------------------------*/
    //test copying of tensors
    std::cout << "Checking copying of vectors." << std::endl;

    KinematicVector bv;
    bv = av;
    assert(av.data() != bv.data() && av(1) == bv(1) && "KinematicVector copy failed.");

    MaterialVector tv;
    tv = sv;
    assert(sv.data() != tv.data() && sv(0) == sv(0) && "MaterialVector copy failed.");


    /*------------------------------------------------------------------------*/
    //test negative
    std::cout << "Checking negative of vectors." << std::endl;

    assert((-bv)(1) == -(bv(1)) && "KinematicVector negative failed.");

    assert((-sv)(1) == -(sv(1)) && "MaterialVector negative failed.");


    /*------------------------------------------------------------------------*/
    //test scalar *=
    std::cout << "Checking *= of vectors." << std::endl;

    bv *= 4; bv*=2.2;
    assert(bv(1) == 4*2.2*av(1) && "KinematicVector scalar multiply failed.");

    tv *= 2; tv*=1.4;
    assert(tv(2) == 2*1.4*sv(2) && "MaterialVector scalar multiply failed.");


    /*------------------------------------------------------------------------*/
    //test scalar division
    std::cout << "Checking scalar division of vectors." << std::endl;

    KinematicVector cv = bv;
    cv /= 2; cv /= 0.25;
    assert(cv(0) == bv(0)/(2*0.25) && "KinematicVector scalar divide failed.");

    MaterialVector uv = tv;
    uv /= 3; uv /= 0.2;
    assert(uv(2) == (tv(2)/3)/0.2 && "MaterialVector scalar divide failed.");


    /*------------------------------------------------------------------------*/
    //test cross product
    std::cout << "Checking cross product of vectors." << std::endl;

    MaterialVector cxb = cv.cross(bv);
    assert(cxb(2) == cv(0)*bv(1) - cv(1)*bv(0) && "KinematicVector cross product failed.");
    MaterialVector uxt = uv.cross(tv);
    assert(uxt(0) == uv(1)*tv(2) - uv(2)*tv(1) && "MaterialVector cross product failed.");


    /*------------------------------------------------------------------------*/
    //test tensor product
    std::cout << "Checking tensor product of vectors." << std::endl;

    KinematicTensor CXB = cv.tensor(bv);
    assert(CXB(0,1) == cv(0)*bv(1) && "KinematicVector tensor product failed.");
    MaterialTensor UXT = uv.tensor(tv);
    MaterialTensor UXB = uv.tensor(bv);
    MaterialTensor BXT = bv.tensor(tv);
    assert(UXT(1,1) == uv(1)*tv(1) && UXB(1,0) == uv(1)*bv(0) && BXT(1,1) == bv(1)*tv(1) && "MaterialVector tensor product failed.");


    /*------------------------------------------------------------------------*/
    //test dot products
    std::cout << "Checking contraction of vectors." << std::endl;

    assert(av.dot(uv) == uv.dot(av) && av.dot(bv) == bv.dot(av) && tv.dot(uv) == uv.dot(tv) && "Vector contraction failed.");

    /*------------------------------------------------------------------------*/
    //test addition
    std::cout << "Checking addition of vectors." << std::endl;

    KinematicVector dv = av + cv + bv;
    assert(dv(0) == av(0) + cv(0) + bv(0) && "KinematicVector addition failed.");
    MaterialVector wv = uv + tv + dv;
    assert(wv(1) == uv(1) + tv(1) + dv(1) && "MaterialVector addition failed.");


    /*------------------------------------------------------------------------*/
    //test subtraction
    std::cout << "Checking subtraction of vectors." << std::endl;

    KinematicVector ev = -dv - av;
    assert(ev(0) == -(dv(0)) - av(0) && "KinematicVector subtraction failed.");
    MaterialVector xv = -wv - tv;
    assert(xv(2) == -(wv(2)) - tv(2) && "MaterialVector subtraction failed.");


    /*------------------------------------------------------------------------*/
    //test multiplication
    std::cout << "Checking tensor multiplication of vectors." << std::endl;

    KinematicVector fv = D*av + bv*2.0 + 3*cv;
    assert(fv(0) == D(0,0)*av(0) + D(0,1)*av(1) + bv(0)*2.0 + 3*cv(0) && "KinematicVector multiplication fialed.");
    MaterialVector yv = T*tv + S*av + 3.0*uv;
    assert(yv(0) == T(0,0)*tv(0) + T(0,1)*tv(1) + T(0,2)*tv(2) + S(0,0)*av(0) + S(0,1)*av(1) + 3.0*uv(0) && "MaterialVector multiplication failed.");


    /*------------------------------------------------------------------------*/
    //test division
    std::cout << "Checking division of vectors." << std::endl;

    KinematicVector gv = dv/0.2 + fv/3;
    assert(gv(1) == dv(1)/0.2 + fv(1)/3 && "KinematicVector division failed.");
    MaterialVector zv = uv/1.4 + tv/2;
    assert(zv(2) == (uv(2)/1.4 + tv(2)/2) && "MaterialVector division failed.");


    /*------------------------------------------------------------------------*/
    //test mapping and array
    std::cout << "Checking mapping of vectors." << std::endl;

    KinematicVectorArray av_array = KinematicVectorArray(2,vector_type);
    av_array[0] = av; av_array[1] = bv; av_array.push_back(dv);
    av_array(2,1) = 4;
    assert(av_array(0,1) == av(1) && av_array(1,0) == bv(0) && av_array(2,1) == 4 && "KinematicVector::Map and KinematicVectorArray failed.");

    MaterialVectorArray bv_array = MaterialVectorArray(3);
    bv_array[0] = tv; bv_array[1] = uv; bv_array.push_back(wv);
    bv_array(3,1) = W(7);
    assert(bv_array(0,2) == tv(2) && bv_array(1,0) == uv(0) && bv_array(3,1) == W(7) && "MaterialVector::Map and MaterialVectorArray failed.");


    /*------------------------------------------------------------------------*/
    //test eigen mapping
    std::cout << "Checking Eigen::Map of vectors." << std::endl;

    EIGEN_MAP_OF_KINEMATIC_VECTOR eb = EIGEN_MAP_OF_KINEMATIC_VECTOR(bv);
    assert(eb(1) == bv(1) && "KinematicVector to Eigen::Map failed.");

    EIGEN_MAP_OF_MATERIAL_VECTOR et = EIGEN_MAP_OF_MATERIAL_VECTOR(tv);
    assert(et(2) == tv(2) && "KinematicVEctor to Eigen::Map failed.");

    std::cout << "Passed all algebra tests." << std::endl;
    return;
}
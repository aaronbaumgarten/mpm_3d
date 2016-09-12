//
// Created by aaron on 8/26/16.
// spmd.hpp
//

#ifndef MPM_3D_SPMD_HPP
#define MPM_3D_SPMD_HPP

class sparsematrix_double {
public:
    double *vals;
    size_t nnz;
    size_t capacity;

    size_t *column_index;
    size_t *row_pointer;
    size_t rows;
    size_t columns;

    //functions
    sparsematrix_double(size_t, size_t, size_t);
    ~sparsematrix_double();
    void spmd_print(const sparsematrix_double*, int);

    sparsematrix_double * spmd_transpose(sparsematrix_double*);

    void spmd_gaxpy(const sparsematrix_double*, const double*, double*);
    void spmd_gatxpy(const sparsematrix_double*, const double*, double*);

    void spmdv(double, const sparsematrix_double*, const double);

    double spmd_slow_get(const sparsematrix_double*, size_t, size_t);
};

#endif //MPM_3D_SPMD_HPP

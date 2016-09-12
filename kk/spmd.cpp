//
// Created by aaron on 8/26/16.
// spmd.cpp (stolen from mpm-2d-legacy)
//

#include <stdlib.h>
#include <stdio.h>

#include "spmd.hpp"

sparsematrix_double::sparsematrix_double(size_t nnzIn, size_t rowsIn, size_t columnsIn):
        nnz(nnzIn),
        capacity(nnzIn),
        rows(rowsIn),
        columns(columnsIn),
        vals(new double[nnzIn]),
        column_index(new size_t[nnzIn]),
        row_pointer(new size_t[rowsIn+1])
{
    row_pointer[rowsIn] = nnzIn;
}

sparsematrix_double::~sparsematrix_double(){
    delete(vals);
    delete(column_index);
    delete(row_pointer);
}

sparsematrix_double * sparsematrix_double::spmd_transpose(sparsematrix_double *sp)
{
    sparsematrix_double *spt  = new sparsematrix_double(sp->nnz, sp->columns, sp->rows);

    for (size_t i = 0; i < spt->rows; i++) {
        spt->row_pointer[i] = 0;
    }

    // count number of elements **in previous row** in transpose matrix
    for (size_t i = 0; i < sp->nnz; i++) {
        spt->row_pointer[1 + sp->column_index[i]]++;
    }

    for (size_t i = 1; i < spt->rows; i++) {
        spt->row_pointer[i] += spt->row_pointer[i - 1];
    }

    for (size_t i = 0; i < rows; i++) {
        for (size_t p = sp->row_pointer[i]; p < sp->row_pointer[i+1]; p++) {
            size_t transpose_row = sp->column_index[p];
            size_t transpose_linear_index = spt->row_pointer[transpose_row];
            spt->vals[transpose_linear_index] = sp->vals[p];
            spt->column_index[transpose_linear_index] = i;
            spt->row_pointer[transpose_row]++;
        }
    }

    for (size_t i = spt->rows; i > 0; i--) {
        spt->row_pointer[i] = spt->row_pointer[i - 1];
    }

    spt->row_pointer[0] = 0;
    spt->row_pointer[spt->rows] = spt->nnz;

    return spt;
}

void sparsematrix_double::spmd_print(const sparsematrix_double *sp, int full)
{
    if (sp == NULL) {
        return;
    }

    printf("nnz: %zu, rows: %zu, cols: %zu\n", sp->nnz, sp->rows, sp->columns);
    if (full == 0) {
        printf("vals: ");
        for (size_t i = 0; i < sp->nnz; i++) {
            printf("%g ", sp->vals[i]);
        }
        printf("\ncolumn_index: ");
        for (size_t i = 0; i < sp->nnz; i++) {
            printf("%zu ", sp->column_index[i]);
        }
        printf("\nrow_pointer: ");
        for (size_t i = 0; i < sp->rows; i++) {
            printf("%zu ", sp->row_pointer[i]);
        }
        printf("\n");
    } else {
        for (size_t i = 0; i < sp->rows; i++) {
            for (size_t j = 0; j < sp->columns; j++) {
                printf("%4.4g ", sparsematrix_double::spmd_slow_get(sp, i, j));
            }
            printf("\n");
        }
    }

    return;
}

double sparsematrix_double::spmd_slow_get(const sparsematrix_double *sp, size_t i, size_t j)
{
    if (sp == NULL) {
        return 0xdeadbeef;
    }

    double ret = 0;
    for (size_t p = sp->row_pointer[i]; p < sp->row_pointer[i+1]; p++) {
        if (sp->column_index[p] == j) {
            ret += sp->vals[p];
        }
    }

    return ret;
}

void sparsematrix_double::spmd_gaxpy(const sparsematrix_double * A, const double * x, double * y)
{
    // does y <- A * x + y
    if (A == NULL || x == NULL || y == NULL) {
        return;
    }

    for (size_t i = 0; i < A->rows; i++) {
        for (size_t p = A->row_pointer[i]; p < A->row_pointer[i+1]; p++) {
            size_t j = A->column_index[p];
            y[i] += A->vals[p] * x[j];
        }
    }

    return;
}

void sparsematrix_double::spmd_gatxpy(const sparsematrix_double * A, const double * x, double * y)
{
    // does y <- A^T * x + y
    if (A == NULL || x == NULL || y == NULL) {
        return;
    }

    for (size_t i = 0; i < A->rows; i++) {
        for (size_t p = A->row_pointer[i]; p < A->row_pointer[i+1]; p++) {
            size_t j = A->column_index[p];
            y[j] += A->vals[p] * x[i];
        }
    }

    return;
}

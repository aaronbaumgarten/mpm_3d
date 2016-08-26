//
// Created by aaron on 8/26/16.
// spmd.cpp (stolen from mpm-2d-legacy)
//

#include <stdlib.h>
#include <stdio.h>

#include "spmd.hpp"

double spmd_slow_get(const struct sparsematrix_double *sp, size_t i, size_t j);

struct sparsematrix_double *spmd_create(size_t nnz, size_t rows, size_t columns)
{
    struct sparsematrix_double *sp = malloc(sizeof(struct sparsematrix_double));
    if (sp == NULL) {
        return NULL;
    }
    sp->nnz = nnz;
    sp->capacity = nnz;
    sp->rows = rows;
    sp->columns = columns;
    sp->vals = calloc(nnz, sizeof(double));
    if (sp->vals == NULL) {
        free(sp);
        return NULL;
    }
    sp->column_index = calloc(nnz, sizeof(size_t));
    if (sp->column_index == NULL) {
        free(sp->vals);
        free(sp);
        return NULL;
    }

    sp->row_pointer = calloc(rows + 1, sizeof(size_t));
    if (sp->row_pointer == NULL) {
        free(sp->column_index);
        free(sp->vals);
        free(sp);
        return NULL;
    }
    sp->row_pointer[rows] = nnz;

    return sp;
}

struct sparsematrix_double *spmd_transpose(struct sparsematrix_double *sp)
{
    struct sparsematrix_double *spt = spmd_create(sp->nnz, sp->columns, sp->rows);

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

    for (size_t i = 0; i < sp->rows; i++) {
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

void spmd_delete(struct sparsematrix_double *sp)
{
    if (sp == NULL) {
        return;
    }

    free(sp->row_pointer);
    free(sp->column_index);
    free(sp->vals);
    free(sp);
    return;
}

void spmd_print(const struct sparsematrix_double *sp, int full)
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
                printf("%4.4g ", spmd_slow_get(sp, i, j));
            }
            printf("\n");
        }
    }

    return;
}

double spmd_slow_get(const struct sparsematrix_double *sp, size_t i, size_t j)
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

void spmd_gaxpy(const struct sparsematrix_double * restrict A, const double * restrict x, double * restrict y)
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

void spmd_gatxpy(const struct sparsematrix_double * restrict A, const double * restrict x, double * restrict y)
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

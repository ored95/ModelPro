#ifndef __ALLOC_H__
#define __ALLOC_H__
#include <stdio.h>
#include <stdlib.h>

double* malloc_array(const size_t);

double** malloc_matrix(const size_t, const size_t);

void free_m(size_t, double **);

void display(const size_t, const double *, const double *);

void display_a(const size_t, const double *);

void display_m(const size_t, const size_t, double **);

#endif
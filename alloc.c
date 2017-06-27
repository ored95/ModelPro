#include "alloc.h"

double* 
malloc_array(const size_t n)
{
	double* temp = (double*)malloc(n * sizeof(double));
	if (!temp)
		return NULL;
	return temp;
}

double** 
malloc_matrix(const size_t m, const size_t n)
{
	double **temp = (double **)malloc(m * sizeof(double*));
	if (!temp)
	{
		printf("Error #1: Not enough memory!\n");
		return NULL;
	}
	for (size_t i = 0; i < m; i++)
	{
		temp[i] = malloc_array(n);
	}
	return temp;
}

void 
free_m(size_t n, double **x)
{
	for (size_t i = 0; i < n; i++)
		free(x[i]);
	free(x);
}

void
display(const size_t n, const double *x, const double *y)
{
	printf("\n====== Source ======\n\n");
	for (int i = 0; i < n; i++)
	{
		printf("\t%.2f\t%.2f\n", x[i] ,y[i]);
	}
	printf("\n");
}

void 
display_a(const size_t n, const double *arr)
{
	printf("\nShow array:\n");
	for (int i = 0; i < n; i++)
	{
		printf("%.9f\n", arr[i]);
	}
	printf("\n");
}

void 
display_m(const size_t m, const size_t n, double** arr)
{
	printf("Show matrix:\n\n");
	for (int i = 0; i < m; i++)
	{
		for (int j = 0 ; j < n; j++)
			printf("%.2f\t", arr[i][j]);
		printf("\n\n");
	}
}
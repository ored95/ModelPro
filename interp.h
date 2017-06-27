#ifndef __INTERP_H__
#define __INTERP_H__
#include <stdio.h>
#include <math.h>
/* Update version */

void show_config(const int, const int, const double *, const double *);

void placed(int *, int *, double, const double *, const size_t);

double quick_interpol_linear(double, const double, const double, const double, const double);

double quick_interpol_logarithm(double, const double, const double, const double, const double);

double interp_log(double t_value, const size_t, const double *, const double *);

#endif
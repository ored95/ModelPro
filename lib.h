#ifndef __LIB_H__
#define __LIB_H__
#include <math.h>
#include "alloc.h"

void show_config(const int, const int, const double *, const double *);

/* Update version */
void placed(int *, int *, double, const double *, const size_t);

double quick_interpol_linear(double, const double, const double, const double, const double);

double interp(double t_value, const size_t, const double *, const double *);

double quick_interpol_logarithm(double, const double, const double, const double, const double);

double interp_log(double t_value, const size_t, const double *, const double *);

#define Tw 2000		/* temperature (Kelvin) */
double temper(double, const double, const double);

double* getKi(const size_t, const char *, const double, const double);

#define PI 3.14159265359		/*	[]		*/
#define Hp 6.6260704E-34		/*	[J.s]	*/
#define Ch 3.0E+10				/* 	[cm/s] 	*/
#define Kh 1.3806485279E-23 	/* 	[J/K] 	*/
#define left_freq 0.02E+15		/* 	[Hz]	*/
#define right_freq 3.0E+15		/* 	[Hz]	*/
#define Ra	0.35				/*	[cm]	*/
#define Df (right_freq - left_freq)
#define Mf ((left_freq + right_freq) / 2.0)

double Up(double);

void get_co(const double, const double, const double *, const size_t, 
			double **, double **, double **, double **);

void show_tridiagonal(const size_t, const double *, const double *, const double *, const double *);

double integ(const int, const double, const double, const int, const double *, const double *);

void solve_tridiagonal(const char *, const double, const double, const size_t,
					   double **, double **, double **);

double Ip(double, double, const double);

double interp_k(double, const size_t, const double *);

double phi(double, double, double, const double, const size_t, const double *);

double intensity(double, double, double, const double, const double, const size_t, const double *);

double u_help(double, double, const double, const double, const size_t, const double *);

double Ui(double, const double, const double, const size_t, const double *);

double* save_Ui(const double, const double, const size_t, const double *);

void save_solution(const char *, const double, const double, const size_t);

double Upf(double, double, const double, const double);

double* get_q(const double, const double, const size_t);

#endif
#ifndef __SCHEME_H__
#define __SCHEME_H__
#include "alloc.h"
#include "interp.h"
#include <math.h>

#define Tw 2000		/* temperature (Kelvin) */
#define PI 3.14159265359		/*	[]		*/
#define Hp 6.6260704E-34		/*	[J.s]	*/
#define Ch 3.0E+10				/* 	[cm/s] 	*/
#define Kh 1.3806485279E-23 	/* 	[J/K] 	*/
//#define Ra 0.35					/*	[cm]	*/
#define C1 (8. * PI * Hp * 1.E+30)
#define C2 (Hp * 1.E+15 / Kh)
#define mN (2. / 0.848 / 3.)	/* Old version */
//#define mN (2. / 0.86065 / 3.)	/* New version */
//#define mN 0.75

double temper(double, const double, const double);

double Up(double, double, double);

void get_co(double, double, const double, const double, const double,
			const size_t, const double *, const double *,
			const size_t, double **, double **, double **, double **);

void show_tridiagonal(const size_t, const double *, const double *, const double *, const double *);

void solve_tridiagonal(double, double, const double, const double, const double,
					   const size_t, const double *, const double *,
					   const size_t, double **);
					   
void getq(const double, const double, const size_t, const double, const double,
		  double **, double **, double **);
					   
void save_solution(const double, const double, const size_t);

void save_Qrelation(const size_t);

#endif
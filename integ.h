#ifndef __INTEG_H__
#define __INTEG_H__
#include <math.h>
#include "alloc.h"
#include "scheme.h"
#include <stdio.h>
#define __C1 (2. * Hp * 1.E+40)
#define __C0 (24. * PI * Hp * 1.E+40)

double __Ip(double, double, double);
double __Up(double, double, double);

void save_F(const double, const double);

#endif
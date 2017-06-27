#include "scheme.h"
#include <math.h>

int main()
{	
	double T0 = 10000, m = 4.;
	int N;
	printf("+ Number of intervals: ");
	scanf("%d", &N);
	
	save_solution(T0, m, N);
	
	// //double ki = 3.56E-1;
	// double freq = (0.19463 + 0.19500) / 2.;
	// double dfreq = 0.19500 - 0.19463;
	// double tf = temper(0., T0, m);
	
	// double upi = Up(tf, freq, dfreq);
	
	// double u = 8. * PI * Hp * pow(freq, 3.) * dfreq * 1.E-45 / pow(Ch, 3.) / (exp(Hp * freq * 1.E+15 / Kh * tf) - 1.);
	
	// printf("\n\n%E\t%E\n", upi, u);
	return 0;
}
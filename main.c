#include "scheme.h"
//#include <math.h>
#include "integ.h"

int main()
{	
	double T0 = 10000, m = 8.;
	//int N = 80;
	//printf("+ Number of intervals: ");
	//scanf("%d", &N);
	
	//save_solution(T0, m, N);
	//save_Qrelation(N);
	save_F(T0, m);
	
	return 0;
}
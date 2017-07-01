#include "scheme.h"
#include <math.h>

int main()
{	
	double T0 = 10000, m = 4.;
	int N;
	printf("+ Number of intervals: ");
	scanf("%d", &N);
	
	save_solution(T0, m, N);
	//save_Qrelation(N);
	
	return 0;
}
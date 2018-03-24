#include "scheme.h"
//#include <math.h>
#include "integ.h"
#include "stream.h"

int main()
{	
	//double T0 = 8000, m = 8.;
	//int N = 80;
	//printf("+ Number of intervals: ");
	//scanf("%d", &N);
	
	//save_solution(T0, m, 50);// printf("\n======================== ");
	//save_Qrelation(N);
	//save_F(T0, m);
	printf("Stream: %E\n", get_stream(50));
	return 0;
}
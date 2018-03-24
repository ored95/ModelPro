#include "stream.h"
#include <stdio.h>

double get_stream(const int __nZ)
{
	FILE *fs = fopen("F.txt", "r");
	if (!fs)
		return 0;
	
	double *ref = malloc_array(__nZ + 1);
	for (int i = 0; i <= __nZ; i++)
		fscanf(fs, "%lf", &ref[i]);
	fclose(fs);
	
	#define Radius 0.35
	double stream = 0.;
	double dz2 = 1. / __nZ / __nZ;
	for (int i = 0; i <= __nZ; i++)
	{
		stream += i * ref[i];
	}
	stream *= dz2 * Radius;
	
	#undef Radius
	
	return stream;
}
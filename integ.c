#include "integ.h"

double 
__Ip(double tf, double freq)
{
	double v0 = __C1 * freq * pow((freq / 3.), 2.);
	double v1 = exp(C2 * freq / tf) - 1.;
	
	return v0 / v1;
}

double
__Up(double tf, double freq)
{
	double v0 = __C0 * pow((freq / 3.), 3.);
	double v1 = exp(C2 * freq / tf) - 1.;
	
	return v0 / v1;
}

void 
save_F(const double T0, const double m)
{
	/* ========== READ DATA =========== */
	double *frequency = malloc_array(194);
	double *temperature = malloc_array(16);
	double **KT = malloc_matrix(193, 16);
	
	FILE *fs = fopen("./sources/frequency.txt", "r");
	rewind(fs);
	for (int i = 0; i < 194; i++)
	{
		fscanf(fs, "%lf", &frequency[i]);
	}
	fclose(fs);
	
	fs = fopen("./sources/temperature.txt", "r");
	rewind(fs);
	for (int i = 0; i < 16; i++)
	{
		fscanf(fs, "%lf", &temperature[i]);
	}
	fclose(fs);
	
	fs = fopen("./sources/TK.txt", "r");
	rewind(fs);
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 193; j++)
		{
			fscanf(fs, "%lf", &KT[j][i]);
		}
	}
	fclose(fs);
	
	/* ===========    Split K(T, F) array    ============ */
	
	/* =========== Solution of integral form ============ */
	size_t __nPhi = 100, __nTheta = 50, __nZ = 80, __nKj = 100;
	double **U = malloc_matrix(193, __nZ + 1);
	
	double dPhi = PI / __nPhi,
		   dTheta = PI / 2./ __nTheta;
	double __sPhi = 0.,
		   __sTheta = dTheta / 2.;

	double dz = 1. / __nZ;
	double *kj = malloc_array(__nKj + 1);	// save kj in the main step
	double *IP = malloc_array(__nKj + 1);	// save IP
	double result = 0.; 					// save the main stuff

#define Radius 0.35		/* cm */	
	for (int iz = 0; iz <= __nZ; iz++)
	{
		double z = dz * iz;
		//	   z00 = z * z;				// position
		//double tf = temper(z, T0, m);	// temperature (changed)
		
		for (int k = 0; k < 193; k++)
		{
			U[k][iz] = 0.;		// init
			double freq = (frequency[k+1] + frequency[k]) / 2.;		// frequency (const)

			// k in KT[k]
			for (int i = 0; i <= __nPhi; i++)
			{
				double phi = __sPhi + dPhi * i;
				for (int j = 0; j < __nTheta; j++)
				{
					double theta = __sTheta + dTheta * j;
				
					// Now calculate main integral by method average-rectangle
					double z1 = z * cos(phi),
						   //z11 = z1 * z1,
						   z22 = 1. - z * z * sin(phi) * sin(phi),
						   z2 = sqrt(z22),
						   d = Radius * (z1 + z2);

					// Calculate integral of (k * Ip * exp(-|kdy) * d(x)
					result = 0.;
					
					// But first, let prepare all coefficient k
					int nKj;				// just index
					double d1, t1;			// relation by segment and temperature

					for (nKj = 0; nKj <= __nKj; nKj++)
					{
						d1 = d * nKj / __nKj;
						t1 = temper(d1 / Radius, T0, m);
						IP[nKj] = __Ip(t1, freq);
						kj[nKj] = interp_log(t1, 16, temperature, KT[k]);
					}
					
					// And now calculate our main stuff
					for (int nKj = 0; nKj <= __nKj; nKj++)
					{
						double s = 0.;
						for (int jx = nKj + 1; jx <= __nKj; jx++)
							s += kj[jx];

						s = s * d / __nKj / sin(theta);
						// Result
						result += kj[nKj] * IP[nKj] * exp(-s);
					}

					U[k][iz] += result * d / sin(theta);
				}
			}

			U[k][iz] *= dTheta * dPhi; printf("\rdone: %d : %d", iz, k);
		}
	}

	free(kj);
	free(IP);

	// Complete our problem
	
	double *div = malloc_array(__nZ + 1);
	printf("\n================================\n\n");
	for (int iz = 0; iz <= __nZ; iz++)
	{
		div[iz] = 0.;
		double z = dz * iz;
		double tf = temper(z, T0, m);	// temperature (changed)
		
		for (int k = 0; k < 193; k++)
		{
			double ki = interp_log(tf, 16, temperature, KT[k]);
			double freq = (frequency[k] + frequency[k+1]) / 2.;
			
			double UP = __Up(tf, freq);
			div[iz] += ki * (UP - U[k][iz]);
		}

		printf("%.5f\n", div[iz]);
	}
	free(div);
	free_m(193, U);
	free(temperature); free(frequency); free_m(16, KT);
}
#include "integ.h"

double 
__Ip(double tf, double freq, double dfreq)
{
	double v0 = __C1 * dfreq * freq * pow((freq / 3.), 2.);
	double v1 = exp(C2 * freq / tf) - 1.;
	
	return v0 / v1;
}

double
__Up(double tf, double freq, double dfreq)
{
	double v0 = __C0 * dfreq * pow((freq / 3.), 3.);
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
	
	/* =========== Solution of integral form ============ */
	size_t __nPhi = 40, __nTheta = 20, __nZ = 40, __nKj = 100;
	double **U = malloc_matrix(193, __nZ + 1);
	
	double dPhi = PI / __nPhi,
		   dTheta = PI / 2./ __nTheta;
	double __sPhi = dPhi / 2.,
		   __sTheta = dTheta / 2.;

	double dz = 1. / __nZ;
	double *kj = malloc_array(__nKj + 1);	// save kj in the main step
	double *IP = malloc_array(__nKj + 1);	// save IP
	double result = 0.; 					// save the main stuff

#define Radius 0.35		/* cm */	
	for (int iz = 0; iz <= __nZ; iz++)
	{
		double z = dz * iz;				// position (* Radius)
		//double tf = temper(z, T0, m);	// temperature (changed)
		
		for (int k = 0; k < 193; k++)
		{
			U[k][iz] = 0.;		// init
			double freq = (frequency[k+1] + frequency[k]) / 2.;		// frequency (const)
			double dfreq = frequency[k+1] - frequency[k];
			
			// k in KT[k]
			for (int i = 0; i < __nPhi; i++)
			{
				double phi = __sPhi + dPhi * i;
				for (int j = 0; j < __nTheta; j++)
				{
					double theta = __sTheta + dTheta * j;
				
					// Now calculate main integral by method average-rectangle
					double D = z * cos(phi) + sqrt(1. - z * z * sin(phi) * sin(phi));
					
					// Calculate integral of (k * Ip * exp(-|kdy) * d(x)
					result = 0.;
					
					// But first, let prepare all coefficient k
					int nKj;				// just index
					double k1, zt, t1;		// relation by segment and temperature

					for (nKj = 0; nKj <= __nKj; nKj++)
					{
						k1 = 1. - nKj / __nKj;
						zt = sqrt(z * z + k1 * k1 * D * D - 2. * k1 * D  * z * cos(phi));
						//printf("%.2f : %.9f : %.9f : %.9f\n", k1, z, D, zt);
						t1 = temper(zt, T0, m);
						IP[nKj] = __Ip(t1, freq, dfreq);
						kj[nKj] = interp_log(t1, 16, temperature, KT[k]);
					}
					
					// And now calculate our main stuff
					for (int nKj = 0; nKj <= __nKj; nKj++)
					{
						double s = 0.;
						for (int jx = nKj + 1; jx <= __nKj; jx++)
						{
							s += kj[jx];
							if (s * Radius / __nZ > 1.)
							{
								s -= kj[jx];
								break;
							}
						}

						s *= (D * Radius / sin(theta)) / __nKj;
						// Result
						result += kj[nKj] * IP[nKj] * exp(-s);
					}
					
					//U[k][iz] += result * (D * Radius * sin(theta)) / __nKj / sin(theta);
					U[k][iz] += result * D;
				}
			}

			U[k][iz] *= 4. * dTheta * dPhi * (Radius / __nKj); printf("\r%3d : %3d", iz, k);
		} //if (iz == 2) break;
	}
	//printf("\n\n");
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
		
		int k = 0;
		double s = 0.;
		for (; k < 193; k++)
		{
			double ki = interp_log(tf, 16, temperature, KT[k]);
			s += ki;

			if (s * Radius / __nZ > 1.)
				break;
		}

		for (int j = 0; j < k; j++)
		{
			double ki = interp_log(tf, 16, temperature, KT[j]);
			double freq = (frequency[j] + frequency[j+1]) / 2.;
			double dfreq = frequency[j+1] - frequency[j];
			
			double UP = __Up(tf, freq, dfreq);
			div[iz] += ki * (UP - U[j][iz]);

			//printf("%3d:\t%.5f\t%E\t%E\t%E\n", k, ki, UP, U[k][iz], div[iz]);
		}
		
		printf("%E\n", div[iz]);
	}

	// fs = fopen("res03.xls", "w");
	// for (int iz = 0; iz <= __nZ; iz++)
	// {
	// 	double z = dz * iz;
	// 	fprintf(fs, "%.5f\t%E\n", z, div[iz]);
	// }
	// fclose(fs);

	free(div);
	free_m(193, U);
	free(temperature); free(frequency); free_m(16, KT);
#undef Radius
}
#include "scheme.h"

double
temper(double x, const double T0, const double m)
{
	return T0 + (Tw - T0) * pow(x, m);
}

double 
Up(double tf, double freq, double dfreq)
{
	double v0 = C1 * pow((freq / 3.), 3.) * dfreq;
	double v1 = exp(C2 * freq / tf) - 1.;
	
	return v0 / v1;
}

void 
get_co(double freq, double dfreq, const double T0, const double m, const double radius,
	   const size_t Ntk, const double *T, const double *k, 
	   const size_t N, double **A, double **B, double **C, double **D)
{
	*A = malloc_array(N + 1);
	*B = malloc_array(N + 1);
	*C = malloc_array(N + 1);
	*D = malloc_array(N + 1);
	
	double Ra = radius;
	double 	h = 1. / N;					// step
	double r3h = 3. * Ra * Ra * h * h;
	
	double  t00 = temper(         0., T0, m),
			t01 = temper(     h / 2., T0, m),
			tN1 = temper(1. - h / 2., T0, m),
			tN0 = temper(         1., T0, m);
	double 	k00 = interp_log(t00, Ntk, T, k),
			k01 = interp_log(t01, Ntk, T, k),
			kN1 = interp_log(tN1, Ntk, T, k),
			kN0 = interp_log(tN0, Ntk, T, k);
	double 	u00 = Up(t00, freq, dfreq),
			u01 = Up(t01, freq, dfreq),
			uN1 = Up(tN1, freq, dfreq),
			uN0 = Up(tN0, freq, dfreq);
	
	for (int j = 1; j <= N; j++)		/* filling arrays A, C */
	{
		double xf = j * h - h / 2.;
		double tf = temper(xf, T0, m);
		double kf = interp_log(tf, Ntk, T, k);
		
		(*A)[j] = xf / kf;
		(*C)[j-1] = (*A)[j];
	}
	
	for (int j = 1; j < N; j++)			/* filling arrays B, D */
	{
		double xf = j * h;
		double tf = temper(xf, T0, m);
		double kf = interp_log(tf, Ntk, T, k);
		(*B)[j] = (*A)[j] + (*C)[j] + r3h * xf * kf;
		
		double Uf = Up(tf, freq, dfreq);
		(*D)[j] = r3h * xf * kf * Uf;
	}
	
	(*A)[0] = 0.;
	(*C)[0] = 8. / (r3h * k01) - k01 / 2.;
	(*B)[0] = (*A)[0] + (*C)[0] + k00 + k01;
	(*D)[0] = k00 * u00 + k01 * u01;
	
	double alpha = r3h * (2. - h / 2.);
	double beta = Ra * h * (2. - h / 2.);
	(*A)[N] = 8. * (1. - h / 2.) / (alpha * kN1) - kN1 / 2.;
	(*C)[N] = 0.;
	(*B)[N] = (*A)[N] + (*C)[N] + kN0 + kN1 + (4. * mN) / beta;
	(*D)[N] = kN1 * uN1 + kN0 * uN0;
}

void 
show_tridiagonal(const size_t N, const double *A, const double *B, const double *C, const double *D)
{
	printf("\n======= TRIDIAGONAL MATRIX =======\n\n");
	for (int j = 0; j <= N; j++)
	{
		printf("\t%9.5f\t%9.5f\t%9.5f\t%E\n", A[j], B[j], C[j], D[j]);
	}
}

void
solve_tridiagonal(double freq, double dfreq, const double T0, const double m, const double radius,
				  const size_t Ntk, const double *T, const double *k,
				  const size_t N, double **y)
{
	double *A = NULL, *B = NULL, *C = NULL, *D = NULL;
	get_co(freq, dfreq, T0, m, radius, Ntk, T, k, N, &A, &B, &C, &D);
		
	// show_tridiagonal(N, A, B, C, D);		// SHOW MATRIX 
	
	double *p = malloc_array(N + 1),
		   *q = malloc_array(N + 1);

	/// Step #1
	p[1] = C[0] / B[0];
	q[1] = D[0] / B[0];
	
	for (int j = 1; j < N; j++)
	{
		p[j+1] = C[j] / (B[j] - A[j] * p[j]);
		q[j+1] = (D[j] + A[j] * q[j]) / (B[j] - A[j] * p[j]);
	}
	
	/// Step #2
	*y = malloc_array(N + 1);
	(*y)[N] = (A[N] * q[N] + D[N]) / (B[N] - A[N] * p[N]);
	for (int j = N; j > 0; j--)
	{
		(*y)[j - 1] = p[j] * (*y)[j] + q[j];
	}
	
	free(p); free(q);
	free(A); free(B); free(C); free(D);
}

void
get_q(const double T0, const double m, const size_t N, const double radius, const double Tkey,
	  double **Qd, double **Qv, double **T)
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
	
	/* ========== SOLVE SYSTEM U-F ========== */
	double **y = malloc_matrix(193, N+1);	// save Ui
	//double *k_tab = malloc_array(16);		// copies from given data
	double *tmp = NULL;
	for (int i = 0; i < 193; i++)
	{
		/* solve tridiagonal system */	
		double freq = (frequency[i+1] + frequency[i]) / 2.;
		double dfreq = frequency[i+1] - frequency[i];
		
		tmp = NULL;
		solve_tridiagonal(freq, dfreq, T0, m, radius, 16, temperature, KT[i], N, &tmp);
		
		/* save to tab */
		for (int j = 0; j <= N; j++)
			y[i][j] = tmp[j];
		free(tmp);
	}
	
	/* APPLYING BORDER TAU */
	*Qd = malloc_array(N + 1);
	*Qv = malloc_array(N + 1);
	*T = malloc_array(N + 1);
	
	double dx = 1. / N;
#define Ra 0.35 /* cm */
	for (int i = 0; i <= N; i++)
	{
		double xf = i * dx;
		double tf = temper(xf, T0, m);
		
		(*T)[i] = 0.;
		int j = 0;
		
		for (; j < 193; j++)
		{
			double ki = interp_log(tf, 16, temperature, KT[j]);
			(*T)[i] += ki;
			
			if ((*T)[i] * Ra / N > Tkey)
			{
				(*T)[i] -= ki;
				break;
			}
		}
		
		(*Qd)[i] = 0.;
		(*Qv)[i] = 0.;
		(*T)[i] *= Ra / N;
		
		for (int k = 0; k < j; k++)
		{
			double ki = interp_log(tf, 16, temperature, KT[k]);
			double freq = (frequency[k] + frequency[k+1]) / 2.;
			double dfreq = frequency[k+1] - frequency[k];
			
			double upi = Up(tf, freq, dfreq);
			(*Qd)[i] += ki * (upi - y[k][i]);
			(*Qv)[i] += ki * upi;
		}
		
		(*Qd)[i] *= Ch;
		(*Qv)[i] *= Ch;
	}
	
#undef Ra
	
	free_m(193, y);
	free(temperature); free(frequency); free_m(16, KT);
}

void 
save_solution(const double T0, const double m, const size_t N)
{
	double *Qd = NULL, *Qv = NULL, *T = NULL;
	#define Rad 0.35 /* cm */
	get_q(T0, m, N, Rad, 1., &Qd, &Qv, &T);
	
	
	//FILE *fs = fopen("result.xls", "w");
	//fprintf(fs, "x\tQv\tQd\th\tTau\n");
	for (int j = 0; j <= N; j++)
	{
		//double x = (double)j / N;
		//double relate = Qv[j] / Qd[j];
		//fprintf(fs, "%.5f\t%E\t%E\t%E\t%E\n", x, Qv[j], Qd[j], relate, T[j]);
		//fprintf(fs, "%.5f\t%E\n", x, Qd[j]);
		printf("%E\n", Qd[j]);
	}
	#undef Rad
	
	//fclose(fs);
	
	free(Qd);
	free(Qv);
	free(T);
}

void
save_Qrelation(const size_t N)
{
	double *Qd = NULL, *Qv = NULL, *T = NULL;
	double T0 = 10000., m = 4.;
	
	FILE *fs = fopen("Qd-Qv.xls", "w");
	
	for (int i = 0; i <= 8; i++)
		fprintf(fs, "\t%d", (int)T0 + i * 2000);
	fprintf(fs, "\n");
	
	for (int j = 1; j <= 16; j++)
	{
		double radius = 0.1 * j;
		fprintf(fs, "%.5f", radius);
		
		for (int i = 0; i <= 8; i++)
		{
			Qd = NULL; Qv = NULL; T = NULL; 
			double tmp = T0 + i * 2000.;
			
			get_q(tmp, m, N, radius, 1.0, &Qd, &Qv, &T);
				
			double relation = Qv[0] / Qd[0];
			fprintf(fs, "\t%E", relation);
			
			free(Qd); free(Qv);
		}
		
		fprintf(fs, "\n");
	}
	
	fclose(fs);
	free(Qd);
	free(Qv);
	free(T);
}
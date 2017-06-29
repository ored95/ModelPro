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
			kN0 = interp_log(tN1, Ntk, T, k);
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
get_q(const double T0, const double m, const size_t N, const double radius,
	  double **Qd, double **Qv)
{
	/* ========== READ DATA =========== */
	double *frequency = malloc_array(194);
	double *temperature = malloc_array(16);
	double **KT = malloc_matrix(16, 193);
	
	FILE *fs = fopen("frequency.txt", "r");
	rewind(fs);
	for (int i = 0; i < 194; i++)
	{
		fscanf(fs, "%lf", &frequency[i]);
	}
	fclose(fs);
	
	fs = fopen("temperature.txt", "r");
	rewind(fs);
	for (int i = 0; i < 16; i++)
	{
		fscanf(fs, "%lf", &temperature[i]);
	}
	fclose(fs);
	
	fs = fopen("TK.txt", "r");
	rewind(fs);
	for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 193; j++)
		{
			fscanf(fs, "%lf", &KT[i][j]);
		}
	}
	fclose(fs);
	
	/* ========== SOLVE SYSTEM U-F ========== */
	double **y = malloc_matrix(193, N+1);	// save Ui
	double *k_tab = malloc_array(16);		// copies from given data
	double *tmp = NULL;
	for (int i = 0; i < 193; i++)
	{
		/* get copy of value (k) from tab */
		for (int j = 0; j < 16; j++)
			k_tab[j] = KT[j][i];
		
		/* solve tridiagonal system */
		double freq = (frequency[i+1] + frequency[i]) / 2.;
		double dfreq = frequency[i+1] - frequency[i];
		
		tmp = NULL;
		solve_tridiagonal(freq, dfreq, T0, m, radius, 16, temperature, k_tab, N, &tmp);
		
		/* save to tab */
		for (int j = 0;  j <= N; j++)
			y[i][j] = tmp[j];
		free(tmp);
	}
	
	/* ========== CALCULATE Q ========== */
	*Qd = malloc_array(N + 1);
	*Qv = malloc_array(N + 1);
	
	double dx = 1. / N;
		   
	for (int i = 0; i <= N; i++)
	{
		double xf = i * dx;
		double tf = temper(xf, T0, m);
		
		(*Qd)[i] = 0.;
		(*Qv)[i] = 0.;
		
		for (int j = 0; j < 193; j++)
		{
			/* get copy of value (k) from tab */
			for (int jp = 0; jp < 16; jp++)
				k_tab[jp] = KT[jp][j];
			
			double ki = interp_log(tf, 16, temperature, k_tab);
			double freq = (frequency[j] + frequency[j+1]) / 2.;
			double dfreq = frequency[j+1] - frequency[j];
			
			double upi = Up(tf, freq, dfreq);
			//if (j == 22)
				//printf("%E\t%E\t%E\t%E\t%.8f\t%.8f\n", xf, ki, upi, y[j][i], freq, dfreq);
				//printf("%E\t%E\t%E\t%E\n", xf, ki, upi, y[j][i]);
			(*Qd)[i] += ki * (upi - y[j][i]);
			(*Qv)[i] += ki * upi;
		}
		  
		(*Qd)[i] *= Ch;
		(*Qv)[i] *= Ch;
	}

	free(k_tab); free_m(193, y);
	free(temperature); free(frequency); free_m(16, KT);
}

void 
save_solution(const double T0, const double m, const size_t N)
{
	double *Qd = NULL, *Qv = NULL;
	#define Rad 0.35 /* cm */
	get_q(T0, m, N, Rad, &Qd, &Qv);
	#undef Rad
	
	FILE *fs = fopen("result.xls", "w");
	for (int j = 0; j <= N; j++)
	{
		double r = (double)j / N;
		double relate = Qd[j] / Qv[j];
		fprintf(fs, "%.5f\t%E\t%E\t%E\n", r, Qd[j], Qv[j], relate);
	}
	
	fclose(fs);
	free(Qd);
	free(Qv);
}

void
save_Qrelation(const size_t N)
{
	double *Qd = NULL, *Qv = NULL;
	
	double T0 = 10000, m = 4.;
	
	FILE *fs = fopen("Qd-Qv.xls", "w");
	for (int j = 1; j <= 16; j++)
	{
		double radius = 0.1 * j;
		get_q(T0, m, N, radius, &Qd, &Qv);
		
		double relation = Qd[0] / Qv[0];
		fprintf(fs, "%.5f\t%E\t\n", radius, relation);
	}
	
	fclose(fs);
	free(Qd);
	free(Qv);
}
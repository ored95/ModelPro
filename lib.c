#include "lib.h"

void 
show_config(const int IL, const int IR, const double *x, const double *y)
{
    printf("\n===== CONFIGURATION ======\n");
    for (int i = IL; i <= IR; i++)
    {
		printf("%.2f\t%.2f\n", x[i], y[i]);
    }
	printf("\n");
}

void 
placed(int *left, int *right, double value,
	   const double *arr, const size_t n)
{
	*left = *right = 0;
	if (value >= arr[n - 1])
	{
		*left = *right = n - 1;
	}
	else
	{
		for (int i = 0; i < n-1; i++)
		{
			if ((value >= arr[i]) && (value < arr[i+1]))
			{
				*left = i;
				*right = i + 1;
				return;
			}
		}
	}
}

double 
quick_interpol_linear(double x, const double x0, const double x1, const double y0, const double y1)
{
	return ((x1 - x) * y0 + (x - x0) * y1) / (x1 - x0);
}

double 
interp(double t_value,
	   const size_t n, const double *tem, const double *k)
{
	int Lt, Rt;
	
	placed(&Lt, &Rt, t_value, tem, n);
	//show_config(Lt, Rt, tem, k);
	if (Lt == Rt)
		return k[Lt];
	
	return quick_interpol_linear(t_value, tem[Lt], tem[Rt], k[Lt], k[Rt]);
}

/// apply for only positive number
double
quick_interpol_logarithm(double x, const double x0, const double x1, const double y0, const double y1)
{
	double 	_x = log(x),
			_x0 = log(x0), _x1 = log(x1),
			_y0 = log(y0), _y1 = log(y1);
	return exp(quick_interpol_linear(_x, _x0, _x1, _y0, _y1));
}

double 
interp_log(double t_value,
		   const size_t n, const double *tem, const double *k)
{
	int Lt, Rt;
	
	placed(&Lt, &Rt, t_value, tem, n);
	//show_config(Lt, Rt, tem, k);
	if (Lt == Rt)
		return k[Lt];
	
	return quick_interpol_logarithm(t_value, tem[Lt], tem[Rt], k[Lt], k[Rt]);
}

double
temper(double x, const double T0, const double m)
{
	return T0 + (Tw - T0) * pow(x, m);
}

double*
getKi(const size_t N, const char *src, const double T0, const double m)
{
	FILE *fs = fopen(src, "r");
	if (!fs)
		return NULL;
	
	double *tem = malloc_array(14);
	if (!tem)
		return NULL;
	
	double *k = malloc_array(14);
	if (!k)
		return NULL;
	
	rewind(fs);
	
	int i = 0;
	while (!feof(fs))
	{
		fscanf(fs, "%lf%lf", &tem[i], &k[i]);
		i++;
	}
	fclose(fs);
	
	int n0 = 2 * N + 1;
	double *result = malloc_array(n0);
	if (!result)
	{
		free(tem); free(k);
		return NULL;
	}
	
	double h = 0.5 / N;
	
	// get ki when x in [0, 1]
	for (i = 0; i < n0; i++)
	{
		double t = temper(h * i, T0, m);
		result[i] = interp_log(t, 14, tem, k);
		// result[i] = 1.0;
	}
	
	free(tem); free(k);
	return result;
}

double 
Up(double tf)
{
	double v0 = 8. * PI * Hp * pow((Mf / Ch), 3.) * Df;
	double v1 = exp(Hp * Mf / (Kh * tf)) - 1.;
	// printf("T: %.2f\tV0: %E\tV1: %E\n", tf, v0, v1);
	return v0 / v1;
}

void 
get_co(const double T0, const double m, const double *k, const size_t N, 
	   double **A, double **B, double **C, double **D)
{
	*A = malloc_array(N + 1);
	*B = malloc_array(N + 1);
	*C = malloc_array(N + 1);
	*D = malloc_array(N + 1);
	
	double h = 1. / N;
	double 	r3 = 3. * Ra * Ra;
	double 	u00 = Up(temper(0., T0, m)),
			u01 = Up(temper(h / 2., T0, m)),
			uN1 = Up(temper(1. - h / 2., T0, m)),
			uN0 = Up(temper(1., T0, m));
	
	for (int j = 1; j <= N; j++)		/* filling arrays A, C */
	{
		double x = j * h;
		(*A)[j] = (x - h / 2.) / k[2 * j - 1];
		(*C)[j-1] = (*A)[j];
	}
	
	for (int j = 1; j < N; j++)			/* filling arrays B, D */
	{
		double x = j * h;
		(*B)[j] = (*A)[j] + (*C)[j] + r3 * h * x * k[2 * j];
		
		double U = Up(temper(x, T0, m));
		(*D)[j] = r3 * h * x * k[2 * j] * U;
	}
	
	(*A)[0] = 0.;
	(*B)[0] = r3 * h * h * (2. * k[0] + k[1]) / 8. + 1./ k[1];
	(*C)[0] = 1./ k[1] - r3 * h * h * k[1] / 8.;
	(*D)[0] = r3 * h * h * (k[1] * u01 + k[0] * u00) / 4.;
	
	double tmp = r3 * h * h * (2. - h/2.) / 32.;
	(*A)[N] = (*C)[N-1] - tmp * k[2*N - 1];
#define mN 0.5
	(*B)[N] = 3. * Ra * mN + (*C)[N-1] + tmp * (2 * k[2*N] + k[2*N - 1]);
#undef mN
	(*C)[N] = 0.;
	(*D)[N] = 2. * tmp * (k[2 * N - 1] * uN1 + k[2 * N] * uN0);
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

double 
integ(const int index, const double T0, const double m, const int N, const double *k, const double *U)
{
	double h = 1. / N;
	double res = 0.;
	double x = - h / 2;
	for (int j = 0; j < index; j++)
	{
		x += h;
		res += k[2 * j + 1] * (Up(temper(x, T0, m)) - (U[j] + U[j+1]) / 2.) * x * x * 0.5;
	}
	return Ch * res / index;
}

void
solve_tridiagonal(const char *src, const double T0, const double m, const size_t N,
				  double **k, double **f, double **y)
{
	*k = getKi(N, src, T0, m);
	if (*k)
	{
		// display_a(2 * N + 1, k);
		double *A = NULL, *B = NULL, *C = NULL, *D = NULL;
		get_co(T0, m, *k, N, &A, &B, &C, &D);
		
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
		
		*f = malloc_array(N + 1);
		(*f)[0] = 0.;
		
		for (int j = 0; j < N; j++)
		{
			(*f)[j+1] = integ(j+1, T0, m, N, *k, *y);
			// printf("\t%E\t%E\n", (*f)[j+1], (*y)[j+1]);
			// f[j] = - Ch * (y[j+1] - y[j-1]) / (6. * Ra * k[2 * j] * h);	/* OLD VERSION */
		}
		
		free(p); free(q);
		free(A); free(B); free(C); free(D);
	}
}

double 
Ip(double tf, double freq, const double delta_freq)
{
	double v0 = 2. * freq * delta_freq * Hp * pow((freq / Ch), 2.);
	double v1 = exp(Hp * freq / (Kh * tf)) - 1.;

	return v0 / v1;
}

/// x in range [0,1]
double 
interp_k(double x, const size_t N, const double *k)
{
	int j_left = (int)(N * x);
	
	if (j_left == N)
		return k[N];

	double Lx = (double)j_left / N,
		   Rx = (double)(j_left + 1) / N;
	return quick_interpol_linear(x, Lx, Rx, k[j_left], k[j_left + 1]);
}

/// x1, x2 - projection relation, range [0,1]
double
phi(double theta, double x1, double x2, const double dL, const size_t N, const double *k)
{
	if (abs(sin(theta)) < 1E-2)
		return 1.;
	
	double res = (interp_k(x1, 2 * N, k) + 
				  interp_k(x2, 2 * N, k)) / 2.;
	
	double dx = (x2 - x1) / N;
	for (int j = 1; j < N; j++)
	{
		double x = x1 + j * dx;
		res += interp_k(x, 2 * N, k);
	}
	
	return exp(- res * dL);
}

double
intensity(double r, double p, double theta, const double T0, const double m, const size_t N, const double *k)
{
	// NOTICE when theta = 0;
	double x = r * cos(p) + sqrt(pow(Ra, 2.) - pow(r * sin(p), 2.));	// length of projection
	double dx = x / N;
	
	double z0 = r / Ra, zN = 1.;
	double dL = dx / sin(theta);
	
	double dfreq = Df / N;
	
	double res = (interp_k(z0, 2 * N, k) * Ip(temper(z0, T0, m),  left_freq, dfreq) * phi(theta, z0, zN, dL, N, k) + 
				  interp_k(zN, 2 * N, k) * Ip(temper(zN, T0, m), right_freq, dfreq) * phi(theta, zN, zN, dL, N, k)) / 2.;
	
	for (int j = 1; j < N; j++)
	{
		double x = j * dx;
		double zX = sqrt(r * r + x * x - 2. * x * r * cos(p)) / Ra;
		double freq = left_freq + j * dfreq;
		
		res += interp_k(zX, 2 * N, k) * Ip(temper(zX, T0, m), freq, dfreq) * phi(theta, zX, zN, dL, N, k);
	}
	
	return dL * res;
}

double 
u_help(double r, double theta, const double T0, const double m, const size_t N, const double *k)
{
	// if (abs(sin(theta)) < 1E-2)
		// return 0.;
	// printf("\n===================\n");
	double dp = PI / N;
	double res = (intensity(r, 0., theta, T0, m, N, k) + 
				  intensity(r, PI, theta, T0, m, N, k)) / 2.;
	
	for (int j = 1; j < N; j++)
	{
		double p = j * dp;
		res += intensity(r, p, theta, T0, m, N, k);
		// printf("\t%.5f\t%E\n", p, intensity(r, p, theta, T0, m, N, k));
	}
	
	return dp * res;
}

double
Ui(double r, const double T0, const double m, const size_t N, const double *k)
{
	double dtheta = PI / 2. / N;
	double t0 = dtheta / 2.;
	
	double res = 0.;
	for (int j = 0; j < N; j++)
	{
		double theta = t0 + j * dtheta;
		res += u_help(r, theta, T0, m, N, k);
	}
	
	return 4. * dtheta * res;
}

double* 
save_Ui(const double T0, const double m, const size_t N, const double *k)
{
	double *u = malloc_array(N + 1);
	double dr = Ra / N;
	
	for (int j = 0; j <= N; j++)
	{
		double r = j * dr;
		u[j] = Ui(r, T0, m, N, k);
	}
	
	return u;
}

void 
save_solution(const char *src, const double T0, const double m, const size_t N)
{
	//double *k = NULL, *f = NULL, *y = NULL;
	//solve_tridiagonal(src, T0, m, N, &k, &f, &y);
	//double *q = get_q(T0, m, N, k, y);
	
	// display_a(2 * N + 1, k);
	// DO MORE STUFF HERE!!!
	// printf("\n+ Q: %E\n", q(T0, m, N, k, y));
	
	/* Update version */
	//double *U = save_Ui(T0, m, N, k);
	
	// for (int j = 0; j <= N; j++)
	//	printf("\t%E\t%E\n", U[j], y[j]);
	//	printf("\t%.5f\t%E\t%E\n", interp_k((double)j / N, 2 * N, k), Ip((double)j / N, T0, m), phi((double)j / N, 1., 1., N, k));
	// display_a(N + 1, U);
	
	// FILE *fs = fopen("result.xls", "w");
	// for (int j = 0; j <= N; j++)
	// {
		// double r = (double)j / N;
		// fprintf(fs, "%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n", r, f[j], r, y[j], r, q[j], r, U[j]);
	// }
	
	// fclose(fs);

	// free(k); free(f); free(y); free(q); free(U);
	double *q = get_q(T0, m, N);
	// display_a(N+1, q);
	free(q);
}

/* 
   =============================================
				NEW PART - Q(z)
   =============================================
 */

double 
Upf(double tf, double freq, const double C1, const double C2)
{
	double v0 = C1 * pow((freq / 3.), 3.);
	double v1 = exp(C2 * freq / tf) - 1.;
	
	return v0 / v1;
}

double*
get_q(const double T0, const double m, const size_t N)
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
	int n0 = 2 * N + 1;
	double *k = malloc_array(n0);			// real ki 
	double **y = malloc_matrix(193, N+1);	// save Ui
	double *k_tab = malloc_array(16);		// copies from given data
	
	double *A = NULL, *B = NULL, *C = NULL, *D = NULL;
	double *p, *q;
	
	for (int i = 0; i < 193; i++)
	{
		A = B = C = D = NULL;
		/* get copy of value (k) from tab */
		for (int j = 0; j < 16; j++)
			k_tab[j] = KT[j][i];
		
		double h = 0.5 / N;
		// get ki when x in [0, 1]
		for (int j = 0; j < n0; j++)
		{
			double t = temper(h * j, T0, m);
			k[j] = interp_log(t, 16, temperature, k_tab);
		}
		
		/* solve tridiagonal system */
		get_co(T0, m, k, N, &A, &B, &C, &D);
		
		p = malloc_array(N + 1);
		q = malloc_array(N + 1);
		
		/// Step #1
		p[1] = C[0] / B[0];
		q[1] = D[0] / B[0];
		
		for (int j = 1; j < N; j++)
		{
			p[j+1] = C[j] / (B[j] - A[j] * p[j]);
			q[j+1] = (D[j] + A[j] * q[j]) / (B[j] - A[j] * p[j]);
		}
		
		/// Step #2: calculate U[i](z)
		y[i][N] = (A[N] * q[N] + D[N]) / (B[N] - A[N] * p[N]);
		for (int j = N; j > 0; j--)
		{
			y[i][j - 1] = p[j] * y[i][j] + q[j];
		}
		
		free(p); free(q);
		free(A); free(B); free(C); free(D);
	}
	
	free(k);
	
	/* ========== CALCULATE Q ========== */
	double *result = malloc_array(N + 1);
		
	double C0 = 8. * PI * Hp * pow(10., 30),
		   C2 = Hp * pow(10., 15)/ Kh;
	
	double dx = 1. / N;
		   
	for (int i = 0; i <= N; i++)
	{
		double x = i * dx;
		double tf = temper(x, T0, m);
		
		result[i] = 0.;
		for (int j = 0; j < 193; j++)
		{
			/* get copy of value (k) from tab */
			for (int jp = 0; jp < 16; jp++)
				k_tab[jp] = KT[jp][j];
			
			double ki = interp_log(tf, 16, temperature, k_tab);
			double freq = (frequency[j] + frequency[j+1]) / 2.;
			double C1 = C0 * (frequency[j+1] - frequency[j]);
			
			double ui = Upf(tf, freq, C1, C2);
			//if (i == 0)
				//printf("%.9f\t%E\t%E\t%E\n", ki, ui, y[j][i]);
			// double ui = Upf2(tf, freq, )
			if (j == 1)
				printf("%.9f\n", ki);
			result[i] += ki * (ui - y[j][i]);
		}
		result[i] *= Ch;
	}

	free(k_tab); free_m(193, y);
	free(temperature); free(frequency); free_m(16, KT);
	return result;
}

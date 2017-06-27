#include "interp.h"

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
#include <math.h>
#include <stdlib.h>

#define addr(i, j) ((i)*n+j)

char dbg_mode, err_mode;

const double eps=0.000000001;

double sum_arr(double *X, int N, int m)
{
	int i;
	double sum =0 ;
	for(i=0; i<N; i++)
	{
		if(m==1) sum += fabs(X[i]);
		 else sum += X[i];
	}
	return sum;
	/*
	int i, j, n;
	// N - размер суммируемого массива X
	for (n = N; n > 1;)
	{ // n - текущая длина суммируемого массива
	  for (i = 0; i < n; i+=2)
	  {
	     j=i/2;
	     if (i+1 == n) X[j] = X[i];
	     else { if(m) X[j] = fabs(X[i]) + fabs(X[i+1]); else X[j] = X[i] + X[i+1]; }
	  }
	  if (n % 2) n = 1 + n/2;
	  else n = n/2;
	}
	return X[0];
	*/
}

double dot_prod(double *v1, double *v2, double *buf, int n_k, int n)
{
	int i, j;
	for(i=n_k, j=0; i<n; i++, j++) buf[j] = v1[i]*v2[i];
	return sum_arr(buf, j, 0);
}

double norm_v(double *v, double *buf, int n_k, int n)
{
	int i, j;
	double res = 0;
	for(i=n_k, j=0; i<n; i++, j++) buf[j] = v[i]*v[i];
	return sqrt(sum_arr(buf, j, 0));
}

void comp_x_k(double *a_k, double *res, double *buf, int n_k, int n)
{
	int i, f=0;
	double _norm;

	for(i=n_k; i<n; i++)
	{
   		if(i==n_k) res[i] = a_k[i] - norm_v(a_k, buf, n_k, n);
		else res[i] = a_k[i];
	}
	_norm = norm_v(res, buf, n_k, n);

	if(_norm>eps) for(i=n_k; i<n; i++) res[i]/=_norm;
}

void matrix_mult_U_v(double *x, double *y, double *res, int n_k, int n)
{
	int i;
	double r = dot_prod(x, y, res, n_k, n);
	
	for(i=n_k; i<n; i++) res[i] = y[i] - 2*r*x[i];
}

void matrix_mult_U_A(double *x, double *res, double *buf, double *A, int n_k, int n)
{
	int i, j;
		for(j=n_k; j<n; j++)
		{
			for(i=n_k; i<n; i++) res[i] = A[addr(i, j)];

			matrix_mult_U_v(x, res, buf, n_k, n);

			for(i=n_k; i<n; i++) A[addr(i, j)] = buf[i];
		}
}

int reverse_move(int *f, double *A, double *B, double *X, int n)
{
	int i, j, k, idx;
	double sum, s, v;

	for(i=n-1; i>=0; i--)
	{
		sum=0; s=0;
		idx=-1;

		for(j=0; j<n; j++)
		{
			if(fabs(A[addr(i, j)])>=eps)
			{
				sum += fabs(A[addr(i, j)]);
				if(idx==-1) idx=j;
				else {
					if(f[j]==0) { f[j]=1; X[j] = 0; }
					else s += X[j]*A[addr(i, j)];
				}
			}
		}
		if(sum<eps && fabs(B[i])<eps) continue;
		if(sum<eps && fabs(B[i])>=eps) return 1;
		
		v = (B[i] - s)/A[addr(i, idx)];
		if(f[idx]!=0) if(fabs(v-X[idx])>=eps) return 1;
		X[idx] = v; f[idx]=1;
	}
	return 0;
}

size_t lss_memsize_01_14(int n)
{
	return 4*n*sizeof(double);
}

int lss_01_14(int n, double *A, double *B, double *X, double *tmp)
{
	double *x_k, *a_k, *buf, sum;
	int i, j, k, n_k = 0;
	int *f;
	//if(dbg_mode) printf("Debug mode");

	buf = &tmp[0];
	x_k =  &tmp[n];
	a_k = &tmp[2*n];
	f = (int *)&tmp[3*n];
	
	for(i=0; i<n; i++) f[i]=0;

	for(j=0; n_k<n-1; j++, n_k++)
	{
		sum = 0;	
		for(i=n_k; i<n; i++) a_k[i] = A[addr(i, j)];
		comp_x_k(a_k, x_k, buf, n_k, n);
		matrix_mult_U_A(x_k, a_k, buf, A, n_k, n);
		matrix_mult_U_v(x_k, B, a_k, n_k, n);
		
		for(i=n_k; i<n; i++) B[i] = a_k[i];
		//for(i=n_k; i<n; i++) sum += A[addr(j, i)];
		for(i=n_k, k=0; i<n; i++, k++) buf[k] = A[addr(j, i)];
		sum=sum_arr(buf, k, 1);
	}
	
	return reverse_move(f, A, B, X, n	);
}


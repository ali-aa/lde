#include <math.h>
#include <time.h>

#define addr(i, j) ((i)*n+(j))

int vec_memsize_01_14(int n)
{
	return (2*n*n+n)*sizeof(double);
}

double norm_vec(double *X, int n)
{
	double sum = 0;
	int i;
	for(i=0; i<n; i++) sum += X[i]*X[i];

	return sqrt(sum);
}

void inverse(double *A, double *E, int n)
{
 	 int i, j, k;
     double tmp;

	// Инициализация единичной матрицы
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			E[addr(i, j)] = (i==j ? 1 : 0);

	for(i=0; i<n; i++)
	{
		// Нормировка строки (первый элемент станет =1)
		tmp = A[addr(i, i)];
		for(j=n-1; j>=0; j--)
		{
			E[addr(i, j)] /= tmp;
			A[addr(i, j)] /= tmp;
		}
		// Исключаем i-ый элемент с каждой строки кроме i-ой
		for(j=0; j<n; j++)
		if (j!=i)
		{
			tmp = A[addr(j, i)];
			for(k=n-1; k>=0; k--)
			{
				E[addr(j, k)] -= E[addr(i, k)]*tmp;
				A[addr(j, k)] -= A[addr(i, k)]*tmp;
			}
		}
	}
}

void make_A(double *A, double lambda, int n)
{
	int i,j;
	for(i=0; i<n; i++) A[addr(i, i)] -= lambda;
}

int mult_vec(double *A, double *Xk, double *Xk_1, int n)
{
	int i, j, k;

	for(i=0; i<n; i++)
	{
		Xk_1[i] = 0;
		for(k=0; k<n; k++) Xk_1[i] += A[addr(i, k)] * Xk[k];
	}
}

int get_vec(double *A, double *Xk_1, double epsilon, double *Xk, int n, int max_iterations)
{
	int j, f, l;
	double norm, delta = 0;
	
	for(j=0; j<n; j++) { Xk[j]= rand()%3 + delta; delta += 0.1; }
	
	norm = norm_vec(Xk, n);
	for(l=0; l<n; l++) Xk[l] /= norm;

	for(j=1; j<=max_iterations; j++)
	{
		mult_vec(A, Xk, Xk_1, n);
		norm = norm_vec(Xk_1, n);
		// Нормируем полученный вектор
		for(l=0; l<n; l++) Xk_1[l] /= norm;

		// Проверяем на сходимость
		f=1;
		for(l=0; l<n; l++) if(fabs(Xk[l] - Xk_1[l])>epsilon) { f=0; break; }
		if(f==1) return 1;
		for(l=0; l<n; l++) Xk[l] = Xk_1[l];
	}
	return -1;
}

int vec_01_14(int n, int max_iterations, double epsilon, double* A, double* E, double* V, double* tmp, double precision)
{
	// Обратный степенной метод
	int i, j, k;
	double lambda, *C1, *C2, *X;
	
	X = &tmp[0];
	C1 = &tmp[n];
	C2 = &tmp[n*n+n];

	for(k=0; k<n; k++)
	{
		for(i=0; i<n; i++)
			for(j=0; j<n; j++) C1[addr(i, j)] = A[addr(i, j)];
		
		lambda = E[k];
		// Производим сдвиг на собственное значение
		make_A(C1, lambda-epsilon*0.1, n);

		// Находим обратную матрицу
		inverse(C1, C2, n);

		get_vec(C2, &C1[addr(k, 0)], epsilon, X, n, max_iterations);
		for(i=0; i<n; i++) V[addr(i, k)] = C1[addr(k, i)];

	}

}


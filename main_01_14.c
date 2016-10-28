#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "task_01_14.h"

#define EXIT -1
#define ERR_MEM -2
#define ERR_FILE -3
#define ERR_FILE_FORMAT -4
#define OK 0
#define addr(i, j) ((i)*n+(j))
#define MAX(a, b) ((a)>(b) ? (a) : (b) )

extern char dbg_mode = 0, err_mode = 0;

void _(char *msg, int line)
{
	if(dbg_mode) printf("%s:%d: %s\n", "file", line, msg);
}

void err(char *msg)
{
	if(err_mode) fprintf(stderr, "Error: %s\n", msg);
}

void free_res(double *A, double *B, double *E, double *V, double *Y, double *tmp, FILE *fin, FILE *fout)
{
	_("Deallocating resources", __LINE__);
	if(A!=NULL) free(A);
	if(B!=NULL) free(B);
	if(E!=NULL) free(E);
	if(V!=NULL) free(V);
	if(Y!=NULL) free(Y);
	if(tmp!=NULL) free(tmp);
	if(fin!=NULL) free(fin);
	if(fout!=NULL) free(fout);
}

int get_param(char *param, int *max_evc_iter, int *max_vec_iter, double *eps, double *prec)
{
	int i = 0, iter_n = 13, iter_v_n = 13, eps_n = 4, prec_n = 5;
	char iter_s[] = "max_evc_iter=";
	char iter_v_s[] = "max_vec_iter=";
	char eps_s[] = "eps=";
	char prec_s[] = "prec=";

	while(param[i]!=0 && iter_s[i]!=0)
		if(param[i]!=iter_s[i]) break; else i++;
	if(i==iter_n) { *max_evc_iter = atoi(param+iter_n); return OK; }

	while(param[i]!=0 && iter_v_s[i]!=0)
		if(param[i]!=iter_v_s[i]) break; else i++;
	if(i==iter_v_n) { *max_vec_iter = atoi(param+iter_v_n); return OK; }

	i=0;
	while(param[i]!=0 && eps_s[i]!=0)
		if(param[i]!=eps_s[i]) break; else i++;
	if(i==eps_n) { *eps = atof(param+eps_n); return OK; }

	i=0;
	while(param[i]!=0 && prec_s[i]!=0)
		if(param[i]!=prec_s[i]) break; else i++;
	if(i==prec_n) { *prec = atof(param+prec_n); return OK; }

	return EXIT;
}

void print_matrix(int n, double *A)
{
	int i, j;
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++) printf ("%lf ", A[addr(i, j)]);
		printf("\n");
	}
	printf("\n");
}

void print_eig_val(int n, double *E)
{
	int i;
	for(i=0; i<n; i++) printf("%lf ", E[i]);
	printf("\n");
}

int parse_cmdline(int argc, char **argv, int *in_file, int *out_file, char *out_mode, char *stat_mode, int *max_evc_iter, int *max_vec_iter, double *eps, double *prec)
{
	char help[] = "Usage: lde input_file_name output_file_name [options]\n \
Where options include:\n \
  -d    print debug messages [default OFF]\n \
  -e    print errors [default OFF]\n \
  -p    print matrix [default OFF]\n \
  -t    print execution time [default OFF]\n \
  -max_evc_iter=<num>   limit number of iterations for evc module\n \
  			[default - 0, i.e. either not limit or compute limit automatically]\n \
  -max_vec_iter=<num>   limit number of iterations for vec module\n \
			[default - 0, i.e. either not limit or compute limit automatically]\n \
  -h, -?    print this and exit\n \
Default input_file_name value is 01_14_in.txt, default output_file_name value is 01_14_out.txt.\n";

	int i, count=0;

	for(i=1; i<argc; i++)
	{
		if(argv[i][0]=='-' && argv[i][2]==0) {
			char t = argv[i][1];
			if(t=='d') dbg_mode=1;
			else if(t=='e') err_mode=1;
			else if(t=='p') *out_mode=1;
			else if(t=='t') *stat_mode=1;
			else if(t=='h' || t=='?') { printf(help); return OK; }
			else { printf("Error: Wrong usage!\n\n%s", help); return EXIT; }
		}
		else if(argv[i][0]=='-' && argv[i][2]!=0) {
			if(get_param(argv[i]+1, max_evc_iter, max_vec_iter, eps, prec)==-1) {
				printf("Error: Wrong usage!\n\n%s", help); return EXIT;
			}
		}
		else {
			if(count==0) { *in_file=i; count++; }
			else if(count==1) { *out_file=i; count++; }
			else { printf("Error: Wrong usage!\n\n%s", help); return EXIT; }
		}
	}
	return OK;
}

int main(int argc, char **argv)
{
	double *A=NULL, *B=NULL, *E=NULL, *tmp=NULL, *V=NULL, *Y=NULL;
	FILE *fin, *fout;
	int n, i, j, sz_tmp, max_evc_iter=0, max_vec_iter=120, res=0, in_file=-1, out_file=-1;
	double prec=1e-14, eps=1e-6;
	char out_mode = 0, stat_mode = 0;
	clock_t t;	

	if(parse_cmdline(argc, argv, &in_file, &out_file, &out_mode, &stat_mode, &max_evc_iter, &max_vec_iter, &eps, &prec)==EXIT) return EXIT;

	_("Opening files", __LINE__);
	fin = fopen(in_file==-1 ? "01_14_in.txt" : argv[in_file], "r");
	fout = fopen(out_file==-1 ? "01_14_out.txt" : argv[out_file], "w");
	
	if(fin==NULL || fout==NULL) {
		_("Can't open input or/and output file!", __LINE__);
		err("Can't open input or/and output file!");
		free_res(A, B, E, V, Y, tmp, fin, fout);				
		return ERR_FILE;
	}
	
	_("Reading input file", __LINE__);
	if(fscanf(fin, "%d", &n)==-1) {
		err("Wrong file format!");
		free_res(A, B, E, V, Y, tmp, fin, fout);
		return ERR_FILE_FORMAT;
	}

	sz_tmp = vec_memsize_01_14(n);

	_("Allocating memory", __LINE__);
	A=(double *)malloc(sizeof(double)*n*n);
	B=(double *)malloc(sizeof(double)*n*n);
	tmp=(double *)malloc(sz_tmp);
	E=(double *)malloc(sizeof(double)*n);
	V=(double *)malloc(sizeof(double)*n*n);
	Y=(double *)malloc(sizeof(double)*n);

	if(A==NULL || E==NULL || tmp==NULL)
	{
		_("Memory allocation failed", __LINE__);
		err("Can't allocate memory!");
		free_res(A, B, E, V, Y, tmp, fin, fout);
		return ERR_MEM;
	}

	_("Reading rest of input file", __LINE__);
	for(i=0; i<n; i++)
		for(j=0; j<n; j++)
			if(fscanf(fin, "%lf", &A[i*n+j])==-1) 
			{ 
				err("Wrong file format!"); 
				free_res(A, B, E, V, Y, tmp, fin, fout);
				return ERR_FILE_FORMAT;
			}
	for(i=0; i<n; i++)
		if(fscanf(fin, "%lf", &Y[i])==-1)
		{
				err("Wrong file format!"); 
				free_res(A, B, E, V, Y, tmp, fin, fout);
				return ERR_FILE_FORMAT;
		}

	for(i=0; i<n; i++)
		for(j=0; j<n; j++) B[addr(i, j)] = A[addr(i, j)];

	if(out_mode==1) print_matrix(n, A);
	
	_("Simplifying matrix", __LINE__);
	t = clock();
	res = sim_01_14(n, A, tmp, prec);
	t = clock() - t;
	if(stat_mode) printf("Execution time of simplification part %f seconds\n", ((float)t)/CLOCKS_PER_SEC);

	if(res==-1) { fprintf(fout, "0\n"); free_res(A, B, E, V, Y, tmp, fin, fout); return OK; }
	if(out_mode==1) print_matrix(n, A);
	
	_("Computing eigen values", __LINE__);
	t = clock();
	res = evc_01_14(n, max_evc_iter, prec, A, E, tmp, eps);
	t = clock() - t;
	if(stat_mode) printf("Execution time of computing eigen values %f seconds\n", ((float)t)/CLOCKS_PER_SEC);

	
	if(res==1) { fprintf(fout, "0\n"); free_res(A, B, E, V, Y, tmp, fin, fout); return OK; }
	if(out_mode==1) print_eig_val(n, E);
	
	_("Computing eigen vectors", __LINE__);
	t = clock();
	res = vec_01_14(n, max_vec_iter, eps, B, E, V, tmp, prec);
	t = clock() - t;
	if(stat_mode) printf("Execution time of computing eigen vectors %f seconds\n", ((float)t)/CLOCKS_PER_SEC);

	fprintf(fout, "%d\n", n);
	for(i=0; i<n; i++) fprintf(fout, "%1.9lf\n", E[i]);
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) fprintf(fout, "%1.9lf ", V[addr(i, j)]);
		fprintf(fout, "\n");
	}

	free(tmp);
	sz_tmp = lss_memsize_01_14(n);
	tmp=(double *)malloc(sz_tmp);

	_("Solving LSE", __LINE__);
	for(i=0; i<n; i++) { E[i]=0; }
	t = clock();
	lss_01_14(n, V, Y, E, tmp);
	t = clock() - t;
	if(stat_mode) printf("Execution time of solving LSE %f seconds\n", ((float)t)/CLOCKS_PER_SEC);

	_("Writing to output file", __LINE__);
	for(i=0; i<n; i++) fprintf(fout, "%1.9lf\n", E[i]);

	free_res(A, B, E, V, Y, tmp, fin, fout);

	return OK;
}



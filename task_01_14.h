int sim_01_14(int n, double* A, double* tmp, double precision);
int evc_01_14(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision);
int vec_01_14(int n, int max_iterations, double epsilon, double* A, double* E, double* V, double* tmp, double precision);
int lss_01_14(int n, double *A, double *B, double *X, double *tmp);
void print_matrix(int n, double *A);
void print_eig_val(int n, double *E);

#ifndef TEST_H
#define TEST_H

void generate_random_sdp_matrix(double *A, int n);
void dense_to_csr(const double *A, int n, int *rows, int *cols, double *val, int *nnz_out);
int cholesky_decomposition(double *A, int n);
int solve_linear_system(double *A, double *b, int n);
void test_cg(int n, double tol);


#endif

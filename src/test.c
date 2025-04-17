#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/devoir_2.h"
#include "../include/test.h"
#include "time.h"


void generate_random_sdp_matrix(double *A, int n) {
    for (int i = 0; i < n * n; i++)
        A[i] = ((double)rand() / RAND_MAX);

    for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < n; k++)
                sum += A[i * n + k] * A[j * n + k];
            A[i * n + j] = sum;
            A[j * n + i] = sum;
        }

    for (int i = 0; i < n; i++)
        A[i * n + i] += n;
}

void dense_to_csr(const double *A, int n, int *rows, int *cols, double *val, int *nnz_out) {
    int count = 0;
    rows[0] = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(A[i * n + j]) > 1e-12) {
                val[count] = A[i * n + j];
                cols[count] = j;
                count++;
            }
        }
        rows[i + 1] = count;
    }
    *nnz_out = count;
}


int cholesky_decomposition(double *A, int n) {
    int i, j, k;
    for (i = 0; i < n; i++) {
        if (A[i * (i + 1) / 2 + i] < 0) {
            fprintf(stderr, "Error: matrix is not positive definite\n");
            return -1;
        }
        A[i * (i + 1) / 2 + i] = sqrt(A[i * (i + 1) / 2 + i]);
        for (j = i + 1; j < n; j++) {
            A[j * (j + 1) / 2 + i] /= A[i * (i + 1) / 2 + i];
        }
        for (k = i + 1; k < n; k++) {
            for (j = k; j < n; j++) {
                A[j * (j + 1) / 2 + k] -= A[j * (j + 1) / 2 + i] * A[k * (k + 1) / 2 + i];
            }
        }
    }
    return 0;
}

int solve_linear_system(double *A, double *b, int n) {

    if (cholesky_decomposition(A, n) != 0) {
        return -1;
    }
    
    double *x = (double *)malloc(n * sizeof(double));
    if (x == NULL) {
        fprintf(stderr, "Error: cannot allocate memory for x\n");
        return -1;
    }

    int i, k;
    double sum;

    for (i = 0; i < n; i++) {
        sum = b[i];
        for (k = 0; k < i; k++) {
            sum -= A[i * (i + 1) / 2 + k] * x[k];
        }
        x[i] = sum / A[i * (i + 1) / 2 + i];
    }

    for (i = n - 1; i >= 0; i--) {
        sum = x[i];
        for (k = i + 1; k < n; k++) {
            sum -= A[k * (k + 1) / 2 + i] * x[k];
        }
        x[i] = sum / A[i * (i + 1) / 2 + i];
    }

    for (i = 0; i < n; i++) {
        b[i] = x[i];
    }

    for (int i = 1; i <= n; i++) {
        if (isinf(b[i]) || isnan(b[i])) {
            fprintf(stderr, "Error: no solution possible\n");
            return -1;
        }
    }

    free(x);
    return 0;
}

void test_cg(int n, double tol)
{
    srand(time(NULL));

    // Allocation des structures
    double *A_dense = malloc(n * n * sizeof(double));
    double *A_chol  = malloc((n * (n + 1)) / 2 * sizeof(double));
    double *b       = malloc(n * sizeof(double));
    double *b_ref   = malloc(n * sizeof(double));
    double *x_cg    = calloc(n, sizeof(double));

    // 1. Générer A SDP et b
    generate_random_sdp_matrix(A_dense, n);
    for (int i = 0; i < n; i++) {
        b[i] = ((double)rand()) / RAND_MAX;
        b_ref[i] = b[i];
    }

    // 2. Copier la matrice au format compact (triangulaire inférieure stockée en ligne)
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= i; j++)
            A_chol[i * (i + 1) / 2 + j] = A_dense[i * n + j];

    // 3. Résolution directe avec Cholesky (écrase b_ref avec la solution)
    if (solve_linear_system(A_chol, b_ref, n) != 0) {
        printf("Échec de la résolution directe\n");
        goto cleanup;
    }

    // 4. Conversion de la matrice vers CSR
    int *rows = malloc((n + 1) * sizeof(int));
    int *cols = malloc(n * n * sizeof(int));       
    double *val = malloc(n * n * sizeof(double));   
    int nnz = 0;

    dense_to_csr(A_dense, n, rows, cols, val, &nnz);

    int iterations = CG(n, nnz, rows, cols, val, b, x_cg, tol);

    double err = 0;
    for (int i = 0; i < n; i++) {
        double diff = x_cg[i] - b_ref[i];
        err += diff * diff;
    }
    err = sqrt(err);

    printf("✅ Résolu en %d itérations\n", iterations);
    printf("🔍 Norme de l'erreur (CG vs Cholesky) : %.8e\n", err);

cleanup:
    free(A_dense);
    free(A_chol);
    free(b);
    free(b_ref);
    free(x_cg);
    free(rows);
    free(cols);
    free(val);
}

void test_pcg( int n, double tol){
    printf("\nTesting PCG ...\n");
    srand(time(NULL));

    // Allocation des structures
    double *A_dense = malloc(n * n * sizeof(double));
    double *A_chol  = malloc((n * (n + 1)) / 2 * sizeof(double));
    double *b       = malloc(n * sizeof(double));
    double *b_ref   = malloc(n * sizeof(double));
    double *x_cg    = calloc(n, sizeof(double));
    double *x_pcg   = calloc(n, sizeof(double));

    // 1. Générer A SDP et b
    generate_random_sdp_matrix(A_dense, n);
    for (int i = 0; i < n; i++) {
        b[i] = ((double)rand()) / RAND_MAX;
        b_ref[i] = b[i];
    }

    // 2. Copier la matrice au format compact (triangulaire inférieure stockée en ligne)
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= i; j++)
            A_chol[i * (i + 1) / 2 + j] = A_dense[i * n + j];

    // 3. Résolution directe avec Cholesky (écrase b_ref avec la solution)
    if (solve_linear_system(A_chol, b_ref, n) != 0) {
        printf("Échec de la résolution directe\n");
        goto cleanup;
    }

    // 4. Conversion de la matrice vers CSR
    int *rows = malloc((n + 1) * sizeof(int));
    int *cols = malloc(n * n * sizeof(int));       
    double *val = malloc(n * n * sizeof(double));   
    int nnz = 0;

    dense_to_csr(A_dense, n, rows, cols, val, &nnz);

    int iterations_pcg = PCG(n, nnz, rows, cols, val, b, x_pcg, tol);
    int iterations_cg = CG(n, nnz, rows, cols, val, b, x_cg, tol);

    double err = 0;
    double err1 = 0;
    for (int i = 0; i < n; i++) {
        // CG vs Cholesky
        double diff = x_cg[i] - b_ref[i];
        err += diff * diff;
        // PCG vs CG
        double diff1 = x_pcg[i] - x_cg[i];
        err1 += diff1 * diff1;
    }
    err = sqrt(err);
    err1 = sqrt(err1);


    printf("✅ Résolu en %d itérations avec CG\n", iterations_cg);
    printf("🔍 Norme de l'erreur (CG vs Cholesky) : %.8e\n", err);
    printf("✅ Résolu en %d itérations avec PCG\n", iterations_pcg);
    printf("🔍 Norme de l'erreur (PCG vs CG) : %.8e\n", err1);    

cleanup:
    free(A_dense);
    free(A_chol);
    free(b);
    free(b_ref);
    free(x_cg);
    free(rows);
    free(cols);
    free(val);
}



void print_csr(int n, const int *rows_idx, const int *cols, const double *values) {
    for (int i = 0; i < n; ++i) {
        printf("Row %d: ", i);
        for (int idx = rows_idx[i]; idx < rows_idx[i + 1]; ++idx) {
            printf("(%d, %.4f) ", cols[idx], values[idx]);
        }
        printf("\n");
    }
}

void test_ILU() {
    printf("Test 1: Matrice 3x3\n");
    int n1 = 3;
    int nnz1 = 7;
    int rows_idx1[] = {0, 2, 5, 7};
    int cols1[]     = {0, 1, 0, 1, 2, 1, 2};
    double A1[]     = {4, -1, -1, 4, -1, -1, 4};
    double L1[7];
    /*
     4  -1   0
    -1   4  -1
     0  -1   4
    */

    ILU(n1, nnz1, rows_idx1, cols1, A1, L1);
    printf("ILU0 factorization (CSR format):\n");
    print_csr(n1, rows_idx1, cols1, L1);

    printf("\nTest 2: Matrice 4x4\n");
    int n2 = 4;
    int nnz2 = 12;
    int rows_idx2[] = {0, 3, 6, 9, 12};
    int cols2[]     = {0, 1, 3, 0, 1, 2, 1, 2, 3, 0, 2, 3};
    double A2[]     = {10, -1, 2, -1, 11, -1, -1, 10, -1, 2, -1, 10};
    double L2[12];
    /*
     10  -1   0   2
     -1  11  -1   0
      0  -1  10  -1
      2   0  -1  10
    */

    ILU(n2, nnz2, rows_idx2, cols2, A2, L2);
    printf("ILU0 factorization (CSR format):\n");
    print_csr(n2, rows_idx2, cols2, L2);

    printf("\nTest 3: Matrice diagonale 3x3\n");
    int n3 = 3;
    int nnz3 = 3;
    int rows_idx3[] = {0, 1, 2, 3};
    int cols3[]     = {0, 1, 2};
    double A3[]     = {5, 8, 3};
    double L3[3];
    /*
     5   0   0
     0   8   0
     0   0   3
    */

    ILU(n3, nnz3, rows_idx3, cols3, A3, L3);
    printf("ILU0 factorization (CSR format):\n");
    print_csr(n3, rows_idx3, cols3, L3);

    printf("\nTest 4: Matrice creuse 4x4\n");
    int n4 = 4;
    int nnz4 = 6;
    int rows_idx4[] = {0, 2, 3, 5, 6};
    int cols4[]     = {0, 2, 1, 0, 2, 3};
    double A4[]     = {4, 2, 3, 2, 5, 6};
    double L4[6];
    /*
     4   0   2   0
     0   3   0   0
     2   0   5   0
     0   0   0   6
    */

    ILU(n4, nnz4, rows_idx4, cols4, A4, L4);
    printf("ILU0 factorization (CSR format):\n");
    print_csr(n4, rows_idx4, cols4, L4);
}
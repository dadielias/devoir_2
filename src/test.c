#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "devoir_2.h"
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
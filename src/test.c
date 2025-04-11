#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "devoir_2.h"

// Insérer ici ta fonction CG(...)

void denseSolve3x3(const double *A_dense, const double *b, double *x) {
    // Résolution directe de A x = b avec 3x3 inversée à la main (juste pour test)
    double a = A_dense[0], b1 = A_dense[1], c = A_dense[2];
    double d = A_dense[3], e = A_dense[4], f = A_dense[5];
    double g = A_dense[6], h = A_dense[7], i = A_dense[8];

    double det = a*(e*i - f*h) - b1*(d*i - f*g) + c*(d*h - e*g);

    if (fabs(det) < 1e-12) {
        printf("Matrix is singular.\n");
        return;
    }

    // Calcul de l'inverse de A à la main
    double inv[9];
    inv[0] =  (e*i - f*h) / det;
    inv[1] = -(b1*i - c*h) / det;
    inv[2] =  (b1*f - c*e) / det;
    inv[3] = -(d*i - f*g) / det;
    inv[4] =  (a*i - c*g) / det;
    inv[5] = -(a*f - c*d) / det;
    inv[6] =  (d*h - e*g) / det;
    inv[7] = -(a*h - b1*g) / det;
    inv[8] =  (a*e - b1*d) / det;

    for (int j = 0; j < 3; j++) {
        x[j] = 0;
        for (int k = 0; k < 3; k++) {
            x[j] += inv[3*j + k] * b[k];
        }
    }
}

void test_CG_vs_direct() {
    // Matrice A (CSR) : 3x3
    int n = 3;
    int nnz = 7;
    int rows_idx[] = {0, 2, 5, 7};
    int cols[] =     {0, 1, 0, 1, 2, 1, 2};
    double A[] =     {4.0, 1.0, 1.0, 3.0, 1.0, 1.0, 2.0};

    double b[] = {1.0, 2.0, 0.0};
    double x_cg[3];
    double x_ref[3];

    int iters = CG(n, nnz, rows_idx, cols, A, b, x_cg, 1e-8);

    // Version dense pour comparaison
    double A_dense[] = {
        4.0, 1.0, 0.0,
        1.0, 3.0, 1.0,
        0.0, 1.0, 2.0
    };
    denseSolve3x3(A_dense, b, x_ref);

    // Affichage
    printf("Résultat CG en %d itérations:\n", iters);
    for (int i = 0; i < 3; i++) {
        printf("x_cg[%d] = %.8f\n", i, x_cg[i]);
    }

    printf("\nSolution exacte (inversion) :\n");
    for (int i = 0; i < 3; i++) {
        printf("x_ref[%d] = %.8f\n", i, x_ref[i]);
    }

    // Erreur
    double err = 0.0;
    for (int i = 0; i < 3; i++) {
        err += (x_cg[i] - x_ref[i]) * (x_cg[i] - x_ref[i]);
    }
    printf("\nNorme de l'erreur CG : %.8e\n", sqrt(err));
}

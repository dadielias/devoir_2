#include "../include/devoir_2.h"
#include "../include/analyse.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

void create_random_sdp_matrix(double *A, int n, double sparsity) {
    for (int i = 0; i < n * n; i++) A[i] = 0.0;

    for (int i = 0; i < n; i++){
        for (int j = i; j < n; j++) {
            if ((double)rand() / RAND_MAX < sparsity) { 
                double value = ((double)rand() / RAND_MAX);
                A[i * n + j] = value;
                A[j * n + i] = value;
            }
        }
    }

        for (int i = 0; i < n; i++) {
            double diagonal_sum = 0.0;
            for (int j = 0; j < n; j++) {
                diagonal_sum += fabs(A[i * n + j]);
            }
            A[i * n + i] = diagonal_sum + 1.0;
        }
}

void create_random_b(double *b, int n) {
    for (int i = 0; i < n; i++) {
        b[i] = ((double)rand()) / RAND_MAX;
    }
}

void denseTOcsr(const double *A, int n, int *rows, int *cols, double *val, int *nnz_out) {
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


int PCG_CG_temp_compare() {
    printf("Début de la comparaison taille/temps CG et PCG\n");


    FILE *file = fopen("./data/PCG_vs_CG.txt", "w");
    if (!file) {
        perror("Erreur lors de l'ouverture du fichier");
        return -1;
    }


    int sizes[] = {5, 10, 15, 20, 50, 75, 100, 150, 200, 500, 1000, 2000}; // Tailles des matrices
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    double sparsity = 0.1;
    int nbr_iter_CG = 0;
    int nbr_iter_PCG = 0;

    for (int s = 0; s < num_sizes; s++) {
        int n = sizes[s];
        int nnz =  n * n; //(int) (sparsity * n) + 1;  // pr être sur de pas avoir de problème de division par 0

        // Allocation mémoire
        int *rows_idx = malloc((n + 1) * sizeof(int));
        int *cols = malloc(nnz * sizeof(int));
        double *A = malloc(nnz * sizeof(double));
        double *b = malloc(n * sizeof(double));
        double *x = malloc(n * sizeof(double));

        // Génération des données
        create_random_sdp_matrix(A, n, sparsity);
        create_random_b(b, n);
        denseTOcsr(A, n, rows_idx, cols, A, &nnz);

        printf("\nTaille de la matrice : %d x %d\n", n, n);
        printf("Nombre de valeurs non nulles : %d\n", nnz);

        // Mesure du temps pour CG
        clock_t start = clock();
        nbr_iter_CG = CG(n, nnz, rows_idx, cols, A, b, x, 1e-10);
        clock_t end = clock();
        double time_CG = (double)(end - start) / CLOCKS_PER_SEC;

        printf("Temps CG : %.6f secondes\n", time_CG);
        printf("Nombre d'itérations CG : %d\n", nbr_iter_CG);

        // Mesure du temps pour PCG
        start = clock();
        nbr_iter_PCG = PCG(n, nnz, rows_idx, cols, A, b, x, 1e-15);
        end = clock();
        double time_PCG = (double)(end - start) / CLOCKS_PER_SEC;

        printf("Temps PCG : %.6f secondes\n", time_PCG);
        printf("Nombre d'itérations PCG : %d\n\n", nbr_iter_PCG);

        // Écriture des résultats dans le fichier
        fprintf(file, "%d %d %.6f %d %.6f %d\n", n, nnz, time_CG, nbr_iter_CG, time_PCG, nbr_iter_PCG);

        // Libération mémoire
        free(rows_idx);
        free(cols);
        free(A);
        free(b);
        free(x);
    }

    fclose(file);
    printf("Résultats enregistrés dans benchmark_results.txt\n");
    return 0;    
}


int main_analyse(){
    if (PCG_CG_temp_compare() != 0) {
        fprintf(stderr, "Erreur lors de la comparaison CG et PCG\n");
        return -1;
    }
    printf("Comparaison CG et PCG terminée avec succès.\n");
    return 0;
}
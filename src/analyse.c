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

void create_tridiag_sdp_matrix(double *A, int n) {
    for (int i = 0; i < n; i++) {
        A[i * n + i] = 2.0;
        if (i > 0) {
            A[i * n + (i - 1)] = -1.0;
            A[(i - 1) * n + i] = -1.0;
        }
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
    printf("Début de la comparaison taille/temps CG, PCG et ILU\n");

    FILE *file = fopen("./data/PCG_vs_CG.txt", "w");
    if (!file) {
        perror("Erreur lors de l'ouverture du fichier");
        return -1;
    }

    // Ajouter un en-tête au fichier pour plus de clarté
    //fprintf(file, "n nnz avg_t_CG avg_iter_CG avg_t_PCG avg_iter_PCG avg_t_ILU\n");

    int sizes[] = {5, 10, 15, 20, 50, 75, 100, 150, 200, 500, 1000, 2000}; // Tailles des matrices
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    double sparsity = 0.1;

    for (int s = 0; s < num_sizes; s++) {
        int n = sizes[s];
        int nnz = n * n; // Nombre de valeurs non nulles (approximatif)

        // Allocation mémoire
        int *rows_idx = malloc((n + 1) * sizeof(int));
        int *cols = malloc(nnz * sizeof(int));
        double *A = malloc(nnz * sizeof(double));
        double *b = malloc(n * sizeof(double));
        double *x = malloc(n * sizeof(double));
        double *L = malloc(nnz * sizeof(double)); // Pour la factorisation ILU

        // Génération des données
        create_random_sdp_matrix(A, n, sparsity);
        create_random_b(b, n);
        denseTOcsr(A, n, rows_idx, cols, A, &nnz);

        printf("\nTaille de la matrice : %d x %d\n", n, n);
        printf("Nombre de valeurs non nulles : %d\n", nnz);

        // Variables pour les moyennes
        double total_time_CG = 0.0, total_time_PCG = 0.0, total_time_ILU = 0.0;
        int total_iter_CG = 0, total_iter_PCG = 0;

        // Répéter 10 fois pour calculer les moyennes
        for (int repeat = 0; repeat < 10; repeat++) {
            // Mesure du temps pour CG
            clock_t start = clock();
            int iter_CG = CG(n, nnz, rows_idx, cols, A, b, x, 1e-15);
            clock_t end = clock();
            total_time_CG += (double)(end - start) / CLOCKS_PER_SEC;
            total_iter_CG += iter_CG;

            // Mesure du temps pour ILU
            start = clock();
            ILU(n, nnz, rows_idx, cols, A, L);
            end = clock();
            total_time_ILU += (double)(end - start) / CLOCKS_PER_SEC;

            // Mesure du temps pour PCG
            start = clock();
            int iter_PCG = PCG(n, nnz, rows_idx, cols, A, b, x, 1e-15);
            end = clock();
            total_time_PCG += (double)(end - start) / CLOCKS_PER_SEC;
            total_iter_PCG += iter_PCG;
        }

        // Calcul des moyennes
        double avg_time_CG = total_time_CG / 10.0;
        double avg_time_PCG = total_time_PCG / 10.0;
        double avg_time_ILU = total_time_ILU / 10.0;
        int avg_iter_CG = total_iter_CG / 10;
        int avg_iter_PCG = total_iter_PCG / 10;

        printf("Temps moyen CG : %.6f secondes, Itérations moyennes CG : %d\n", avg_time_CG, avg_iter_CG);
        printf("Temps moyen ILU : %.6f secondes\n", avg_time_ILU);
        printf("Temps moyen PCG : %.6f secondes, Itérations moyennes PCG : %d\n\n", avg_time_PCG, avg_iter_PCG);

        // Écriture des résultats dans le fichier
        fprintf(file, "%d %d %.6f %d %.6f %d %.6f\n", n, nnz, avg_time_CG, avg_iter_CG, avg_time_PCG, avg_iter_PCG, avg_time_ILU);

        // Libération mémoire
        free(rows_idx);
        free(cols);
        free(A);
        free(b);
        free(x);
        free(L);
    }

    fclose(file);
    printf("Résultats enregistrés dans PCG_vs_CG.txt\n");
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
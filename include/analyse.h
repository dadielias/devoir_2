#ifndef ANALYSE_H
#define ANALYSE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// Déclaration des fonctions
void create_random_sdp_matrix(double *A, int n, double sparsity);
void create_random_b(double *b, int n);
void denseTOcsr(const double *A, int n, int *rows, int *cols, double *val, int *nnz_out);
int PCG_CG_temp_compare();
int main_analyse();

#endif
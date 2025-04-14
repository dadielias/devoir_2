#include "../include/devoir_2.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


/**
 * @brief effectue le produit matrice-vecteur d'une matrice en format CSR
 * @param n taille de la matrice
 * @param rows_idx tableau d'index de lignes
 * @param cols tableau d'index de colonnes
 * @param A tableau de valeurs non nulles
 * @param x vecteur d'entrée
 * @param y vecteur de sortie
 */
static inline void spmv(int n, const int *rows_idx, const int *cols, const double *A, const double *x, double *y)
{
    for (int i = 0; i < n; i++)
    {
        double s = 0;
        for (int j = rows_idx[i]; j < rows_idx[i + 1]; j++)
        {
            s += A[j] * x[cols[j]];
        }
        y[i] = s;
    }
}

/**
 * @brief calcule le résidu de l'équation Ax = b
 * @param n taille de la matrice
 * @param rows_idx tableau d'index de lignes
 * @param cols tableau d'index de colonnes
 * @param A tableau de valeurs non nulles
 * @param x vecteur solution
 * @param b vecteur de droite
 */
static inline void residual(int n, const int *rows_idx, const int *cols, const double *A, const double *x, const double *b, double *r)
{
    spmv(n, rows_idx, cols, A, x, r);
    for (int i = 0; i < n; i++)
    {
        r[i] = b[i] - r[i];
    }
}

/**
 * @brief effectue le produit scalaire de deux vecteurs de taille n
 * @param n taille des vecteurs
 * @param x premier vecteur
 * @param y second vecteur
 */
static inline double dot(int n, const double *x, const double *y)
{
    double s = 0.0;
    for (int i = 0; i < n; i++)
    {
        s += x[i] * y[i];
    }
    return s;
}

/**
 * @brief effectue l'addition de deux vecteurs avec multiplication par un scalaire
 * @param n taille des vecteurs
 * @param x premier vecteur
 * @param y second vecteur
 * @param a scalaire
 */
static inline void axpy(int n, double *x, const double *y, double a)
{
    for (int i = 0; i < n; i++)
    {
        x[i] += a * y[i];
    }
}

void Matvec(
    int n,
    int nnz,
    const int *rows_idx,
    const int *cols,
    const double *A,
    const double *v,
    double *Av
) {}

/**
 * Résout le système linéaire Lx = b où L est une matrice résultant d'une factorisation incomplète.
 *
 * @param L Tableau de sortie de même taille que A, contenant la factorisation incomplète :
 *           - Si csr_sym() renvoie 0 : L contient L dans sa partie strictement inférieure
 *             et U dans sa partie supérieure.
 *           - Si csr_sym() renvoie 1 : L contient L dans sa partie strictement inférieure
 *             et D dans sa diagonale.
 * @param b Vecteur second membre du système linéaire.
 * @param x Vecteur solution du système linéaire, calculé en sortie.
 */
void solve(
    int n,
    int nnz,
    const int *rows_idx,
    const int *cols,
    const double *L,
    const double *b,
    double *x
) {
    // Résoudre Lz = b
    for (int i = 0; i < n; i++) {
        double sum = b[i];
        for (int j = rows_idx[i]; j < rows_idx[i + 1]; j++) {
            if (cols[j] < i) {
                sum -= L[j] * x[cols[j]];
            }
        }
        x[i] = sum;
    }

    // Résoudre Ux = z
    for (int i = n - 1; i >= 0; i--) {
        double sum = x[i];
        for (int j = rows_idx[i]; j < rows_idx[i + 1]; j++) {
            if (cols[j] > i) {
                sum -= L[j] * x[cols[j]];
            }
        }
        x[i] = sum / L[rows_idx[i + 1] - 1];
    }
}

int CG(
    int n,
    int nnz,
    const int *rows_idx,
    const int *cols,
    const double *A,
    const double *b,
    double *x,
    double eps
) {
    double *r = malloc(n * sizeof(double));
    double *p = malloc(n * sizeof(double));
    double *Ap = malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        x[i] = 0.0;
        r[i] = b[i];  // car x = 0 initialement
        p[i] = r[i];
    }

    int max_iter = 10000;
    int iter;
    double r_norm2 = 0.0, r_norm2_new = 0.0, r0_norm2 = 0.0;

    // Norme initiale : r0^T r0
    r0_norm2 = dot(n, r, r);
    r_norm2 = r0_norm2;

    for (iter = 0; iter < max_iter; iter++) {
        
        // Ap = A * p
        spmv(n, rows_idx, cols, A, p, Ap);

        // Calcul alpha : alpha = r^T r / (p^T Ap)
        double alpha_num = r_norm2;  // PK PAS d^T r --> dot(n, r, d);
        double alpha_den = dot(n, p, Ap);

        double alpha = alpha_num / alpha_den;

        // x = x + alpha * p
        axpy(n, x, p, alpha);


        // r = r - alpha * Ap
        axpy(n, r, Ap, -alpha);
        //residual(n, rows_idx, cols, A, x, b, r); -> la mm a priori

        // r^T r
        r_norm2_new = dot(n, r, r);

        // Critère d'arrêt relatif
        if (sqrt(r_norm2_new / r0_norm2) < eps) {
            break;
        }

        double beta = r_norm2_new / r_norm2;

        for (int i = 0; i < n; i++) {
            p[i] = r[i] + beta * p[i];
        }

        r_norm2 = r_norm2_new;
    }

    free(r);
    free(p);
    free(Ap);

    return iter + 1;
}

double get_element_csr(int i, int j, const int *rows_idx, const int *cols, const double *A) {
    // Parcourt les colonnes de la ligne i
    for (int idx = rows_idx[i]; idx < rows_idx[i + 1]; idx++) {
        if (cols[idx] == j) {
            return A[idx]; // Retourne la valeur si la colonne correspond
        }
    }
    return 0.0; // Retourne 0 si l'élément est nul
}


int get_index_csr(int i, int j, const int *rows_idx, const int *cols) {
    for (int idx = rows_idx[i]; idx < rows_idx[i + 1]; idx++) {
        if (cols[idx] == j)
            return idx;
    }
    return -1; 
}

/**
 * @brief Calcule la factorisation incomplète LU = LDL* d'une matrice A en format CSR.
 *
 * @param n Taille de la matrice (nombre de lignes/colonnes).
 * @param nnz Nombre de valeurs non nulles dans la matrice A.
 * @param rows_idx Tableau d'index des lignes (format CSR).
 * @param cols Tableau d'index des colonnes (format CSR).
 * @param A Tableau des valeurs non nulles de la matrice A (format CSR).
 * @param L Tableau de sortie de même taille que A, contenant la factorisation incomplète :
 *           - Si csr_sym() renvoie 0 : L contient L dans sa partie strictement inférieure
 *             et U dans sa partie supérieure.
 *           - Si csr_sym() renvoie 1 : L contient L dans sa partie strictement inférieure
 *             et D dans sa diagonale.
 */
void ILU(
    int n,
    int nnz,
    const int *rows_idx,
    const int *cols,
    const double *A,
    double *L
) {    
    // Copier A dans L en format symétrique
    for (int i = 0; i < nnz; i++) {
        L[i] = A[i];
    }

    double Akk, Lik, Lkj;
    int colj, idx;
    for (int k = 0; k < n; k++){
    
        // Cherche Akk = L[k][k]
        idx = get_index_csr(k, k, rows_idx, cols);
        Akk = L[idx];

        // Si Akk == 0, on ne peut pas continuer
        if (fabs(Akk) < 1e-14 || idx == -1) continue;
        
        // On divise la ligne k par Akk
        for (int j = k + 1; j < n; j++) {
            // L[k][j] = L[k][j] / Akk;
            idx = get_index_csr(j, k, rows_idx, cols);
            if (idx == -1) continue;
            L[idx] /= Akk;

            idx = get_index_csr(k, j, rows_idx, cols);
            if (idx == -1) continue;
            Lkj = L[idx];

            for (int i = k + 1; i < n; i++) {
                idx = get_index_csr(i, k, rows_idx, cols);
                if (idx == -1) continue;
                Lik = L[idx];

                // On met à jour L[i][j] = L[i][j] - L[i][k] * L[k][j]
                idx = get_index_csr(i, j, rows_idx, cols);
                if (idx == -1) continue;
                L[idx] -= Lik * Lkj;
            }
        }
    }
}




int PCG(
    int n,
    int nnz,
    const int *rows_idx,
    const int *cols,
    const double *A,
    const double *b,
    double *x,
    double eps
) {
    double *r = malloc(n * sizeof(double));
    double *p = malloc(n * sizeof(double));
    double *z = malloc(n * sizeof(double));
    double *Ap = malloc(n * sizeof(double));
    double *M = malloc(nnz * sizeof(double));

    // Préconditionneur
    // M = LU  = A = L D L^T
    ILU(n, nnz, rows_idx, cols, A, M);
    
    // M z = r
    solve(n, nnz, rows_idx, cols, A, b, z);

    for (int i = 0; i < n; i++) {
        x[i] = 0.0;
        r[i] = b[i];  // car x = 0 initialement
        p[i] = z[i];
    }

    int max_iter = 10000;
    int iter;
    double r_norm2 = 0.0, r_norm2_new = 0.0, r0_norm2 = 0.0, zpr_old = 0.0;

    // Norme initiale : r0^T r0
    r0_norm2 = dot(n, r, r);
    r_norm2 = r0_norm2;

    for (iter = 0; iter < max_iter; iter++) {
        
        // Ap = A * p
        spmv(n, rows_idx, cols, A, p, Ap);

        // Calcul alpha : alpha = r^T z / (p^T Ap)
        double alpha_num = dot(n, r, z);
        double alpha_den = dot(n, p, Ap);

        double alpha = alpha_num / alpha_den;

        // x = x + alpha * p
        axpy(n, x, p, alpha);


        // r = r - alpha * Ap
        axpy(n, r, Ap, -alpha);
        //residual(n, rows_idx, cols, A, x, b, r); -> la mm a priori

        // r^T r
        r_norm2_new = dot(n, r, r);

        // Critère d'arrêt relatif
        if (sqrt(r_norm2_new / r0_norm2) < eps) {
            break;
        }

        zpr_old = dot(n, z, z);

        // M z = r
        solve(n, nnz, rows_idx, cols, A, r, z);

        double beta = dot(n, r, z) / zpr_old;

        for (int i = 0; i < n; i++) {
            p[i] = z[i] + beta * p[i];
        }

        r_norm2 = r_norm2_new;
    }

    free(r);
    free(p);
    free(Ap);
    free(z);
    return iter + 1;
}

int csr_sym() {
    return 0; //  Both parts
    // return 1; // Lower part only
}

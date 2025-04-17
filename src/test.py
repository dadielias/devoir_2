import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spilu

def print_csr(n, rows_idx, cols, values):
    """Affiche une matrice au format CSR."""
    for i in range(n):
        print(f"Row {i}: ", end="")
        for idx in range(rows_idx[i], rows_idx[i + 1]):
            print(f"({cols[idx]}, {values[idx]:.4f}) ", end="")
        print()

def test_ilu():
    print("Test 1: Matrice 3x3")
    A1_dense = np.array([
        [4, -1,  0],
        [-1, 4, -1],
        [0, -1, 4]
    ], dtype=np.float64)
    A1_csr = csr_matrix(A1_dense)
    ilu1 = spilu(A1_csr, fill_factor=1.0, drop_tol=0.0)
    print("L =\n", ilu1.L.toarray())
    print("U =\n", ilu1.U.toarray())

    print("\nTest 2: Matrice 4x4")
    A2_dense = np.array([
        [10, -1,  0,  2],
        [-1, 11, -1,  0],
        [0, -1, 10, -1],
        [2,  0, -1, 10]
    ], dtype=np.float64)
    A2_csr = csr_matrix(A2_dense)
    ilu2 = spilu(A2_csr, fill_factor=1.0, drop_tol=0.0)
    print("L =\n", ilu2.L.toarray())
    print("U =\n", ilu2.U.toarray())

    print("\nTest 3: Matrice diagonale 3x3")
    A3_dense = np.array([
        [5, 0, 0],
        [0, 8, 0],
        [0, 0, 3]
    ], dtype=np.float64)
    A3_csr = csr_matrix(A3_dense)
    ilu3 = spilu(A3_csr, fill_factor=1.0, drop_tol=0.0)
    print("L =\n", ilu3.L.toarray())
    print("U =\n", ilu3.U.toarray())

    print("\nTest 4: Matrice creuse 4x4")
    A4_dense = np.array([
        [4, 0, 2, 0],
        [0, 3, 0, 0],
        [2, 0, 5, 0],
        [0, 0, 0, 6]
    ], dtype=np.float64)
    A4_csr = csr_matrix(A4_dense)
    ilu4 = spilu(A4_csr, fill_factor=1.0, drop_tol=0.0)
    print("L =\n", ilu4.L.toarray())
    print("U =\n", ilu4.U.toarray())

test_ilu()
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spilu

# Matrice 3x3 comme dans ton exemple
A_dense = np.array([
    [4, -1,  0],
    [-1, 4, -1],
    [0, -1, 4]
], dtype=np.float64)

# Conversion CSR
A_csr = csr_matrix(A_dense)

# Factorisation ILU0 (remplissage nul)
ilu = spilu(A_csr, fill_factor=1.0, drop_tol=0.0)

# Récupérer les matrices L et U (format CSC)
L = ilu.L 
U = ilu.U

print("L =\n", ilu.L.toarray())
print("U =\n", ilu.U.toarray())


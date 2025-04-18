import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Fonction pour charger un fichier CSV (2 colonnes : itérations, erreurs)
def load_csv(filename):
    data = np.loadtxt(filename, delimiter=',')
    k = data[:, 0]  # Itérations
    err = data[:, 1]  # Erreurs relatives
    return k, err


sns.set_theme(style="darkgrid")

# Charger les données
k_cg, err_cg = load_csv("cg_residuals.csv")
k_pcg, err_pcg = load_csv("pcg_residuals.csv")

# Créer le graphique
plt.figure(figsize=(8, 5))
plt.semilogy(k_cg, err_cg, label='CG', linestyle='-', color='blue', linewidth=1)
plt.semilogy(k_pcg, err_pcg, label='PCG', linestyle='-', color='orange', linewidth=1)

plt.xlabel("Itérations k", fontname="Arial")
plt.ylabel("Erreur relative ‖r_k‖ / ‖r_0‖", fontname="Arial")
plt.title("Convergence de CG vs PCG", fontname="Arial")
plt.grid(True, which='both', linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig("convergence_graph.png")  # Enregistre l’image
plt.show()

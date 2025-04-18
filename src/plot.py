import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Charger les données depuis le fichier
data_file = "./data/PCG_vs_CG.txt"
data = np.loadtxt(data_file)


# Extraire les colonnes
matrix_sizes = data[:, 0]  # Taille des matrices
iterations_CG = data[:, 3]  # Nombre d'itérations pour CG
iterations_PCG = data[:, 5]  # Nombre d'itérations pour PCG
time_CG = data[:, 2]  # Temps pour CG
time_PCG = data[:, 4]  # Temps pour PCG
time_ILU = data[:, 6]  # Temps pour ILU
time_PCG_minus_ILU = time_PCG - time_ILU  # Temps PCG sans le coût d'ILU

# Plot 1 : Comparaison du nombre d'itérations entre CG et PCG
plt.figure(figsize=(10, 6))
plt.plot(matrix_sizes, iterations_CG, label="CG", marker='o', linestyle='-', color='blue')
plt.plot(matrix_sizes, iterations_PCG, label="PCG", marker='s', linestyle='--', color='red')

# Ajouter des labels et une légende
plt.xlabel("Taille de la matrice (n)", fontsize=12)
plt.ylabel("Nombre d'itérations", fontsize=12)
plt.title("Comparaison du nombre d'itérations : CG vs PCG", fontsize=14)
plt.legend(fontsize=12)
plt.grid(which="both", linestyle="--", linewidth=0.5)

# Sauvegarder et afficher le graphique
plt.savefig("./data/iterations_CG_vs_PCG.png", dpi=300)

# Plot 2 : Évolution temporelle de ILU, CG, PCG et PCG sans ILU
plt.figure(figsize=(10, 6))
plt.loglog(matrix_sizes, time_CG, label="CG")
plt.loglog(matrix_sizes, time_PCG, label="PCG")
plt.loglog(matrix_sizes, time_ILU, label="ILU", linestyle='dotted', color='green')

# Ajouter des labels et une légende
plt.xlabel("Taille de la matrice log(n)", fontsize=12)
plt.ylabel("Temps d'exécution log(secondes)", fontsize=12)
plt.title("Évolution temporelle : CG, PCG, ILU)", fontsize=14)
plt.legend(fontsize=12)
plt.grid(which="both", linestyle="--", linewidth=0.5)

# Sauvegarder et afficher le graphique
plt.savefig("./data/times_CG_PCG_ILU.png", dpi=300)

# Plot 3 : Comparaison CG vs PCG sans ILU
plt.figure(figsize=(10, 6))
plt.loglog(matrix_sizes, time_CG, label="CG")
plt.loglog(matrix_sizes, time_PCG_minus_ILU, label="PCG (sans ILU)")

# Ajouter des labels et une légende
plt.xlabel("Taille de la matrice log(n)", fontsize=12)
plt.ylabel("Temps d'exécution log(secondes)", fontsize=12)
plt.title("Comparaison temporelle : CG vs PCG (sans ILU)", fontsize=14)
plt.legend(fontsize=12)
plt.grid(which="both", linestyle="--", linewidth=0.5)

# Sauvegarder et afficher le graphique
plt.savefig("./data/times_CG_vs_PCG_without_ILU.png", dpi=300)

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


import numpy as np
import matplotlib.pyplot as plt

# Charger les données depuis le fichier
data_file = "./data/PCG_vs_CG.txt"
data = np.loadtxt(data_file)

# Extraire les colonnes
matrix_sizes = data[:, 0]  # Taille des matrices
time_CG = data[:, 2]       # Temps pour CG
time_PCG = data[:, 4]      # Temps pour PCG

# Tracer le graphique log-log
plt.figure(figsize=(10, 6))
plt.loglog(matrix_sizes, time_CG, label="CG", marker='o', linestyle='-', color='blue')
plt.loglog(matrix_sizes, time_PCG, label="PCG", marker='s', linestyle='--', color='red')

# Ajouter des labels et une légende
plt.xlabel("Taille de la matrice (n)", fontsize=12)
plt.ylabel("Temps d'exécution (secondes)", fontsize=12)
plt.title("Comparaison de la complexité temporelle : CG vs PCG", fontsize=14)
plt.legend(fontsize=12)
plt.grid(which="both", linestyle="--", linewidth=0.5)

# Sauvegarder et afficher le graphique
plt.savefig("./data/CG_vs_PCG_loglog.png", dpi=300)
plt.show()
# Import functions from the pipeline module
from plot.plotting import plot_pc_by_species

#########################################################################
# To execute correctly this script, use this line for running:
# PYTHONPATH=scripts python scripts/additional/analyze_general.py
#########################################################################

n_genomes = 10000 # Number of genomes per species


outdir = "results/Ecoli_complete"

# csv with info
csv_path = f"{outdir}/final.csv"
# Plots the distribution of promoters among species and Pc variants
count_table = plot_pc_by_species(csv_path = csv_path,
                                 title = f"Distribution of Pc variants by Species - db: E coli, total genomes: {n_genomes}",
                                 outdir = outdir)

print(count_table)
# Saves results to csv
count_table.to_csv(f"{outdir}/count_table.csv")

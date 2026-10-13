# Import functions from the pipeline module
from plot.plotting import plot_pc_by_species

n_genomes = 1500 # Number of genomes per species
# csv with info
path = f"../results/final/{n_genomes}genomes/final.csv"

outdir = f"../results/final/{n_genomes}genomes"
# Plots the distribution of promoters among species and Pc variants
count_table = plot_pc_by_species(csv_path = path,
                                 title = f"Distribution of Pc variants by Species - db: {n_genomes * 5} genomes ESKAPEE",
                                 outdir = outdir)

print(count_table)
# Saves results to csv
count_table.to_csv(f"{outdir}/count_table.csv", index=False)

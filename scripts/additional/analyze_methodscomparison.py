# Import libraries 
import os
import glob
import pandas as pd

# Import functions from the pipeline module
from pipeline.blastn import process_blast_tsv
from pipeline.useful_funct import contig_to_genome, genome_key
from pipeline.analyze_Pc import extract_integron_info_gbk

# -------- Paths with results to analyze
n = 500 # Number of genomes per species
genomes_dir = f"../data_practicas/{n}genomes"
tsv = f"../results_practicas/{n}genomes/blast/blast_Int1_Pc_pattern.fa.tsv"
if2_dir = f"../results_practicas/{n}genomes/IF2"
csv = f"results/plots/{n}genomes/genomes_hits.csv"

out_dir = f"results/plots/{n}genomes/"

# -------- 1) genomes with different number of hits in BLAST and IF2
df = pd.read_csv(csv, 
                 dtype = {"genome": str}) # to keep the acc as string
disc = df[df["BLAST_hits"] != df["IF2_hits"]] # discordant df
df_disc_tocsv = disc[["genome", "BLAST_hits", "IF2_hits"]]

disc_genomes = disc["genome"].tolist() # Transform the column of discordant genomes into a list

print("Discordant genomes:")
print(df_disc_tocsv.to_string(index = False))

# Save into a csv
df_disc_tocsv.to_csv(f"{out_dir}/discordant.csv", index = False)

# -------- 2) BLAST rows of those genomes
blast_df = process_blast_tsv(tsv)  # Adds strand, start and end
mapping = contig_to_genome(genomes_dir)
contig2acc = {} # Dictionary to store mapping contig - acc
for contig, path in mapping.items(): # Iterates over contig ID and its path
    contig2acc[contig] = genome_key(path)

blast_df["genome"] = blast_df["seq_id"].map(contig2acc)

### CONTROL:
# Contigs not found in the mapping have no genome (NaN)
n_total = len(blast_df) # Total BLAST rows
n_no_genome = blast_df["genome"].isna().sum() # Rows whose contig was not in the mapping

print(f"\nBLAST rows without genome: {n_no_genome} of {n_total}")

# Cols of the output of blast
cols = ["genome", "seq_id", "identity", "length", "query_start", 
        "query_end", "start", "end", "strand", "Evalue"]

# Keep only the rows whose genome is in the discordant list
sub = blast_df[blast_df["genome"].isin(disc_genomes)] 
sub = sub.sort_values(["genome", "seq_id", "start"])

print("\nBLAST rows of the discordant genomes:")
print(sub.to_string(index = False))

sub.to_csv(f"{out_dir}/discordant_blast.csv", index = False)

# -------- 3) IF2 integrases of those genomes (0-based start, as in the gbk)
rows = []
for genome in disc_genomes:
    # Collect the .gbk files of the IF2 result folders of this genome
    gbk_files = []
    for folder in glob.glob(os.path.join(if2_dir, f"Results_Integron_Finder_*{genome}*")):
        gbk_files += glob.glob(os.path.join(folder, "*.gbk"))
    info = extract_integron_info_gbk(genome, gbk_files)
    for i in info:
        if not i["calin"]: # Skip CALIN (no integrase)
            rows.append({"genome": genome, 
                         "locus_id": i["locus_id"],
                         "start": i["start"], 
                         "end": i["end"],
                         "strand": i["strand"], 
                         "integron_type": i["integron_type"],
                         "integron_id": i["integron_id"]})
if2_df = pd.DataFrame(rows)
print("\nIF2 integrases of the discordant genomes:")
print(if2_df.to_string(index = False))
if2_df.to_csv(f"{out_dir}/discordant_IF2.csv", index = False)

import os
import time
import subprocess
import glob
import shutil
import pandas as pd
import numpy as np
import shlex
from Bio import SeqIO

from pipeline.useful_funct import contig_to_genome
from pipeline.metadata_info import add_genome_metadata

# --------Definition of functions
#------------------------ Function to run blastn
def run_blastn(query, output, db, max_target, evalue_threshold = 10):
    """
    Executes the blastn command with the specified query, output file, and database.

    Args
    ------
    query (str): The path to the query fasta file.
    output (str): The path to the output file where results will be saved.
    db (str): The path to the blast database.
    max_target (int): Maximum number of aligned sequences keeped into the results.
    evalue_threshold (int): E-value threshold.
                            Default to 10.

    Return
    ------
    output (str): The path to the output file where results are saved.
    total_time (float): The time taken to execute the blastn command in seconds.
    """
    # Creates the output directory to avoid errors
    os.makedirs(os.path.dirname(output), exist_ok=True)
    # Blastn command
    blast_cmd = ["blastn", "-query", query, # Query fasta file
                 "-db", db, # Blast database
                 "-out", output, 
                 "-max_target_seqs", str(max_target), # Max number of aligned seqs to keep
                 "-evalue", str(evalue_threshold),
                 "-outfmt", "6"] # Output format 6 (tabular)
    start = time.time() # Start time to measure execution time
    # Execute the command and handle errors
    print(f"Executing: \n{shlex.join(blast_cmd)}")
    try:
        subprocess.run(blast_cmd, 
                       check = True, capture_output = True, text = True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing blastn:\n{e.stderr}")
        return

    total_time = time.time() - start # End time to measure execution time

    return output, total_time

#------------------------ Function to process a BLAST tabular output (outfmt 6)
def process_blast_tsv(path, min_length = None):
    """
    Reads a BLAST tabular output (outfmt 6) and adds strand and normalized
    subject coordinates.

    Args
    ------
    path (str): The path to the BLAST tabular output file.
    min_length (int): Minimum length of hits to be considered.
                      Default to None.

    Return
    ------
    df (pd.DataFrame): One row per hit, with the 12 standard outfmt 6 columns plus:
                       - strand (str): "+" if seq_start <= seq_end, "-" otherwise.
                       - start (int): Lowest subject coordinate of the hit.
                       - end (int): Highest subject coordinate of the hit.
                       
                       If the file is empty, an empty DataFrame with the same
                       columns is returned.
    """
    # Column names of the BLAST tabular output (outfmt 6)
    columns = ["query_id", "seq_id", "identity", "length", "mismatch", "gapopen", 
               "query_start", "query_end", "seq_start", "seq_end", "Evalue", "BitScore"]
    try: # Reads the csv file
        df = pd.read_csv(path, sep = "\t", names = columns)
    except pd.errors.EmptyDataError:
        # If no hits, keeps the expected columns 
        df = pd.DataFrame(columns = columns)
    if min_length:
        df = df[df["length"] >= min_length] # Filters hits by minimum length
    # If start is lower than end, the strand is positive
    pos_strand = df["seq_start"] <= df["seq_end"]
    # Assign the strand column based on the comparison of seq_start and seq_end
    df["strand"] = np.where(pos_strand, "+", "-")  
    # Assign the start and end columns based on the minimum and maximum of seq_start and seq_end
    df["start"] = df[["seq_start", "seq_end"]].min(axis=1)
    df["end"] = df[["seq_start", "seq_end"]].max(axis=1)
    
    return df

#------------------------ Function to extract the nucleotide sequence of each BLAST hit
def extract_hit_sequences(hits_df):
    """
    Extracts the nucleotide sequence of each BLAST hit from its genome fasta file.
    
    Args
    ------
    hits_df (pd.DataFrame): Hits table with, at least, the columns
                            genome_dir, seq_id, start, end and strand.

    Return
    ------
    hits_df (pd.DataFrame): A copy of the input with an additional "sequence" column. 
    """

    hits_df = hits_df.copy()
    hits_df["sequence"] = ""

    # Group by genome so each .fna file is read only once
    for genome_dir, group in hits_df.groupby("genome_dir"):
        # Get the path to the .fna file in the genome directory
        fna_ext = glob.glob(os.path.join(genome_dir, "*.fna")) # list
        fna = fna_ext[0] # string

        # Parse the fasta file and create a dictionary of contig sequences
        contigs = SeqIO.to_dict(SeqIO.parse(fna, "fasta"))

        for idx, row in group.iterrows(): # Iterate over each hit in the group
            
            contig = contigs[row["seq_id"]] # For the current contig, get the SeqRecord object
            # Extract the sequence from the contig using the start and end coordinates
            seq = contig.seq[row["start"] - 1 : row["end"]] # (BLAST coordinates are 1-based and inclusive)
            # If the strand is negative, reverse complement the sequence
            if row["strand"] == "-":
                seq = seq.reverse_complement()

            # Add the extracted sequence to the hits_df
            hits_df.at[idx, "sequence"] = str(seq)

    return hits_df

def process_hits(tsv_path, dir_genomes, multifasta, 
                 min_length = 1000, extra_columns = {}):
    """
    Parses the BLAST results, removes the genomes without hits and extracts the
    sequences of the remaining hits.

    Args
    ------
    tsv_path (str): The path to the BLAST tabular output file.
    dir_genomes (str): The path to the directory containing one subdirectory per genome.
    multifasta (str): The path to the multifasta file used as BLAST database source.
    min_length (int): Minimum length of hits to be considered.
                      Default to 1000.
    extra_columns (dict, optional): Extra columns to add, as {name: path in the JSON}.
                                    Defaults to {}}.

    Return
    ------
    df (pd.DataFrame): One row per hit with the columns:
                       genome, genome_dir, seq_id, start, end, strand and sequence.
                       If extra_columns is provided, additional columns are added.

    """
    # Read the BLAST results and add strand and normalized coordinates
    hits_df = process_blast_tsv(tsv_path, min_length = min_length)
    # Map of each contig ID to its corresponding genome directory 
    contig2genome = contig_to_genome(dir_genomes)

    # Control to prevent deleting everything by mistake
    if hits_df.empty: # If BLAST returned no hits, something failed
        raise ValueError("BLAST returned no hits: nothing was deleted.")
    
    # Add genome_dir and genome columns to hits_df
    # For each seq_id in hits_df, map gets the corresponding genome directory from contig2genome
    hits_df["genome_dir"] = hits_df["seq_id"].map(contig2genome)
    # For each genome_dir in hits_df, os.path.basename gets the genome accession (filename)
    hits_df["genome"] = hits_df["genome_dir"].map(os.path.basename)

    # Genomes without any hit are removed from disk
    all_genomes = set(contig2genome.values()) # Set of all genome directories
    genomes_with_hit = set(hits_df["genome_dir"]) # Set of genome directories that have at least one hit
    
    to_delete = sorted(all_genomes - genomes_with_hit) # Set of genome directories without hits

    for each_del in to_delete: # Iterates over the genome directories to delete
        shutil.rmtree(each_del) # Remove the directory and all its contents
    # Removes multifasta file
    os.remove(multifasta)
    
    # Select and reorder the columns needed for the sequence extraction
    columns_to_keep = ["genome", "genome_dir", "seq_id", "start", "end", "strand"]
    hits_df = hits_df[columns_to_keep] # Keep only the relevant columns

    # Extracts the nucleotide sequence of each hit from its genome fasta file
    final_df = extract_hit_sequences(hits_df)

    # If extra_columns is provided, add the specified metadata columns to hits_df
    df = add_genome_metadata(final_df, extra_columns)

    return df

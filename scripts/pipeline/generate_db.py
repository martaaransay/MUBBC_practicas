# Import libraries 
import glob
import os
import shlex
import subprocess
import shutil

from pipeline.downloadgenomes import download_genomes
from pipeline.useful_funct import unique_genome_files


# --------Definition of functions
#------------------------ Function to download genomes for specific species
def download_species_genomes(species, outdir, flags, log_dir, quality_filter_flag = False ,target_num = None):
    """
    Uses the download_genomes function to download genomes for each species in the species dictionary.
    
    Args
    ------
    species (dict): A dictionary where keys are short names and values are long names of species.
    outdir (str): The output directory where genomes will be downloaded.
    flags (dict): Flags for the genome download.
    log_dir (str): Directory where logs will be saved.
    quality_filter_flag (bool, opt): If evaluates to True, executes quality_filter_genomes() function.
                                     Default is False.
    target_num (int, opt): The target number of genomes to download for each species.
                           Default to None.

    Returns
    -------
    dict: A dictionary where keys are short names of species and values are the directories of the genomes.
    """
    # Creates the output directory to avoid errors
    os.makedirs(outdir, exist_ok = True)
    dirs = {}
    # Iterates over the dictionary
    for short_name, long_name in species.items():
        print(f"-------------------- Downloading dehydrated genomes of: {long_name} --------------------")
        # Generates the path to save the genomes
        path = os.path.join(outdir, short_name, f"{short_name}.zip")
        # Executes download function
        dir = download_genomes(source = "taxon",
                               name = long_name, # Species to download
                               log_dir = log_dir,
                               quality_filter_flag = quality_filter_flag,
                               filename = path,
                               flags = flags,
                               target_num = target_num # number of genomes to download
            )

        dirs[short_name] = dir

    return dirs

#------------------------ Function to generate a multifasta file
def generate_multifasta(path_to_genomes, outfile, flag_deduplicate = True):
    """
    Generates a multifasta file with the concatenation of the fasta files of a specific path.

    Args
    ------
    path_to_genomes (str): The path to the genomes to be concatenated.
    outfile (str): The output file where the multifasta will be saved.
    flag_deduplicate (bool): If True, keeps only one file per assembly.
                             Defaults to True.
    """
    # Process the path to genomes to evaluate the wildcard correctly
    print("Generating multifasta...\n")
    files_to_db = sorted(glob.glob(path_to_genomes))

    if flag_deduplicate:
        dict_unique_files = unique_genome_files(files_to_db)
        files_to_db = list(dict_unique_files.values())
    
    # Writes the multifasta file
    with open(outfile, "w") as out:
        for file in files_to_db: # For each file
            with open(file, "r") as each_file:
                out.write(each_file.read()) # Writes each fasta in the multifasta
                out.write("\n") # To avoid errors adds a new line between each fasta
    return


#------------------------ Function to generate the database from the multifasta file
def generate_db(multifasta, log_dir, out_db):
    """
    Generates a BLAST database from a multifasta file.

    Args
    ------
    multifasta (str): Path to the multifasta file.
    log_dir (str): Path to the log directory.
    out_db (str): Path where the database will be saved.
    """
    # Path to a log file (to store errors)
    log_file = os.path.join(log_dir, "makeblastdb_errors.log")
    # The output directory is created to avoid errores
    os.makedirs(out_db, exist_ok=True)
    # Command to execute makeblastdb
    makeblastdb_cmd = ["makeblastdb", "-in", multifasta, 
                       "-dbtype", "nucl", # nucleotides database
                       "-out", out_db, 
                       "-parse_seqids"] # To include the sequence identifiers in the database
    print("Generating BLAST database...\n")
    print(f"Executing: \n{shlex.join(makeblastdb_cmd)}\n")
    # Execute the command
    try:
        subprocess.run(makeblastdb_cmd, 
                       check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing makeblastdb.\nDetails in {log_file}.")
        # Error stored into the log file
        with open(log_file, "a") as file:
            file.write("-------------- makeblastdb log --------------\n")
            file.write(f"{e.stderr}\n")
        return
    
    return

#------------------------ Function to remove all the files related to the database
def clean_db(outdir):
    """
    Removes all the files related to the database, including multifasta files, genomes, and the database itself.

    Args
    ------
    outdir (str): The output directory where the files to be removed are located.
    """
    shutil.rmtree(outdir) # Removes the entire directory and its contents

    return

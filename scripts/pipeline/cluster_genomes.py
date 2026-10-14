# Import libraries 
import subprocess
import shlex
import os
import csv
import glob

# --------Definition of functions
#------------------------ Function to run CheckM2 using the conda environment
def run_checkm2(genomes_paths, 
                out_dir,
                genome_ext = "fna", 
                lowmem = None, 
                threads = 4):
    """
    Executes CheckM2 predict to assess the quality of genomes using a conda environment.

    Args:
    ------
    genomes_paths (list): List of paths to the input genome files.
    out_dir (str): Path to the output directory where results will be saved.
    genome_ext (str, optional): File extension of the input genome files. 
                                Defaults to "fna".
    lowmem (bool, optional): If provided and evaluates to True, runs CheckM2 in low memory mode. 
                             Defaults to None.
    threads (int, optional): Number of threads to use for the CheckM2 execution. 
                             Defaults to 4.

    Return:
    ------
    checkm2_report (str or None): The complete path to the generated 'quality_report.tsv' file. 
                                  Returns None if an error occurs during execution.
    """
    #FIXME: poner en el readme o así del script la creación del entorno conda:
    # subprocess.run(["conda", "create", "-n", "checkm2", "-c", "conda-forge", "checkm2", "-y"], check=True)
    ## FIXME: Y TMB LA DATABASE DE DIAMOND
    # subprocess.run(["conda", "run", "-n", "checkm2", "checkm2", "database", "--download"], check=True)
    
    # Command to run CheckM2
    checkm2_cmd = ["checkm2", "predict", 
                   "--threads", str(threads),
                   "--input", *genomes_paths,
                   "-x", genome_ext, # Extension
                   "--output-directory", out_dir,
                   "--force", # to avoid errors and overwrite the files
                   "--remove_intermediates"] # To remove intermediate files

    if lowmem: # If an argument was provided, adds it to the flags list
        checkm2_cmd.append("--lowmem")

    # Command to run CheckM2 using the conda environment
    conda_checkm2_cmd = ["conda", "run", "-n", "checkm2", *checkm2_cmd]
 
    try:
        subprocess.run(conda_checkm2_cmd, capture_output=True, text=True, check=True)
    
    except subprocess.CalledProcessError as e:
        print(f"Error executing checkm2 predict:\n{e.stderr}")
        return None
    
    # Path to the checkm2 report
    checkm2_report = os.path.join(out_dir, "quality_report.tsv")

    return checkm2_report

#------------------------ Function to adapt the quality report to use dRep (genome_info)
def checkm2_to_genomeinfo(quality_report_tsv, genome_info_csv, genome_ext = "fna"):
    """
    Converts a CheckM2 quality report TSV file into a dRep-compatible genome_info CSV file.

    Args:
    ------
    quality_report_tsv (str): Path to the input TSV file containing the CheckM2 quality report.
    genome_info_csv (str): Path where the output genome info CSV file will be saved.
    genome_ext (str, optional): File extension to append to the genome names. 
                                Defaults to "fna".

    Return:
    ------
    genome_info_csv (str): The path to the generated genome info CSV file.
    """
    # Opens both files to read and write (tsv and csv)
    with open(quality_report_tsv) as tsv_in, \
         open(genome_info_csv, "w", newline = "") as csv_out: # newline to avoid \n

        tsv_file = csv.DictReader(tsv_in, delimiter = "\t") # Open file in reader mode

        csv_file = csv.writer(csv_out) # Open file in writer mode

        csv_file.writerow(["genome", "completeness", "contamination"]) # Header file

        for line in tsv_file: # For each line, gets the needed info
            genome_name = line["Name"] # name
            complete_genome_name = genome_name + "." + genome_ext # Adds the extension
            completeness = line["Completeness"] # completeness
            contamination = line["Contamination"] # contamination
            # Writes the info into the new csv file
            csv_file.writerow([complete_genome_name, 
                               completeness, 
                               contamination ])
    
    return genome_info_csv

#------------------------ Function to run dRep 
def cluster_genomes(genomes_paths,
                    adapted_file, 
                    threshold = 0.99,
                    threads = 1,
                    length = 50000,
                    completeness = 75,
                    contamination = 25,
                    skani_extra = None,
                    cov_thresh = 0.1,
                    chunksize = 500):
    """
    Executes genome clustering and dereplication using dRep (skani algorithm).

    Args:
    ------
    genomes_paths (list): List of paths to the input genome files.
    adapted_file (str): Reference file path used to determine the output directory structure.
    threshold (float, optional): Average Nucleotide Identity (ANI) threshold for secondary clustering. 
                                 Defaults to 0.99.
    threads (int, optional): Number of threads to use for CheckM2 and dRep. 
                             Defaults to 4.
    length (int, optional): Minimum genome length in base pairs. 
                            Defaults to 50000.
    completeness (float/int, optional): Minimum genome completeness percentage. 
                                        Defaults to 75.
    contamination (float/int, optional): Maximum genome contamination percentage. 
                                         Defaults to 25.
    skani_extra (str or list, optional): Additional arguments to pass to the skani algorithm. 
                                         Defaults to None.
    cov_thresh (float, optional): Minimum alignment coverage threshold. 
                                  Defaults to 0.1.
    chunksize (int, optional): Chunk size for multiround primary clustering. 
                               Defaults to 500.

    Return:
    ------
    genomes_path (str or None): A pattern path to the clustered genome files. 
                                Returns None if an error occurs during CheckM2 or dRep execution.
    """
    # Parental directory
    parental_dir = adapted_file.replace("/ncbi_dataset/data", "") 
        
    # Name of the folder with the clustering threshold
    folder_name = f"clustering_{threshold}"
    
    # New path
    out_dir = os.path.join(parental_dir, folder_name)

    # Create the output directory to avoid errors
    os.makedirs(out_dir, exist_ok=True)
    
    # Set the output directory of CheckM2
    checkm2_out = os.path.join(out_dir, "checkm2")
    # Execute CheckM2
    quality_report = run_checkm2(genomes_paths, 
                                 checkm2_out, 
                                 threads = threads)
    # In case the quality report is not generated, exits
    if quality_report is None:
        print("CheckM2 failed, aborting clustering.")
        return

    # Adapts the quality report tsv file to the needed csv file
    csv_file = os.path.join(out_dir, "genomeInfo.csv") # Path to the new csv
    genomeinfo_csv = checkm2_to_genomeinfo(quality_report,
                                           csv_file)

    # Text file with paths to genomes 
    # (instead of an expanded wildcard - recommended by dRep)
    genomes_list_file = os.path.join(out_dir, "genomes_list.txt")

    # Write the genome paths
    with open(genomes_list_file, "w") as f:
        for path in genomes_paths:
            f.write(path + "\n")

    # Command to execute dRep
    dRep_cmd = ["dRep", "dereplicate", "-p", str(threads),
                "-l", str(length), "-comp", str(completeness), "-con", str(contamination), 
                "--genomeInfo", genomeinfo_csv, "--S_algorithm", "skani", 
                "-sa", str(threshold), "-nc", str(cov_thresh), 
                "--low_ram_primary_clustering", "--multiround_primary_clustering", 
                "-g", genomes_list_file,  "--primary_chunksize", str(chunksize), 
                "--run_tertiary_clustering", "--skip_plots",  out_dir
                ]
    
    # In case more parameters for skani are needed
    if skani_extra: 
        if isinstance(skani_extra, str): # Argument as a string
            skani_extra = skani_extra.split()
        if isinstance(skani_extra, list): # Argument as a list
            dRep_cmd.extend(skani_extra)

    try:
        subprocess.run(dRep_cmd, capture_output=True, text=True, check=True)

    except subprocess.CalledProcessError as e:
        print(f"Error executing dRep dereplicate:\n{e.stderr}")
        return
    # Get the new paths to the clustered genome files
    genomes_path = os.path.join(out_dir, "dereplicated_genomes","*.fna")
    genomes_path_glob = sorted(glob.glob(genomes_path))
    return genomes_path_glob
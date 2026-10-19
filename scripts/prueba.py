import os
# Import functions from the pipeline module 
from pipeline.downloadgenomes import datasets_flags, download_genomes
from pipeline.generate_db import download_species_genomes, generate_multifasta, generate_db
from pipeline.useful_funct import contig_to_genome, list_fasta_files
from pipeline.parallel_download import process_all
from pipeline.blastn import run_blastn, process_hits
from pipeline.analyze_Pc import classify_pc_regex, info_to_csv
from pipeline.cluster_genomes import cluster_genomes

if __name__ == "__main__":

    # CHANGE THE PARAMETERS BEFORE THE EXECUTION:

    # OPTIONAL: number of genomes to download for each species
    n_genomes = 100
    # Directory where genomes will be stored
    data_dir = "../data_practicas/5genomes/Ecoli"
    # Directory where results will be saved
    results_dir = "results/PRUEBA" 
    # Directory where logs will be saved
    log_dir = f"{results_dir}/logs"

    # Flag to check if genomes are downloaded:
    #### If False, the download is included in the pipeline
    genomes_downloaded = True
    
    # Flag to check if the database is already created:
    #### If False, the generation of the database is included in the pipeline
    db_created = False

    # Flag to remove downloaded files (negative genomes, db and multifasta)
    #### If False, the files are not removed
    remove_files = False

    # Flag to filter genomes by quality metrics
    #### If True, executes quality_filter_genomes() function.
    quality_filter_flag = False

    # Flag to cluster genomes 
    #### If True, cluster genomes by
    cluster_genomes_flag = True

    #### ------------------------------------------------------- PIPELINE
    # genomes_dir = f"{data_dir}/prueba"
    # multifasta_path = f"{genomes_dir}/multifasta_genomes.fasta"

    # Control to avoid errors of non-existing directories
    # os.makedirs(genomes_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    if not genomes_downloaded: ### FIXME : CAMBIAR ESTO
        flags = datasets_flags(assembly_level="complete,chromosome,contig",
                               exclude_atypical=True,
                               mag=False)
    
        # species to download
        # names = {"Ecoli" : "Escherichia coli", 
        #         "Klebsiella" : "Klebsiella pneumoniae",
        #         "Pseudomonas" : "Pseudomonas aeruginosa",
        #         "Acinetobacter" : "Acinetobacter baumannii",
        #         "Enterobacter" : "Enterobacter"}
        
        ecoli = {"Ecoli" : "Escherichia coli",
                 "Klebsiella" : "Klebsiella pneumoniae"}
        ## For downloading normal genomes --> download_genomes() directly
        
        dirs = download_species_genomes(species = ecoli, 
                                        outdir = genomes_dir, 
                                        flags = flags,
                                        log_dir = log_dir,
                                        target_num = n_genomes,
                                        quality_filter_flag = quality_filter_flag)

        failed = process_all(dirs, log_dir, total_workers=30, n_split=500)
        
        if failed: 
            print(f"Failed to rehydrate: {failed}")
            exit(0)

    # path_to_genomes = f"{data_dir}/*/*/ncbi_dataset/data/*/*.fna"

    if not db_created:
        if cluster_genomes_flag:
            # List of the paths were the genomes are downloaded
            paths = list_fasta_files(directory = data_dir,
                                     depth = 4)


            clustered_path = cluster_genomes(genomes_paths = paths,
                                             output_dir = results_dir,
                                             log_dir = log_dir)
                            # threshold enables to change the similarity threshold to cluster
            # If no errrors occured, the path is updated
            if clustered_path:
                path_to_genomes = f"{clustered_path}/*.fna"

            
                

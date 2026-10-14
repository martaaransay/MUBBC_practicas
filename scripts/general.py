import os
# Import functions from the pipeline module 
from pipeline.downloadgenomes import datasets_flags, download_genomes
from pipeline.generate_db import download_species_genomes, generate_multifasta, generate_db
from pipeline.useful_funct import contig_to_genome
from pipeline.parallel_download import process_all
from pipeline.blastn import run_blastn, process_hits
from pipeline.analyze_Pc import classify_pc_regex, info_to_csv


if __name__ == "__main__":

    # CHANGE THE PARAMETERS BEFORE THE EXECUTION:

    # OPTIONAL: number of genomes to download for each species
    n_genomes = 100000 
    # Directory where genomes will be stored
    data_dir = "../data"
    # Directory where results will be saved
    results_dir = "results/Ecoli_complete" 
    # Directory where logs will be saved
    log_dir = f"{results_dir}/log"

    # Flag to check if the database is already created:
    #### If False, the generation of the database is included in the pipeline
    db_created = True

    # Flag to remove downloaded files (negative genomes, db and multifasta)
    #### If False, the files are not removed
    remove_files = False

    #### ------------------------------------------------------- PIPELINE
    genomes_dir = f"{data_dir}/Ecoli_complete"
    multifasta_path = f"{genomes_dir}/multifasta_genomes.fasta"

    # Control to avoid errors of non-existing directories
    os.makedirs(genomes_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    if not db_created:
        flags = datasets_flags(assembly_level="complete,chromosome,contig",
                               exclude_atypical=True,
                               mag=False)
    
        # species to download
        # names = {"Ecoli" : "Escherichia coli", 
        #         "Klebsiella" : "Klebsiella pneumoniae",
        #         "Pseudomonas" : "Pseudomonas aeruginosa",
        #         "Acinetobacter" : "Acinetobacter baumannii",
        #         "Enterobacter" : "Enterobacter"}
        
        ecoli = {"Ecoli" : "Escherichia coli"}
        ## For downloading normal genomes --> download_genomes() directly
        
        dirs = download_species_genomes(species = ecoli, 
                                        outdir = genomes_dir, 
                                        flags = flags,
                                        log_dir = log_dir,
                                        target_num = n_genomes)
        
        failed = process_all(dirs, log_dir, total_workers=30, n_split=500)
        
        if failed: 
            print(f"Failed to rehydrate: {failed}")
            exit(0)

        generate_multifasta(path_to_genomes = f"{genomes_dir}/*/*/ncbi_dataset/data/*/*.fna",
                            outfile = multifasta_path)
        generate_db(multifasta = multifasta_path,
                    log_dir = log_dir,
                    out_db = f"{genomes_dir}/db/db")

            
    query_file = "data/Int1/Int1_Pc_pattern.fa"

    extra_columns = {"bioproject": "assemblyInfo.bioprojectLineage.bioprojects.accession",
                     "biosample": "assemblyInfo.biosample.accession",
                     "organism_name": "organism.organismName",
                     "assembly_level": "assemblyInfo.assemblyLevel"}
    
    result, t = run_blastn(query = query_file, 
                           output = f"{results_dir}/blast/blast.tsv",
                           db = f"{genomes_dir}/db/db",
                           log_dir = log_dir,
                           max_target = len(contig_to_genome(genomes_dir))) # max num of seqs to keep: - mapping length

    hits_df = process_hits(result,
                           dir_genomes = genomes_dir,
                           multifasta = multifasta_path,
                           min_length = 1000,
                           extra_columns = extra_columns)

    for hit in hits_df.to_dict("records"):
        # Classify each sequence to a defined Pc variant
        pc_info = classify_pc_regex(hit["sequence"],
                                    hit_start = hit["start"],
                                    hit_strand = hit["strand"])
        
        # Add info for: genome, seq_id and species_folder
        row_info = {"genome": hit["genome"],   
                    "contig": hit["seq_id"],  
                    "species": hit["species_folder"]}

        row_info.update(pc_info) # Add the info about the classification of Pc

        
        for col in extra_columns.keys(): # Add extra cols
            row_info[col] = hit[col]

        info_to_csv(row_info, f"{results_dir}/final.csv")


